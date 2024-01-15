#!/usr/bin/env nextflow

nextflow.enable.dsl=2
SPLIT_BY_CHROM = true
MAF_THRESHOLD = params.maf
GENCODE_GTF = params.gencode_gtf
GENE_MAPPABILITY = params.gene_mappability
CROSS_MAPPABILITY = params.cross_mappability
SNP_MAPPABILITY = params.snp_mappability
SNP_MAPPABILITY_INDEX = params.snp_mappability_index
//TENSORQTL_CONTAINER = 'docker.io/porchard/tensorqtl_dev:774fb5b'
TENSORQTL_CONTAINER = 'docker.io/porchard/tensorqtl_dev:0124d74'
CLUSTER_OPTIONS = "--account=dcmb_dept --partition=dcmb_dept"
L = 10
CHROMS = (1..22) + ['X']


PHENOTYPE_PCS = [
    'Lung': '10',
    'Whole_blood': '10',
    'T_cell': '10',
    'PBMC': '10',
    'Monocyte': '10',
    'Nasal_epithelial': '10',
    'Lung___EUR': '10',
    'Monocyte___EUR': '5',
    'Nasal_epithelial___EUR': '10',
    'PBMC___AFR': '5',
    'PBMC___EAS': '5',
    'PBMC___EUR': '10',
    'T_cell___EUR': '5',
    'Whole_blood___AFR': '10',
    'Whole_blood___EUR': '10'
]

DROP_PHENOTYPE_PCS = [
    'Dummy': [25] + (51..100)
]


def covariate_file_to_tissue (f) {
        f.getName().tokenize('.')[0]
}

def covariate_file_to_pcs (f) {
        f.getName().tokenize('.')[2]
}

def phenotype_file_to_tissue (f) {
        f.getName().tokenize('.')[0]
}


process prep_covariates {

    tag "${tissue}"
    container 'docker.io/porchard/general:20220406125608'
    clusterOptions "${CLUSTER_OPTIONS}"
    memory '10 GB'

    input:
    tuple val(tissue), path('covariates.txt')

    output:
    tuple val(tissue), path("${tissue}.covariates.tsv")

    script:
    DROP_COVARIATES_CMD = DROP_PHENOTYPE_PCS.containsKey(tissue) ? ("grep -v -w ") + DROP_PHENOTYPE_PCS[tissue].collect({y -> "-e phenotype_PC" + y }).join(' ') + " covariates.txt > ${tissue}.covariates.tsv": "cp covariates.txt ${tissue}.covariates.tsv"

    """
    # drop phenotype_PCs from covariate matrix
    ${DROP_COVARIATES_CMD}
    """

}

process prep_phenotypes {

    tag "${tissue}"
    container 'docker.io/porchard/general:20220406125608'
    clusterOptions "${CLUSTER_OPTIONS}"
    memory '10 GB'

    input:
    tuple val(tissue), path('phenotypes.tmp.bed.gz')

    output:
    tuple val(tissue), path("${tissue}.phenotypes.bed.gz")

    """
    zcat phenotypes.tmp.bed.gz | perl -pe 's/^chr//' | gzip -c > ${tissue}.phenotypes.bed.gz
    """

}

process prep_genotypes_per_chrom {

    tag "${tissue} ${chrom}"
    container 'docker.io/porchard/general:20220406125608'
    clusterOptions "${CLUSTER_OPTIONS}"
    memory '20 GB'

    input:
    tuple val(chrom), path(genotypes), val(tissue), path("${tissue}.phenotypes.bed.gz")

    output:
    tuple val(tissue), val(chrom), path("${chrom}.genotypes.bed"), path("${chrom}.genotypes.bim"), path("${chrom}.genotypes.fam")

    """
    # filter genotypes to the individuals of interest
    zcat ${tissue}.phenotypes.bed.gz | awk 'NR==1' | cut -f5- | perl -pe 's/\\t/\\n/g' | perl -pe 's/^/0\\t/' > keep.txt
    plink -bfile $chrom --output-chr M --keep-allele-order --make-bed --maf $MAF_THRESHOLD --out ${chrom}.genotypes --keep keep.txt --indiv-sort file keep.txt
    """

}


process remove_snps_with_mappability_below_1 {

    tag "${tissue} ${chrom}"
    container 'docker.io/porchard/general:20220406125608'
    clusterOptions "${CLUSTER_OPTIONS}"
    memory '50 GB'

    input:
    tuple val(tissue), val(chrom), path('genotypes.bed'), path('genotypes.bim'), path('genotypes.fam'), path(snp_mappability), path(snp_mappability_index)

    output:
    tuple val(tissue), val(chrom), path("${chrom}.genotypes.bed"), path("${chrom}.genotypes.bim"), path("${chrom}.genotypes.fam")

    """
    cut -f2 genotypes.bim > variants.txt
    get-variant-mappability.py $snp_mappability variants.txt > variant-mappability.txt
    variant-mappability-to-mappable-variants.py variant-mappability.txt > mappable.txt
    plink -bfile genotypes --output-chr M --keep-allele-order --extract mappable.txt --make-bed --out ${chrom}.genotypes
    """

}


process remove_monomorphic {

    tag "${tissue} ${chrom}"
    container 'docker.io/porchard/general:20220406125608'
    clusterOptions "${CLUSTER_OPTIONS}"
    memory '50 GB'

    input:
    tuple val(tissue), val(chrom), path('genotypes.bed'), path('genotypes.bim'), path('genotypes.fam')

    output:
    tuple val(tissue), path("${chrom}.genotypes.bed"), path("${chrom}.genotypes.bim"), path("${chrom}.genotypes.fam")

    """
    list-monomorphic-fast.py genotypes.bed > drop-variants.txt
    plink -bfile genotypes --output-chr M --keep-allele-order --exclude drop-variants.txt --make-bed --out ${chrom}.genotypes
    """

}


process merge_genotypes_across_chrom {

    tag "${tissue}"
    container 'docker.io/porchard/general:20220406125608'
    clusterOptions "${CLUSTER_OPTIONS}"
    memory '60 GB'

    input:
    tuple val(tissue), path(bed), path(bim), path(fam)

    output:
    tuple val(tissue), path("mergedchroms*")

    """
    for c in {1..22} X; do echo chr\${c}.genotypes >> merge.txt; done
    plink --merge-list merge.txt --make-bed --out mergedchroms --output-chr M --keep-allele-order --indiv-sort 'none'
    """

}

process run_trans_permutations {

    tag "${tissue}"
    publishDir "${params.results}/trans"
    clusterOptions "${CLUSTER_OPTIONS} --gres=gpu:1"
    container "${TENSORQTL_CONTAINER}"
    beforeScript "export PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:100; python -m torch.utils.collect_env"
    memory '120 GB'
    time '10h'

    input:
    tuple val(tissue), path(phenotypes), path(covariates), path(genotypes)

    output:
    tuple val(tissue), path("${tissue}.permutations.pickle"), emit: permutations_pickle

    script:
    prefix = "${tissue}"

    """
    python3.8 -m tensorqtl --mode trans_permutations --invnorm --permutations 20000 --seed 2022 --output_text --maf_threshold $MAF_THRESHOLD --batch_size 500 --return_r2 --covariates $covariates mergedchroms $phenotypes $prefix
    """

}

// if have enough memory, can run_trans and then straight to add_mappability_info
// otherwise, use run_trans_by_chrom --> merge_trans_across_chrom --> add_mappability_info
process run_trans {

    tag "${tissue}"
    publishDir "${params.results}/trans"
    clusterOptions "${CLUSTER_OPTIONS} --gres=gpu:1"
    container "${TENSORQTL_CONTAINER}"
    beforeScript "export PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:100; python -m torch.utils.collect_env"
    memory '120 GB'
    time '10h'

    input:
    tuple val(tissue), path(phenotypes), path(covariates), path(genotypes)

    output:
    tuple val(tissue), path("${tissue}.trans_qtl_pairs.txt.gz"), emit: pairs_df

    script:
    prefix = "${tissue}"

    """
    python3.8 -m tensorqtl --mode trans --invnorm --seed 2022 --output_text --maf_threshold $MAF_THRESHOLD --batch_size 2000 --return_r2 --covariates $covariates mergedchroms $phenotypes $prefix
    """

}


process run_trans_by_chrom {

    tag "${tissue} ${chrom}"
    publishDir "${params.results}/trans-per-chrom"
    clusterOptions "${CLUSTER_OPTIONS} --gres=gpu:1"
    container "${TENSORQTL_CONTAINER}"
    beforeScript "export PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:100; python -m torch.utils.collect_env"
    memory '120 GB'
    time '10h'

    input:
    tuple val(tissue), path(phenotypes), path(covariates), path(genotypes)
    each chrom

    output:
    tuple val(tissue), path("${tissue}.${chrom}.phenotypes.bed.gz")
    tuple val(tissue), path("${tissue}.${chrom}.trans_qtl_pairs.txt.gz"), emit: to_merge

    script:
    prefix = "${tissue}.${chrom}"

    """
    zcat $phenotypes | awk '(\$1=="$chrom") || (\$1=="#chr")' | gzip -c > ${tissue}.${chrom}.phenotypes.bed.gz
    python3.8 -m tensorqtl --mode trans --invnorm --seed 2022 --output_text --maf_threshold $MAF_THRESHOLD --batch_size 2000 --return_r2 --covariates $covariates mergedchroms ${tissue}.${chrom}.phenotypes.bed.gz $prefix
    """

}

process merge_trans_across_chrom {

    tag "${tissue}"
    publishDir "${params.results}/trans"
    clusterOptions "${CLUSTER_OPTIONS}"
    container "${TENSORQTL_CONTAINER}"
    memory '10 GB'
    time '1h'

    input:
    tuple val(tissue), path(pairs)

    output:
    tuple val(tissue), path("${tissue}.trans_qtl_pairs.txt.gz"), emit: pairs_df

    """
    zcat ${pairs[0]} | awk 'NR==1' > header.txt
    cp header.txt ${tissue}.trans_qtl_pairs.txt
    zcat ${pairs.join(' ')} | grep -v -f header.txt >> ${tissue}.trans_qtl_pairs.txt
    gzip ${tissue}.trans_qtl_pairs.txt
    """

}


process add_mappability_info {

    tag "${tissue}"
    publishDir "${params.results}/trans-with-mappability-info"
    container 'docker.io/porchard/general:20220406125608'
    clusterOptions "${CLUSTER_OPTIONS}"
    memory '90 GB'
    time '7h'

    input:
    tuple val(tissue), path(pairs), path(gtf), path(cross_mappability), path(gene_mappability)

    output:
    tuple val(tissue), path("${tissue}.trans_qtl_pairs.with_mappability.txt.gz")

    """
    add-mappability-info-to-trans-qtl.py --trans $pairs --gtf $gtf --cross-mappability $cross_mappability --cross-mappability-window 1000000 --gene-mappability $gene_mappability | gzip -c > ${tissue}.trans_qtl_pairs.with_mappability.txt.gz
    """

}


process apply_permutations {

    tag "${tissue}"
    publishDir "${params.results}/trans-top"
    container "${TENSORQTL_CONTAINER}"
    clusterOptions "${CLUSTER_OPTIONS}"
    memory '20 GB'

    input:
    tuple val(tissue), path(pairs), path(permutations_pickle), path(phenotype_groups)

    output:
    tuple val(tissue), path("${tissue}.trans_qtl.top.txt")

    """
    trans-sqtl-filter-apply-permutations-fdr-correct.py --phenotype-groups $phenotype_groups --pairs $pairs --permutations $permutations_pickle > ${tissue}.trans_qtl.top.txt
    """

}

process prep_finemapping {

    tag "${tissue}"
    publishDir "${params.results}/trans-susie-phenotypes"
    container 'docker.io/porchard/general:20220406125608'
    clusterOptions "${CLUSTER_OPTIONS}"
    memory '30 GB'

    input:
    tuple val(tissue), path(top), path(phenotypes)

    output:
    tuple val(tissue), path("${tissue}.phenotypes_for_finemapping.bed.gz"), optional: true

    """
    prep-trans-eqtl-finemapping.py $top $phenotypes ${tissue}.phenotypes_for_finemapping.bed.gz
    """

}

// could reduce L to 5?
process scan_susie {

    tag "${tissue}"
    publishDir "${params.results}/trans-susie"
    clusterOptions "${CLUSTER_OPTIONS} --gres=gpu:1"
    container "${TENSORQTL_CONTAINER}"
    beforeScript "export PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:100; python -m torch.utils.collect_env"
    memory '120 GB'

    input:
    tuple val(tissue), path("phenotypes.tmp.bed.gz"), path(genotypes), path(covariates)

    output:
    tuple val(tissue), path("${tissue}.susie.pickle"), emit: pickle
    path("${tissue}.SuSiE_summary.parquet")

    script:
    prefix = "${tissue}"

    """
    zcat phenotypes.tmp.bed.gz | perl -pe 's/^chr//' | gzip -c > ${tissue}.phenotypes.bed.gz
    # make fake cis-output file to force fine-mapping of everything
    printf "phenotype_id\\tqval\\n" > cis-output.txt
    zcat ${tissue}.phenotypes.bed.gz | awk 'NR>1' | cut -f4 | perl -pe 's/\\n/\\t0.01\\n/' >> cis-output.txt
    python3.8 -m tensorqtl --seed 2021 --mode cis_susie --invnorm --cis_output cis-output.txt --maf_threshold $MAF_THRESHOLD --L $L --covariates $covariates mergedchroms ${tissue}.phenotypes.bed.gz $prefix
    mv ${tissue}.SuSiE.pickle ${tissue}.susie.pickle
    """

}


process susie_pickle_to_summary {

    tag "${tissue}"
    publishDir "${params.results}/trans-susie"
    container 'docker.io/porchard/general:20220406125608'
    clusterOptions "${CLUSTER_OPTIONS}"
    memory '30 GB'

    input:
    tuple val(tissue), path(pickle)

    output:
    tuple val(tissue), path("${tissue}.cs.txt"), path("${tissue}.converged.txt")

    """
    susie-pickle-2 summarize --pickle $pickle --single-effect-pips --deduplicate --prefix ${tissue}.
    """

}

workflow trans {
    COVARIATE_FILES = Channel.fromPath(params.covariate_glob).map({it -> [covariate_file_to_tissue(it), covariate_file_to_pcs(it), it]}).filter({it -> it[1] == PHENOTYPE_PCS[it[0]]}).map({it -> [it[0], it[2]]})
    GENOTYPE_FILES = Channel.fromPath(params.plink_glob).map({it -> [it.getName().tokenize('.')[0], it]}).groupTuple()
    PHENOTYPE_FILES = Channel.fromPath(params.phenotype_glob).map({it -> [phenotype_file_to_tissue(it), it]})
    phenotype_groups = Channel.fromPath(params.phenotype_groups_glob).map({it -> [phenotype_file_to_tissue(it), it]})
    snp_mappability = Channel.fromPath(SNP_MAPPABILITY)
    snp_mappability_index = Channel.fromPath(SNP_MAPPABILITY_INDEX)
    cross_mappability = Channel.fromPath(CROSS_MAPPABILITY)
    gene_mappability = Channel.fromPath(GENE_MAPPABILITY)
    gencode_gtf = Channel.fromPath(GENCODE_GTF)

    covariates_prepped = prep_covariates(COVARIATE_FILES)
    phenotypes_prepped = prep_phenotypes(PHENOTYPE_FILES)
    filtered_genotypes = prep_genotypes_per_chrom(GENOTYPE_FILES.combine(phenotypes_prepped)).combine(snp_mappability).combine(snp_mappability_index) | remove_snps_with_mappability_below_1 | remove_monomorphic

    genotypes_merged = filtered_genotypes.groupTuple() | merge_genotypes_across_chrom // tissue, plink_files

    trans_input = phenotypes_prepped.combine(covariates_prepped, by: 0).combine(genotypes_merged, by: 0)
    trans_permutations = run_trans_permutations(trans_input)
    
    if (SPLIT_BY_CHROM) {
        trans_out = run_trans_by_chrom(trans_input, CHROMS).to_merge.groupTuple() | merge_trans_across_chrom
    } else {
        trans_out = run_trans(trans_input)
    }
    
    trans_out_with_mappability = trans_out.pairs_df.combine(gencode_gtf).combine(cross_mappability).combine(gene_mappability) | add_mappability_info
    top = trans_out_with_mappability.combine(trans_permutations.permutations_pickle, by: 0).combine(phenotype_groups, by: 0) | apply_permutations

    finemapping_phenotypes = top.combine(phenotypes_prepped, by: 0) | prep_finemapping
    finemapping_genotypes = genotypes_merged
    susie_in = finemapping_phenotypes.combine(finemapping_genotypes, by: 0).combine(covariates_prepped, by: 0)
    susie_out = scan_susie(susie_in).pickle | susie_pickle_to_summary

}
