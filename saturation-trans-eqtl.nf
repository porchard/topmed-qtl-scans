#!/usr/bin/env nextflow

nextflow.enable.dsl=2
MAF_THRESHOLD = params.maf
MAF_THRESHOLD_FINEMAPPING = params.maf_finemapping
GENCODE_GTF = params.gencode_gtf
GENE_MAPPABILITY = params.gene_mappability
CROSS_MAPPABILITY = params.cross_mappability
SNP_MAPPABILITY = params.snp_mappability
SNP_MAPPABILITY_INDEX = params.snp_mappability_index
//TENSORQTL_CONTAINER = 'docker.io/porchard/tensorqtl_dev:774fb5b'
TENSORQTL_CONTAINER = 'docker.io/porchard/tensorqtl_dev:0124d74'
CLUSTER_OPTIONS = "--account=dcmb_dept --partition=dcmb_dept"
L = 10


DROP_PHENOTYPE_PCS = [
    'Whole_blood_500': (51..100),
    'Whole_blood_1500': [25] + (51..100),
    'Whole_blood_2000': [24] + (51..100),
    'Whole_blood_2500': [24] + (51..100),
    'Whole_blood_3000': [24, 23] + (51..100),
    'Whole_blood_3500': [23, 24, 25] + (51..100),
    'Whole_blood_4000': [24, 23] + (51..100),
    'Whole_blood_4500': [25, 24, 23] + (51..100),
    'Whole_blood_5000': [25] + (51..100),
    'Whole_blood_5500': [25, 24] + (51..100),
    'Whole_blood_6000': [25] + (51..100),
    'Whole_blood_6500': [25] + (51..100)
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

// first, prep covariates
// next, prep each chrom independently
// then merge chroms
// then scan

process prep_covariates {

    tag "${tissue}"
    container 'docker.io/porchard/general:20220406125608'
    clusterOptions "${CLUSTER_OPTIONS}"
    memory '10 GB'

    input:
    tuple val(tissue), val(phenotype_pcs), path('covariates.txt')

    output:
    tuple val(tissue), val(phenotype_pcs), path("${tissue}.${phenotype_pcs}.covariates.tsv")

    script:
    DROP_COVARIATES_CMD = DROP_PHENOTYPE_PCS.containsKey(tissue) ? ("grep -v -w ") + DROP_PHENOTYPE_PCS[tissue].collect({y -> "-e phenotype_PC" + y }).join(' ') + " covariates.txt > ${tissue}.${phenotype_pcs}.covariates.tsv": "cp covariates.txt ${tissue}.${phenotype_pcs}.covariates.tsv"

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
    tuple val(tissue), path(phenotypes), val(phenotype_pcs), path(covariates), path(genotypes)

    output:
    tuple val(tissue), val(phenotype_pcs), path("${tissue}.${phenotype_pcs}.permutations.pickle"), emit: permutations_pickle

    script:
    prefix = "${tissue}.${phenotype_pcs}"

    """
    python3.8 -m tensorqtl --mode trans_permutations --invnorm --permutations 20000 --seed 2022 --output_text --maf_threshold $MAF_THRESHOLD --batch_size 500 --return_r2 --covariates $covariates mergedchroms $phenotypes $prefix
    """

}


process run_trans {

    tag "${tissue}"
    publishDir "${params.results}/trans"
    clusterOptions "${CLUSTER_OPTIONS} --gres=gpu:1"
    container "${TENSORQTL_CONTAINER}"
    beforeScript "export PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:100; python -m torch.utils.collect_env"
    memory '120 GB'
    time '10h'

    input:
    tuple val(tissue), path(phenotypes), val(phenotype_pcs), path(covariates), path(genotypes)

    output:
    tuple val(tissue), val(phenotype_pcs), path("${tissue}.${phenotype_pcs}.trans_qtl_pairs.txt.gz"), emit: pairs_df

    script:
    prefix = "${tissue}.${phenotype_pcs}"

    """
    python3.8 -m tensorqtl --mode trans --invnorm --seed 2022 --output_text --maf_threshold $MAF_THRESHOLD --batch_size 500 --return_r2 --covariates $covariates mergedchroms $phenotypes $prefix
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
    tuple val(tissue), val(phenotype_pcs), path(pairs), path(gtf), path(cross_mappability), path(gene_mappability)

    output:
    tuple val(tissue), val(phenotype_pcs), path("${tissue}.${phenotype_pcs}.trans_qtl_pairs.with_mappability.txt.gz")

    """
    add-mappability-info-to-trans-qtl.py --trans $pairs --gtf $gtf --cross-mappability $cross_mappability --cross-mappability-window 1000000 --gene-mappability $gene_mappability | gzip -c > ${tissue}.${phenotype_pcs}.trans_qtl_pairs.with_mappability.txt.gz
    """

}


process apply_permutations {

    tag "${tissue}"
    publishDir "${params.results}/trans-top"
    container "${TENSORQTL_CONTAINER}"
    clusterOptions "${CLUSTER_OPTIONS}"
    memory '20 GB'

    input:
    tuple val(tissue), val(phenotype_pcs), path(pairs), path(permutations_pickle)

    output:
    tuple val(tissue), val(phenotype_pcs), path("${tissue}.${phenotype_pcs}.trans_qtl.top.txt")

    """
    trans-eqtl-filter-apply-permutations-fdr-correct.py --pairs $pairs --permutations $permutations_pickle > ${tissue}.${phenotype_pcs}.trans_qtl.top.txt
    """

}

// process prep_finemapping {

//     tag "${tissue}"
//     publishDir "${params.results}/trans-susie-phenotypes"
//     container 'docker.io/porchard/general:20220406125608'
//     clusterOptions "${CLUSTER_OPTIONS}"
//     memory '30 GB'

//     input:
//     tuple val(tissue), path(top), path(phenotypes)

//     output:
//     tuple val(tissue), path("${tissue}.phenotypes_for_finemapping.bed.gz"), optional: true

//     """
//     prep-trans-eqtl-finemapping.py $top $phenotypes ${tissue}.phenotypes_for_finemapping.bed.gz
//     """

// }


// process prep_genotypes_per_chrom_for_finemapping {

//     tag "${tissue} ${chrom}"
//     container 'docker.io/porchard/general:20220406125608'
//     clusterOptions "${CLUSTER_OPTIONS}"
//     memory '20 GB'

//     input:
//     tuple val(chrom), path(genotypes), val(tissue), path("${tissue}.phenotypes.bed.gz")

//     output:
//     tuple val(tissue), path("${chrom}.genotypes.bed"), path("${chrom}.genotypes.bim"), path("${chrom}.genotypes.fam")

//     """
//     # filter genotypes to the individuals of interest
//     zcat ${tissue}.phenotypes.bed.gz | awk 'NR==1' | cut -f5- | perl -pe 's/\\t/\\n/g' | perl -pe 's/^/0\\t/' > keep.txt
//     plink -bfile $chrom --output-chr M --keep-allele-order --make-bed --maf $MAF_THRESHOLD_FINEMAPPING --out ${chrom}.genotypes --keep keep.txt --indiv-sort file keep.txt
//     """

// }


// process merge_genotypes_across_chrom_for_finemapping {

//     tag "${tissue}"
//     container 'docker.io/porchard/general:20220406125608'
//     clusterOptions "${CLUSTER_OPTIONS}"
//     memory '60 GB'

//     input:
//     tuple val(tissue), path(bed), path(bim), path(fam)

//     output:
//     tuple val(tissue), path("mergedchroms*")

//     """
//     for c in {1..22} X; do echo chr\${c}.genotypes >> merge.txt; done
//     plink --merge-list merge.txt --make-bed --out mergedchroms --output-chr M --keep-allele-order --indiv-sort 'none'
//     """

// }


// // could reduce L to 5?
// process scan_susie {

//     tag "${tissue}"
//     publishDir "${params.results}/trans-susie"
//     clusterOptions "${CLUSTER_OPTIONS} --gres=gpu:1"
//     container "${TENSORQTL_CONTAINER}"
//     beforeScript "export PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:100; python -m torch.utils.collect_env"
//     memory '120 GB'

//     input:
//     tuple val(tissue), path("phenotypes.tmp.bed.gz"), path(genotypes), path(covariates)

//     output:
//     tuple val(tissue), path("${tissue}.susie.pickle"), emit: pickle
//     path("${tissue}.SuSiE_summary.parquet")

//     script:
//     prefix = "${tissue}"

//     """
//     zcat phenotypes.tmp.bed.gz | perl -pe 's/^chr//' | gzip -c > ${tissue}.phenotypes.bed.gz
//     # make fake cis-output file to force fine-mapping of everything
//     printf "phenotype_id\\tqval\\n" > cis-output.txt
//     zcat ${tissue}.phenotypes.bed.gz | awk 'NR>1' | cut -f4 | perl -pe 's/\\n/\\t0.01\\n/' >> cis-output.txt
//     python3.8 -m tensorqtl --seed 2021 --mode cis_susie --invnorm --cis_output cis-output.txt --maf_threshold $MAF_THRESHOLD_FINEMAPPING --L $L --covariates $covariates mergedchroms ${tissue}.phenotypes.bed.gz $prefix
//     mv ${tissue}.SuSiE.pickle ${tissue}.susie.pickle
//     """

// }


// process susie_pickle_to_summary {

//     tag "${tissue}"
//     publishDir "${params.results}/trans-susie"
//     container 'docker.io/porchard/general:20220406125608'
//     clusterOptions "${CLUSTER_OPTIONS}"
//     memory '30 GB'

//     input:
//     tuple val(tissue), path(pickle)

//     output:
//     tuple val(tissue), path("${tissue}.cs.txt"), path("${tissue}.converged.txt")

//     """
//     susie-pickle-2 summarize --pickle $pickle --single-effect-pips --deduplicate --prefix ${tissue}.
//     """

// }


workflow trans {
    COVARIATE_FILES = Channel.fromPath(params.covariate_glob).map({it -> [covariate_file_to_tissue(it), covariate_file_to_pcs(it), it]})
    GENOTYPE_FILES = Channel.fromPath(params.plink_glob).map({it -> [it.getName().tokenize('.')[0], it]}).groupTuple()
    PHENOTYPE_FILES = Channel.fromPath(params.phenotype_glob).map({it -> [phenotype_file_to_tissue(it), it]})
    snp_mappability = Channel.fromPath(SNP_MAPPABILITY)
    snp_mappability_index = Channel.fromPath(SNP_MAPPABILITY_INDEX)
    cross_mappability = Channel.fromPath(CROSS_MAPPABILITY)
    gene_mappability = Channel.fromPath(GENE_MAPPABILITY)
    gencode_gtf = Channel.fromPath(GENCODE_GTF)

    covariates_prepped = prep_covariates(COVARIATE_FILES)
    phenotypes_prepped = prep_phenotypes(PHENOTYPE_FILES)
    filtered_genotypes = prep_genotypes_per_chrom(GENOTYPE_FILES.combine(phenotypes_prepped)).combine(snp_mappability).combine(snp_mappability_index) | remove_snps_with_mappability_below_1 | remove_monomorphic

    genotypes_merged = filtered_genotypes.groupTuple() | merge_genotypes_across_chrom // tissue, plink_files
    
    trans_in = phenotypes_prepped.combine(covariates_prepped, by: 0).combine(genotypes_merged, by: 0)
    trans_out = run_trans(trans_in)
    trans_permutations = run_trans_permutations(trans_in)
    trans_out_with_mappability = trans_out.pairs_df.combine(gencode_gtf).combine(cross_mappability).combine(gene_mappability) | add_mappability_info
    top = trans_out_with_mappability.combine(trans_permutations.permutations_pickle, by: [0, 1]) | apply_permutations







    // finemapping_phenotypes = top.combine(phenotypes_prepped, by: 0) | prep_finemapping
    // finemapping_genotypes = genotypes_merged
    // //finemapping_genotypes = prep_genotypes_per_chrom(GENOTYPE_FILES.combine(phenotypes_prepped)).groupTuple() | merge_genotypes_across_chrom // tissue, plink_files
    // susie_in = finemapping_phenotypes.combine(finemapping_genotypes, by: 0).combine(covariates_prepped, by: 0)
    // susie_out = scan_susie(susie_in).pickle | susie_pickle_to_summary

}
