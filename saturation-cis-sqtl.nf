#!/usr/bin/env nextflow

nextflow.enable.dsl=2
MAF_THRESHOLD = params.maf
//L = params.L
LS = [10]
PERMUTATION_GLOB = params.permutation_glob // {tissue}.{chrom}.{phenotype_pcs}.*
TENSORQTL_CONTAINER = 'docker.io/porchard/tensorqtl_dev:5ea048f_pyarrow_11'
CLUSTER_OPTIONS = "--account=dcmb_dept --partition=dcmb_dept"
COVARIATE_GLOB = params.covariate_glob // {tissue}.tensorqtl-in.{pcs}.covariates.tsv
PHENOTYPE_GLOB = params.phenotype_glob // {tissue}.{chrom}.phenotypes.bed.gz
PLINK_GLOB = params.plink_glob // {tissue}.{chrom}.{bim,bam,fam}
PHENOTYPE_GROUP_GLOB = params.phenotype_group_glob // {tissue}.{chrom}.phenotype-groups.txt

COVARIATE_FILES = Channel.fromPath(COVARIATE_GLOB).map({it -> [covariate_file_to_tissue(it), covariate_file_to_pcs(it), it]}) // tissue, phenotype PCs, covariate file
PHENOTYPE_FILES = Channel.fromPath(PHENOTYPE_GLOB).map({it -> it.getName().tokenize('.')[0..1] + [it]}) // tissue, chrom, phenotype file
GENOTYPE_FILES = Channel.fromPath(PLINK_GLOB).map({it -> it.getName().tokenize('.')[0..1] + [it]}).groupTuple(by: [0, 1]) // tissue, chrom, plink files

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

SUSIE_MAX_L = [
    'Lung': 10,
    'Whole_blood': 20,
    'T_cell': 10,
    'PBMC': 10,
    'Monocyte': 10,
    'Nasal_epithelial': 10,
    'Lung___EUR': 10,
    'Monocyte___EUR': 10,
    'Nasal_epithelial___EUR': 10,
    'PBMC___AFR': 10,
    'PBMC___EAS': 10,
    'PBMC___EUR': 10,
    'T_cell___EUR': 10,
    'Whole_blood___AFR': 10,
    'Whole_blood___EUR': 10
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


process scan_permutations {

    tag "${tissue} ${chrom} ${phenotype_pcs}"
    publishDir "${params.results}/cis/permutations"
    clusterOptions "${CLUSTER_OPTIONS} --gres=gpu:1"
    container "${TENSORQTL_CONTAINER}"
    beforeScript "export PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:100"
    memory '30 GB'
    time '7h'

    input:
    tuple val(tissue), val(chrom), path(phenotypes), path(phenotype_groups), path(genotypes), val(phenotype_pcs), path(covariates)

    output:
    tuple val(tissue), val(phenotype_pcs), path("${tissue}.${chrom}.${phenotype_pcs}.cis_qtl.txt.gz")

    when:
    ['0', '5', '10', '15', '20', '25', '30'].contains(phenotype_pcs)

    script:
    prefix = "${tissue}.${chrom}.${phenotype_pcs}"

    """
    python3.8 -m tensorqtl --qvalue_lambda 0 --seed 2021 --mode cis --invnorm --maf_threshold $MAF_THRESHOLD --phenotype_groups $phenotype_groups --covariates $covariates ${tissue}.${chrom} $phenotypes $prefix
    """

}



process scan_nominal {

    tag "${tissue} ${chrom} ${phenotype_pcs}"
    publishDir "${params.results}/cis/nominal"
    clusterOptions "${CLUSTER_OPTIONS} --gres=gpu:1"
    container "${TENSORQTL_CONTAINER}"
    beforeScript "export PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:100; python -m torch.utils.collect_env"
    memory '70 GB'
    time '4h'

    input:
    tuple val(tissue), val(chrom), path(phenotypes), path(phenotype_groups), path(genotypes), val(phenotype_pcs), path(covariates), path(permutation_results)

    output:
    tuple val(tissue), val(chrom), path("${tissue}.${chrom}.${phenotype_pcs}.cis_qtl_pairs.*.parquet"), emit: all_pairs
    tuple val(tissue), val(chrom), path("${tissue}.${chrom}.${phenotype_pcs}.cis_qtl.signif_pairs.parquet"), emit: signif_pairs

    when:
    PHENOTYPE_PCS[tissue] == phenotype_pcs

    script:
    prefix = "${tissue}.${chrom}.${phenotype_pcs}"

    """
    python3.8 -m tensorqtl --seed 2021 --mode cis_nominal --invnorm --cis_output $permutation_results --phenotype_groups $phenotype_groups --maf_threshold $MAF_THRESHOLD --covariates $covariates ${tissue}.${chrom} $phenotypes $prefix
    """

}

process merge_significant_pairs {

    tag "${tissue} ${phenotype_pcs}"
    publishDir "${params.results}/cis/nominal"
    clusterOptions "${CLUSTER_OPTIONS}"
    //container 'library://porchard/default/general:20220107'
    container "${TENSORQTL_CONTAINER}"
    memory '50 GB'
    time '2h'

    input:
    tuple val(tissue), val(chroms), path(signif_pairs)

    output:
    path("${tissue}.cis_qtl.signif_pairs.parquet")

    """
    merge-parquets.py ${tissue}.cis_qtl.signif_pairs.parquet ${signif_pairs.join(' ')}
    """

}


process postprocess_nominal_get_significant_pairs {

    tag "${tissue} ${phenotype_pcs}"
    publishDir "${params.results}/cis/nominal"
    clusterOptions "${CLUSTER_OPTIONS}"
    container 'library://porchard/default/general:20220107'
    memory '50 GB'
    time '2h'

    input:
    tuple val(tissue), path(nominal_parquets), path(phenotype_groups), path(permutation_results)

    output:
    path("${tissue}.cis_qtl.signif_pairs.parquet")

    """
    post-get-significant-pairs.py --permutations $permutation_results --groups $phenotype_groups --nominal-files ${nominal_parquets.join(' ')} --out ${tissue}.cis_qtl.signif_pairs.parquet
    """

}


process scan_susie {

    tag "${tissue} ${chrom} ${L}"
    publishDir "${params.results}/cis/susie"
    clusterOptions "${CLUSTER_OPTIONS} --gres=gpu:1"
    container "${TENSORQTL_CONTAINER}"
    beforeScript "export PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:100"
    memory '40 GB'
    time '20h'

    input:
    tuple val(tissue), val(chrom), path(phenotypes), path(genotypes), val(phenotype_pcs), path(covariates), path(significant_pairs)
    each L

    output:
    tuple val(tissue), val(L), path("${tissue}.${chrom}.${phenotype_pcs}.${L}.SuSiE.pickle"), path("${tissue}.${chrom}.${phenotype_pcs}.${L}.SuSiE_summary.parquet")

    when:
    (PHENOTYPE_PCS[tissue] == phenotype_pcs) && (L <= SUSIE_MAX_L[tissue])

    script:
    prefix = "${tissue}.${chrom}.${phenotype_pcs}.${L}"

    """
    python3.8 -m tensorqtl --seed 2021 --mode cis_susie --cis_output $significant_pairs --invnorm --maf_threshold $MAF_THRESHOLD --L $L --covariates $covariates ${tissue}.${chrom} $phenotypes $prefix
    """

}

process postprocess_susie_each_L {

    tag "${tissue}"
    publishDir "${params.results}/cis/susie-postprocessed-per-L/${L}"
    clusterOptions "${CLUSTER_OPTIONS}"
    container "${TENSORQTL_CONTAINER}"
    memory { 75.GB * task.attempt }
    time '50h'
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(tissue), val(L), path(susie_pickles), path(susie_parquets)
    each by_cluster_or_gene

    output:
    path("${tissue}*.txt")
    path("${tissue}*.pickle")

    script:
    prefix = by_cluster_or_gene == 'by_cluster' ? (tissue + '.by-cluster.') : (tissue + '.by-gene.')
    by_cluster_flag = by_cluster_or_gene == 'by_cluster' ? '--by-cluster' : ''

    """
    postprocess-cis-sqtl.py --susie-pickle ${susie_pickles.join(' ')} --susie-parquet ${susie_parquets.join(' ')} --prefix $prefix $by_cluster_flag
    """

}


process combine_permutations {

    tag "${tissue}"
    clusterOptions "${CLUSTER_OPTIONS}"
    publishDir "${params.results}/cis/merged-permutations"
    memory '10 GB'
    container 'library://porchard/default/general:20220107'

    input:
    tuple val(tissue), val(pcs), path(permutation_files)

    output:
    tuple val(tissue), val(pcs), path("${tissue}.${pcs}.cis_qtl.txt.gz")

    """
    combine-permutations-and-compute-qvalue-bh-and-nominal-p-thresholds.py ${tissue}.${pcs} ${permutation_files}
    """

}

// if it's a chrom with lots of variants, I think lots of memory is used in outputting table, and this can lead to (1) long run times and (2) OOM errors. Hence chunk_size
process nominal_calculate_p_in_log10_space {

    memory { 140.GB * task.attempt }
    cache 'lenient'
    container 'library://porchard/default/r-arrow:20220117'
    maxForks 30
    tag "${tissue} ${chrom}"
    clusterOptions "${CLUSTER_OPTIONS}"
    time '20h'
    maxRetries 3
    errorStrategy 'retry'

    input:
    tuple val(tissue), val(chrom), path('old.parquet'), path(covariates)

    output:
    tuple val(tissue), val(chrom), path("${tissue}.${chrom}.nominal.parquet")

    """
    #!/usr/bin/env Rscript

    library(arrow)

    message('Reading parquet')
    tmp <- read_parquet('old.parquet')
    message('Reading covariates')
    covariates <- read.table('${covariates}', sep='\\t', head=T, as.is=T, row.names=1)
    dof <- ncol(covariates) - nrow(covariates) - 2
    message('Calculating p-values in log space')
    tmp\$negative_log10_p <- abs((pt(-abs(tmp\$slope / tmp\$slope_se),dof,log.p=TRUE)+log(2))/log(10))
    message('Outputting parquet')
    write_parquet(tmp, '${tissue}.${chrom}.nominal.parquet', chunk_size=1000000)
    """

}

// use characters for start/end to avoid scientific notation
process nominal_calculate_p_in_log10_space_output_text {

    memory { 50.GB * task.attempt }
    cache 'lenient'
    container 'library://porchard/default/r-arrow:20220117'
    maxForks 30
    tag "${tissue} ${chrom}"
    clusterOptions "${CLUSTER_OPTIONS}"
    time '20h'
    maxRetries 3
    errorStrategy 'retry'

    input:
    tuple val(tissue), val(chrom), path('old.parquet'), path(covariates)

    output:
    tuple val(tissue), val(chrom), path("${tissue}.${chrom}.nominal.txt")

    """
    #!/usr/bin/env Rscript

    options(scipen=999)
    library(arrow)

    message('Reading parquet')
    tmp <- read_parquet('old.parquet')
    message('Reading covariates')
    covariates <- read.table('${covariates}', sep='\\t', head=T, as.is=T, row.names=1)
    dof <- ncol(covariates) - nrow(covariates) - 2
    message('Calculating p-values in log space')
    tmp\$negative_log10_p <- abs((pt(-abs(tmp\$slope / tmp\$slope_se),dof,log.p=TRUE)+log(2))/log(10))
    message('Adding new columns')
    COLS <- colnames(tmp)
    variant_id_split <- strsplit(tmp\$variant_id, '_')
    tmp\$chrom <- sapply(variant_id_split, function(x){x[1]})
    tmp\$end <- as.numeric(sapply(variant_id_split, function(x){x[2]}))
    tmp\$start <- tmp\$end - 1
    tmp\$start <- as.character(tmp\$start)
    tmp\$end <- as.character(tmp\$end)
    tmp <- tmp[,c('chrom', 'start', 'end', COLS)]
    colnames(tmp) <- c('#chrom', colnames(tmp)[2:ncol(tmp)])
    message('Outputting txt')
    write.table(tmp, '${tissue}.${chrom}.nominal.txt', sep='\t', quote=F, row.names=F)
    """

}

/*

library(arrow)

message('Reading parquet')
tmp <- read_parquet('old.parquet')
message('Reading covariates')
covariates <- read.table('${covariates}', sep='\\t', head=T, as.is=T, row.names=1)
dof <- ncol(covariates) - nrow(covariates) - 2
message('Calculating p-values in log space')
tmp\$negative_log10_p <- abs((pt(-abs(tmp\$slope / tmp\$slope_se),dof,log.p=TRUE)+log(2))/log(10))
message('Adding new columns')
COLS <- colnames(tmp)
variant_id_split <- strsplit(tmp$variant_id, '_')
tmp$chrom <- sapply(variant_id_split, function(x){x[1]})
tmp$end <- as.numeric(sapply(variant_id_split, function(x){x[2]}))
tmp$start <- tmp$end - 1
tmp <- tmp[,c('chrom', 'start', 'end', COLS)]
colnames(tmp) <- c('#chrom', colnames(tmp)[2:ncol(tmp)])
message('Outputting parquet')
#write.table(tmp, '${tissue}.${chrom}.nominal.txt')


*/




process make_nominal_sorted_bed {

    memory { 50.GB * task.attempt }
    maxForks 20
    container 'docker.io/porchard/general:20220406125608'
    clusterOptions "${CLUSTER_OPTIONS}"
    tag "${tissue} ${chrom}"
    time '10h'
    maxRetries 2
    errorStrategy 'retry'

    input:
    tuple val(tissue), val(chrom), path(parquet)

    output:
    tuple val(tissue), path("${tissue}.${chrom}.tsv")

    """
    tensorqtl-nominal-to-bed.py $parquet > tmp.txt
    head -n 1 tmp.txt > ${tissue}.${chrom}.tsv
    cat tmp.txt | awk 'NR>1' | sort -k1,1 -k2n,2 >> ${tissue}.${chrom}.tsv
    rm tmp.txt
    """

}


process make_nominal_sorted_bed_from_text {

    memory { 50.GB * task.attempt }
    maxForks 20
    container 'docker.io/porchard/general:20220406125608'
    clusterOptions "${CLUSTER_OPTIONS}"
    tag "${tissue} ${chrom}"
    time '10h'
    maxRetries 2
    errorStrategy 'retry'

    input:
    tuple val(tissue), val(chrom), path('tmp.txt')

    output:
    tuple val(tissue), path("${tissue}.${chrom}.tsv")

    """
    head -n 1 tmp.txt > ${tissue}.${chrom}.tsv
    cat tmp.txt | awk 'NR>1' | sort -k1,1 -k2n,2 >> ${tissue}.${chrom}.tsv
    rm tmp.txt
    """

}


process merge_and_tabix_nominal {

    publishDir "${params.results}/cis/nominal/tabix"
    container 'docker.io/porchard/general:20220406125608'
    clusterOptions "${CLUSTER_OPTIONS}"
    tag "${tissue}"
    time '10h'

    input:
    tuple val(tissue),  path(x)

    output:
    path("${tissue}.tsv.gz")
    path("${tissue}.tsv.gz.tbi")

    script:
    fs = x.sort({y -> y.getName().tokenize('.')[1]})

    """
    head -n 1 ${x[0]} > ${tissue}.tsv
    cat ${fs.join(' ')} | grep -v "^#chrom" >> ${tissue}.tsv
    bgzip ${tissue}.tsv
    tabix -p bed ${tissue}.tsv.gz
    """

}


process scan_independent {

    tag "${tissue} ${chrom} ${phenotype_pcs}"
    publishDir "${params.results}/independent"
    clusterOptions "${CLUSTER_OPTIONS} --gres=gpu:1"
    container "${TENSORQTL_CONTAINER}"
    beforeScript "export PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:100; python -m torch.utils.collect_env"
    time '168h'
    memory '60 GB'

    input:
    tuple val(tissue), val(chrom), path(phenotypes), path(phenotype_groups), path(genotypes), val(phenotype_pcs), path(covariates), path(permutation_results)

    output:
    tuple val(tissue), path("${tissue}.${chrom}.${phenotype_pcs}.cis_independent_qtl.txt.gz")

    when:
    PHENOTYPE_PCS[tissue] == phenotype_pcs

    script:
    prefix = "${tissue}.${chrom}.${phenotype_pcs}"

    """
    python3.8 -m tensorqtl --seed 2021 --mode cis_independent --invnorm --cis_output $permutation_results --phenotype_groups $phenotype_groups --maf_threshold $MAF_THRESHOLD --covariates $covariates ${tissue}.${chrom} $phenotypes $prefix
    """

}


process merge_conditional {

    publishDir "${params.results}/independent/merged"
    tag "${tissue} ${phenotype_pcs}"
    container 'docker.io/porchard/general:20220406125608'
    clusterOptions "${CLUSTER_OPTIONS}"

    input:
    tuple val(tissue), path(x)

    output:
    path("${tissue}.tsv")

    """
    zcat ${x[0]} | awk 'NR==1' > header.txt
    cp header.txt ${tissue}.tsv
    zcat ${x.join(' ')} | grep -v -f header.txt >> ${tissue}.tsv
    """

}


workflow permutations {
    PHENOTYPE_GROUP_FILES = Channel.fromPath(PHENOTYPE_GROUP_GLOB).map({it -> it.getName().tokenize('.')[0..1] + [it]}) // tissue, chrom, phenotype groups file
    //PHENOTYPE_FILES.combine(PHENOTYPE_GROUP_FILES, by: [0, 1]).combine(GENOTYPE_FILES, by: [0, 1]).combine(COVARIATE_FILES, by: 0).view()
    tissue_chrom_pheno_geno = PHENOTYPE_FILES.combine(PHENOTYPE_GROUP_FILES, by: [0, 1]).combine(GENOTYPE_FILES, by: [0, 1]).combine(COVARIATE_FILES, by: 0) | scan_permutations
    tissue_chrom_pheno_geno.groupTuple(by: [0, 1]) | combine_permutations
}


workflow nominal {
    permutation_results = Channel.fromPath(PERMUTATION_GLOB).map({it -> it.getName().tokenize('.')[0..1] + [it]}).filter({it -> PHENOTYPE_PCS[it[0]] == it[1]}).map({it -> [it[0], it[2]]}) // tissue, file
    PHENOTYPE_GROUP_FILES = Channel.fromPath(PHENOTYPE_GROUP_GLOB).map({it -> it.getName().tokenize('.')[0..1] + [it]}) // tissue, chrom, phenotype groups file
    tissue_chrom_pheno_geno = PHENOTYPE_FILES.combine(PHENOTYPE_GROUP_FILES, by: [0, 1]).combine(GENOTYPE_FILES, by: [0, 1]).combine(COVARIATE_FILES, by: 0).combine(permutation_results, by: 0) | scan_nominal
    tissue_chrom_pheno_geno.signif_pairs.groupTuple() | merge_significant_pairs // replaces the line below. I believe it's valid to do this on a per-chromosome basis; compare to previous results, should be similar
    // tissue_chrom_pheno_geno.transpose().map({it -> [it[0], it[2]]}).groupTuple().combine(PHENOTYPE_GROUP_FILES, by: 0).combine(permutation_results, by: 0) | postprocess_nominal_get_significant_pairs
    covariate_files = COVARIATE_FILES.filter({it -> it[1] == PHENOTYPE_PCS[it[0]]}).map({it -> [it[0], it[2]]})
    //(tissue_chrom_pheno_geno.transpose().combine(covariate_files, by: 0) | nominal_calculate_p_in_log10_space | make_nominal_sorted_bed).groupTuple() | merge_and_tabix_nominal
    (tissue_chrom_pheno_geno.all_pairs.transpose().combine(covariate_files, by: 0) | nominal_calculate_p_in_log10_space_output_text | make_nominal_sorted_bed_from_text).groupTuple() | merge_and_tabix_nominal
}


workflow susie {
    signif_pairs = Channel.fromPath(params.signif_pairs_glob).map({it -> [it.getName().tokenize('.')[0], it]}) // tissue, signif pairs
    x = scan_susie(PHENOTYPE_FILES.combine(GENOTYPE_FILES, by: [0, 1]).combine(COVARIATE_FILES, by: 0).combine(signif_pairs, by: 0), LS)
    //x.groupTuple(by: [0, 1]).combine(Channel.from(['by_cluster', 'by_gene'])) | postprocess_susie_each_L
    postprocess_susie_each_L(x.groupTuple(by: [0, 1]), ['by_cluster', 'by_gene'])
}


workflow conditional {
    permutation_results = Channel.fromPath(PERMUTATION_GLOB).map({it -> it.getName().tokenize('.')[0..1] + [it]}).filter({it -> PHENOTYPE_PCS[it[0]] == it[1]}).map({it -> [it[0], it[2]]}) // tissue, file
    PHENOTYPE_GROUP_FILES = Channel.fromPath(PHENOTYPE_GROUP_GLOB).map({it -> it.getName().tokenize('.')[0..1] + [it]}) // tissue, chrom, phenotype groups file
    independent_out = scan_independent(PHENOTYPE_FILES.combine(PHENOTYPE_GROUP_FILES, by: [0, 1]).combine(GENOTYPE_FILES, by: [0, 1]).combine(COVARIATE_FILES, by: 0).combine(permutation_results, by: 0))
    independent_out.groupTuple() | merge_conditional
}
