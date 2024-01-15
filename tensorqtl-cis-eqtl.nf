#!/usr/bin/env nextflow

nextflow.enable.dsl=2
MAF_THRESHOLD = params.maf
LS = [10, 20, 30]
TENSORQTL_CONTAINER = 'docker.io/porchard/tensorqtl_dev:5ea048f'
TENSORQTL_SUSIE_CONTAINER = 'docker.io/porchard/tensorqtl_dev:5ea048f_pyarrow_11' // need this for parquet
CLUSTER_OPTIONS = "--account=dcmb_dept --partition=dcmb_dept"
COVARIATE_GLOB = params.covariate_glob // {tissue}.tensorqtl-in.{pcs}.covariates.tsv
PHENOTYPE_GLOB = params.phenotype_glob // {tissue}.{chrom}.phenotypes.bed.gz
PLINK_GLOB = params.plink_glob // {tissue}.{chrom}.{bim,bam,fam}
PERMUTATION_GLOB = params.permutation_glob // ${tissue}..*

COVARIATE_FILES = Channel.fromPath(COVARIATE_GLOB).map({it -> [covariate_file_to_tissue(it), covariate_file_to_pcs(it), it]}) // tissue, phenotype PCs, covariate file
PHENOTYPE_FILES = Channel.fromPath(PHENOTYPE_GLOB).map({it -> it.getName().tokenize('.')[0..1] + [it]}) // tissue, chrom, phenotype file
GENOTYPE_FILES = Channel.fromPath(PLINK_GLOB).map({it -> it.getName().tokenize('.')[0..1] + [it]}).groupTuple(by: [0, 1]) // tissue, chrom, plink files


PHENOTYPE_PCS = [
    'Lung': '75',
    'Whole_blood': '100',
    'T_cell': '30',
    'PBMC': '30',
    'Monocyte': '30',
    'Nasal_epithelial': '30',
    'Lung___EUR': '20',
    'Monocyte___EUR': '15',
    'Nasal_epithelial___EUR': '20',
    'PBMC___AFR': '20',
    'PBMC___EAS': '10',
    'PBMC___EUR': '20',
    'T_cell___EUR': '15',
    'Whole_blood___AFR': '40',
    'Whole_blood___EUR': '40'
]

SUSIE_MAX_L = [
    'Lung': 20,
    'Whole_blood': 30,
    'T_cell': 10,
    'PBMC': 20,
    'Monocyte': 10,
    'Nasal_epithelial': 10,
    'Lung___EUR': 20,
    'Monocyte___EUR': 10,
    'Nasal_epithelial___EUR': 10,
    'PBMC___AFR': 10,
    'PBMC___EAS': 10,
    'PBMC___EUR': 10,
    'T_cell___EUR': 10,
    'Whole_blood___AFR': 20,
    'Whole_blood___EUR': 20
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
    publishDir "${params.results}/permutations"
    clusterOptions "${CLUSTER_OPTIONS} --gres=gpu:1"
    container "${TENSORQTL_CONTAINER}"
    beforeScript "export PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:100; python -m torch.utils.collect_env"
    memory '40 GB'
    time '10h'

    input:
    tuple val(tissue), val(chrom), path(phenotypes), path(genotypes), val(phenotype_pcs), path(covariates)

    output:
    tuple val(tissue), val(phenotype_pcs), path("${tissue}.${chrom}.${phenotype_pcs}.cis_qtl.txt.gz")

    script:
    prefix = "${tissue}.${chrom}.${phenotype_pcs}"

    """
    python3.8 -m tensorqtl --qvalue_lambda 0 --seed 2021 --invnorm --mode cis --maf_threshold $MAF_THRESHOLD --covariates $covariates ${tissue}.${chrom} $phenotypes $prefix
    """

}


process combine_permutations {

    tag "${tissue} ${phenotype_pcs}"
    memory '10 GB'
    container 'docker.io/porchard/general:20220406125608'
    clusterOptions "${CLUSTER_OPTIONS}"
    publishDir "${params.results}/permutations/merge-chroms"


    input:
    tuple val(tissue), val(phenotype_pcs), path(permutation_files)

    output:
    tuple val(tissue), val(phenotype_pcs), path("${tissue}.${phenotype_pcs}.cis_qtl.txt.gz")

    """
    combine-permutations-and-compute-qvalue-bh.py ${tissue}.${phenotype_pcs} ${permutation_files}
    """

}


process plot_permutations {

    tag "${tissue}"
    memory '10 GB'
    container 'docker.io/porchard/general:20220406125608'
    clusterOptions "${CLUSTER_OPTIONS}"
    publishDir "${params.results}/permutations/plot"

    input:
    tuple val(tissue), path(permutations)

    output:
    path("*.png")

    """
    plot-pcs-vs-egenes.py --prefix ${tissue}. --permutation-files ${permutations.join(' ')} --fdr 0.05
    """

}


process scan_nominal {

    tag "${tissue} ${chrom} ${phenotype_pcs}"
    publishDir "${params.results}/nominal"
    clusterOptions "${CLUSTER_OPTIONS} --gres=gpu:1"
    container "${TENSORQTL_SUSIE_CONTAINER}"
    beforeScript "export PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:100; python -m torch.utils.collect_env"
    memory '70 GB'
    time '4h'

    input:
    tuple val(tissue), val(chrom), path(phenotypes), path(genotypes), val(phenotype_pcs), path(covariates)

    output:
    tuple val(tissue), val(chrom), path("${tissue}.${chrom}.${phenotype_pcs}.cis_qtl_pairs.*.parquet")

    when:
    PHENOTYPE_PCS[tissue] == phenotype_pcs

    script:
    prefix = "${tissue}.${chrom}.${phenotype_pcs}"

    """
    python3.8 -m tensorqtl --seed 2021 --mode cis_nominal --invnorm --maf_threshold $MAF_THRESHOLD --covariates $covariates ${tissue}.${chrom} $phenotypes $prefix
    """

}


process nominal_calculate_p_in_log10_space {

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
    tuple val(tissue), val(chrom), path("${tissue}.${chrom}.nominal.parquet")

    """
    #!/usr/bin/env Rscript

    library(arrow)

    tmp <- read_parquet('old.parquet')
    covariates <- read.table('${covariates}', sep='\\t', head=T, as.is=T, row.names=1)
    dof <- ncol(covariates) - nrow(covariates) - 2
    tmp\$negative_log10_p <- abs((pt(-abs(tmp\$slope / tmp\$slope_se),dof,log.p=TRUE)+log(2))/log(10))
    write_parquet(tmp, '${tissue}.${chrom}.nominal.parquet')
    """

}


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


process merge_and_tabix_nominal {

    publishDir "${params.results}/nominal/tabix"
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



process scan_susie {

    tag "${tissue} ${chrom} ${L}"
    publishDir "${params.results}/susie/before-merge"
    clusterOptions "${CLUSTER_OPTIONS} --gres=gpu:1"
    container "${TENSORQTL_SUSIE_CONTAINER}"
    beforeScript "export PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:100; python -m torch.utils.collect_env"
    memory '60 GB'
    time '24h'

    input:
    tuple val(tissue), val(chrom), path(phenotypes), path(genotypes), val(phenotype_pcs), path(covariates), path(permutations)
    each L

    output:
    tuple val(tissue), val(phenotype_pcs), path("${tissue}.${chrom}.${phenotype_pcs}.${L}.susie.pickle"), emit: pickle
    path("${tissue}.${chrom}.${phenotype_pcs}.${L}.SuSiE_summary.parquet")

    when:
    (PHENOTYPE_PCS[tissue] == phenotype_pcs) && (L <= SUSIE_MAX_L[tissue])

    script:
    prefix = "${tissue}.${chrom}.${phenotype_pcs}.${L}"

    """
    python3.8 -m tensorqtl --seed 2021 --mode cis_susie --cis_output $permutations --invnorm --maf_threshold $MAF_THRESHOLD --L $L --covariates $covariates ${tissue}.${chrom} $phenotypes $prefix
    mv ${tissue}.${chrom}.${phenotype_pcs}.${L}.SuSiE.pickle ${tissue}.${chrom}.${phenotype_pcs}.${L}.susie.pickle
    """

}


process select_susie_L {

    tag "${tissue}"
    publishDir "${params.results}/susie/merge-and-select-L"
    container 'docker.io/porchard/general:20220406125608'
    clusterOptions "${CLUSTER_OPTIONS}"
    memory '100 GB'
    time '24h'

    input:
    tuple val(tissue), val(phenotype_pcs), path(pickles)

    output:
    tuple val(tissue), val(phenotype_pcs), path("${tissue}.${phenotype_pcs}.susie.pickle")

    """
    susie-pickle select-l --out ${tissue}.${phenotype_pcs}.susie.pickle --pickles ${pickles.join(' ')}
    """

}


process susie_pickle_to_summary {

    tag "${tissue}"
    publishDir "${params.results}/susie/merge-and-select-L"
    container 'docker.io/porchard/general:20220406125608'
    clusterOptions "${CLUSTER_OPTIONS}"
    memory '30 GB'

    input:
    tuple val(tissue), val(phenotype_pcs), path(pickle)

    output:
    tuple val(tissue), val(phenotype_pcs), path("${tissue}.${phenotype_pcs}.cs.txt"), path("${tissue}.${phenotype_pcs}.converged.txt")

    """
    susie-pickle-2 summarize --pickle $pickle --single-effect-pips --deduplicate --prefix ${tissue}.${phenotype_pcs}.
    """

}


process scan_independent {

    tag "${tissue} ${chrom} ${phenotype_pcs}"
    publishDir "${params.results}/independent"
    clusterOptions "${CLUSTER_OPTIONS} --gres=gpu:1"
    container "${TENSORQTL_SUSIE_CONTAINER}"
    beforeScript "export PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:100; python -m torch.utils.collect_env"
    time '20h'
    memory '60 GB'

    input:
    tuple val(tissue), val(chrom), path(phenotypes), path(genotypes), val(phenotype_pcs), path(covariates), path(permutation_results)

    output:
    tuple val(tissue), path("${tissue}.${chrom}.${phenotype_pcs}.cis_independent_qtl.txt.gz")

    when:
    PHENOTYPE_PCS[tissue] == phenotype_pcs

    script:
    prefix = "${tissue}.${chrom}.${phenotype_pcs}"

    """
    python3.8 -m tensorqtl --seed 2021 --mode cis_independent --invnorm --cis_output $permutation_results --maf_threshold $MAF_THRESHOLD --covariates $covariates ${tissue}.${chrom} $phenotypes $prefix
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


process residualize_phenotypes_and_genotypes {
    
    tag "${tissue} ${chrom} ${phenotype_pcs}"
    publishDir "${params.results}/cis/residualize"
    clusterOptions "${CLUSTER_OPTIONS} --gres=gpu:1"
    container "${TENSORQTL_CONTAINER}"
    beforeScript "export PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:100; python -m torch.utils.collect_env"
    memory '150 GB'
    time '6h'

    input:
    tuple val(tissue), val(chrom), path(phenotypes), path(genotypes), val(phenotype_pcs), path(covariates)

    output:
    path("${tissue}.${chrom}.*.parquet")

    when:
    PHENOTYPE_PCS[tissue] == phenotype_pcs

    script:
    prefix = "${tissue}.${chrom}.${phenotype_pcs}"

    """
    residualize-invnorm.py --covariates $covariates --phenotypes $phenotypes --plink-prefix ${tissue}.${chrom} --prefix ${tissue}.${chrom}.
    """

}


workflow permutations {
    permutations_out = PHENOTYPE_FILES.combine(GENOTYPE_FILES, by: [0, 1]).combine(COVARIATE_FILES, by: 0) | scan_permutations
    combined_permutations = permutations_out.groupTuple(by: [0, 1]) | combine_permutations
    combined_permutations.map({it -> [it[0], it[2]]}).groupTuple() | plot_permutations
}

workflow nominal {
    COVARIATE_FILES = COVARIATE_FILES.filter({it -> it[1] == PHENOTYPE_PCS[it[0]]})
    tissue_chrom_pheno_geno = PHENOTYPE_FILES.combine(GENOTYPE_FILES, by: [0, 1]).combine(COVARIATE_FILES, by: 0) | scan_nominal
    covariate_files = COVARIATE_FILES.map({it -> [it[0], it[2]]})
    (tissue_chrom_pheno_geno.transpose().combine(covariate_files, by: 0) | nominal_calculate_p_in_log10_space | make_nominal_sorted_bed).groupTuple() | merge_and_tabix_nominal
}

workflow susie {
    COVARIATE_FILES = COVARIATE_FILES.filter({it -> it[1] == PHENOTYPE_PCS[it[0]]})
    PERMUTATIONS = Channel.fromPath(PERMUTATION_GLOB).map({it -> [it.getName().tokenize('.')[0], it]})
    tissue_chrom_pheno_geno = PHENOTYPE_FILES.combine(GENOTYPE_FILES, by: [0, 1]).combine(COVARIATE_FILES, by: 0).combine(PERMUTATIONS, by: 0)
    susie_out = scan_susie(tissue_chrom_pheno_geno, LS)
    susie_out.pickle.groupTuple(by: [0, 1]) | select_susie_L | susie_pickle_to_summary
}

workflow residualize {
    COVARIATE_FILES = COVARIATE_FILES.filter({it -> it[1] == PHENOTYPE_PCS[it[0]]})
    tissue_chrom_pheno_geno = PHENOTYPE_FILES.combine(GENOTYPE_FILES, by: [0, 1]).combine(COVARIATE_FILES, by: 0) | residualize_phenotypes_and_genotypes
}

workflow conditional {
    COVARIATE_FILES = COVARIATE_FILES.filter({it -> it[1] == PHENOTYPE_PCS[it[0]]})
    PERMUTATIONS = Channel.fromPath(PERMUTATION_GLOB).map({it -> [it.getName().tokenize('.')[0], it]})
    independent_out = PHENOTYPE_FILES.combine(GENOTYPE_FILES, by: [0, 1]).combine(COVARIATE_FILES, by: 0).combine(PERMUTATIONS, by: 0) | scan_independent
    independent_out.groupTuple() | merge_conditional
}
