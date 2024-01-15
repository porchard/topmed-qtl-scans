#!/usr/bin/env nextflow

nextflow.enable.dsl=2
MAF_THRESHOLD = params.maf
LS = [10, 20]
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
    'Whole_blood_500': (10..100).step(10).collect({x -> x.toString()}),
    'Whole_blood_1000': (10..100).step(10).collect({x -> x.toString()}),
    'Whole_blood_1500': (10..100).step(10).collect({x -> x.toString()}),
    'Whole_blood_2000': (10..100).step(10).collect({x -> x.toString()}),
    'Whole_blood_2500': (10..100).step(10).collect({x -> x.toString()}),
    'Whole_blood_3000': (10..100).step(10).collect({x -> x.toString()}),
    'Whole_blood_3500': (10..100).step(10).collect({x -> x.toString()}),
    'Whole_blood_4000': (10..100).step(10).collect({x -> x.toString()}),
    'Whole_blood_4500': (10..100).step(10).collect({x -> x.toString()}),
    'Whole_blood_5000': (10..100).step(10).collect({x -> x.toString()}),
    'Whole_blood_5500': (10..100).step(10).collect({x -> x.toString()}),
    'Whole_blood_6000': (10..100).step(10).collect({x -> x.toString()}),
    'Whole_blood_6500': (10..100).step(10).collect({x -> x.toString()})
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
    PHENOTYPE_PCS[tissue].contains(phenotype_pcs)

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


workflow permutations {
    permutations_out = PHENOTYPE_FILES.combine(GENOTYPE_FILES, by: [0, 1]).combine(COVARIATE_FILES, by: 0) | scan_permutations
    combined_permutations = permutations_out.groupTuple(by: [0, 1]) | combine_permutations
    combined_permutations.map({it -> [it[0], it[2]]}).groupTuple() | plot_permutations
}

workflow susie {
    PERMUTATIONS = Channel.fromPath(PERMUTATION_GLOB).map({it -> it.getName().tokenize('.')[0..1] + [it]})
    COVARIATES_AND_PERMUTATIONS = COVARIATE_FILES.combine(PERMUTATIONS, by: [0, 1]) // tissue, pcs, covariates, permutations
    tissue_chrom_pheno_geno = PHENOTYPE_FILES.combine(GENOTYPE_FILES, by: [0, 1]).combine(COVARIATES_AND_PERMUTATIONS, by: 0)
    susie_out = scan_susie(tissue_chrom_pheno_geno, LS)
    susie_out.pickle.groupTuple(by: [0, 1]) | select_susie_L | susie_pickle_to_summary
}