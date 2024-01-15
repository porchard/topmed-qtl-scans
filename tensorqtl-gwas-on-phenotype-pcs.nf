#!/usr/bin/env nextflow

nextflow.enable.dsl=2
MAC_THRESHOLD = params.mac
TENSORQTL_CONTAINER = 'docker.io/porchard/tensorqtl_dev:774fb5b'
CLUSTER_OPTIONS = "--account=dcmb_dept --partition=dcmb_dept"


def covariate_file_to_tissue (f) {
        f.getName().tokenize('.')[0]
}

def phenotype_file_to_tissue (f) {
        f.getName().tokenize('.')[0]
}

process prep {

    publishDir "${params.results}/trans/prep/${tissue}/${chrom}"
    tag "${tissue} ${chrom}"
    container 'docker.io/porchard/general:20220406125608'
    clusterOptions "${CLUSTER_OPTIONS}"
    memory '20 GB'

    input:
    tuple val(chrom), path(genotypes), val(tissue), path('covariates.txt')

    output:
    tuple val(tissue), val(chrom), path("${tissue}.${chrom}.phenotypes.bed.gz"), path("${tissue}.covariates.tsv"), path("${tissue}.genotypes.*")

    """
    # make new phenotype file, consisting of phenotype PCs
    covariate-file-to-pca-phenotype-file.py covariates.txt | gzip -c > ${tissue}.${chrom}.phenotypes.bed.gz
    # drop phenotype_PCs from covariate matrix
    grep -v "^phenotype_PC" covariates.txt > ${tissue}.covariates.tsv
    # filter genotypes to the individuals of interest
    zcat ${tissue}.${chrom}.phenotypes.bed.gz | awk 'NR==1' | cut -f5- | perl -pe 's/\\t/\\n/g' | perl -pe 's/^/0\\t/' > keep.txt
    plink -bfile $chrom --output-chr M --keep-allele-order --make-bed --mac $MAC_THRESHOLD --out ${tissue}.genotypes --keep keep.txt --indiv-sort file keep.txt
    """

}


process run_gwas {

    tag "${tissue} ${chrom}"
    publishDir "${params.results}/trans/gwas/nomerge"
    clusterOptions "${CLUSTER_OPTIONS} --gres=gpu:1"
    beforeScript "export PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:100; python -m torch.utils.collect_env"
    container "${TENSORQTL_CONTAINER}"
    memory '100 GB'

    input:
    tuple val(tissue), val(chrom), path(phenotypes), path(covariates), path(genotypes)

    output:
    tuple val(tissue), path("${tissue}.${chrom}.trans_qtl_pairs.txt.gz")

    script:
    prefix = "${tissue}.${chrom}"

    """
    python3.8 -m tensorqtl --mode trans --invnorm --output_text --maf_threshold 0.0000001 --batch_size 5000 --return_r2 --covariates $covariates ${tissue}.genotypes $phenotypes $prefix
    """

}

process merge_chroms {

    tag "${tissue}"
    publishDir "${params.results}/trans/gwas/merge"
    clusterOptions "${CLUSTER_OPTIONS}"
    container 'docker.io/porchard/general:20220406125608'
    memory '10 GB'

    input:
    tuple val(tissue), path(x)

    output:
    path("${tissue}.trans_qtl_pairs.txt.gz")

    """
    zcat ${x[0]} | awk 'NR==1' > header.txt
    cat header.txt > ${tissue}.trans_qtl_pairs.txt
    zcat ${x.join(' ')} | grep -v -f header.txt >> ${tissue}.trans_qtl_pairs.txt
    gzip ${tissue}.trans_qtl_pairs.txt
    """

}


workflow gwas {
    COVARIATE_FILES = Channel.fromPath(params.covariate_glob).map({it -> [covariate_file_to_tissue(it), it]})
    GENOTYPE_FILES = Channel.fromPath(params.plink_glob).map({it -> [it.getName().tokenize('.')[0], it]}).groupTuple()

    gwas_out = (GENOTYPE_FILES.combine(COVARIATE_FILES) | prep) | run_gwas
    gwas_out.groupTuple() | merge_chroms
}
