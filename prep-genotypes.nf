#!/usr/bin/env nextflow

nextflow.enable.dsl=2
MAF_THRESHOLD = params.maf
PHENOTYPE_GLOB = params.phenotype_glob // {tissue}.tensorqtl-in.phenotypes.bed.gz
PLINK_GLOB = params.plink_glob // {chrom}.{bim/bam/fam} plink files 
EXCLUDE_VARIANTS = params.exclude_variants // list of variant IDs to exclude
CLUSTER_OPTIONS = "--account=dcmb_dept --partition=dcmb_dept"


process cis_eqtl_get_genotypes_and_phenotypes {

    tag "${tissue} ${chrom}"
    container 'docker.io/porchard/general:20220406125608'
    clusterOptions "${CLUSTER_OPTIONS}"
    memory '20 GB'
    publishDir "${params.results}/genotypes", pattern: "*.{bim,bed,fam}"
    publishDir "${params.results}/phenotypes", pattern: "*.phenotypes.bed.gz"

    input:
    tuple val(chrom), path(genotypes), val(tissue), path('phenotypes.tmp.bed.gz'), path(exclude_variants)

    output:
    path("${tissue}.${chrom}.phenotypes.bed.gz")
    path("${tissue}.${chrom}.bim")
    path("${tissue}.${chrom}.bed")
    path("${tissue}.${chrom}.fam")

    // when:
    // ['chr10', 'chr13', 'chr5', 'chr7', 'chrX'].contains(chrom) // these ones have long indels that had to be removed

    """
    # Filter phenotypes to the chromosome of interest
    zcat phenotypes.tmp.bed.gz | awk 'NR==1' > ${tissue}.${chrom}.phenotypes.bed
    zcat phenotypes.tmp.bed.gz | grep -w ${chrom} | perl -pe 's/^chr//' >> ${tissue}.${chrom}.phenotypes.bed
    gzip ${tissue}.${chrom}.phenotypes.bed
    # filter genotypes to the individuals of interest
    zcat ${tissue}.${chrom}.phenotypes.bed.gz | awk 'NR==1' | cut -f5- | perl -pe 's/\\t/\\n/g' | perl -pe 's/^/0\\t/' > keep.txt
    plink -bfile $chrom --output-chr M --keep-allele-order --make-bed --maf $MAF_THRESHOLD --out ${tissue}.${chrom} --keep keep.txt --exclude $exclude_variants --indiv-sort file keep.txt
    """

}


process cis_sqtl_get_genotypes_and_phenotypes_and_groups {

    tag "${tissue} ${chrom}"
    container 'library://porchard/default/general:20220107'
    clusterOptions "${CLUSTER_OPTIONS}"
    memory '10 GB'
    publishDir "${params.results}/genotypes", pattern: "*.{bim,bed,fam}"
    publishDir "${params.results}/phenotypes", pattern: "*.phenotypes.bed.gz"
    publishDir "${params.results}/phenotype-groups", pattern: "*.phenotype-groups.txt"

    input:
    tuple val(chrom), path(genotypes), val(tissue), path('phenotypes.tmp.bed.gz'), path('phenotype-groups.txt'), path(exclude_variants)

    output:
    path("${tissue}.${chrom}.phenotypes.bed.gz")
    path("${tissue}.${chrom}.phenotype-groups.txt")
    path("${tissue}.${chrom}.bim")
    path("${tissue}.${chrom}.bed")
    path("${tissue}.${chrom}.fam")

    """
    # Filter phenotypes to the chromosome of interest
    zcat phenotypes.tmp.bed.gz | awk 'NR==1' > ${tissue}.${chrom}.phenotypes.tmp
    zcat phenotypes.tmp.bed.gz | grep -w ${chrom} | perl -pe 's/^chr//' >> ${tissue}.${chrom}.phenotypes.tmp
    resort-tensorqtl-phenotypes.py ${tissue}.${chrom}.phenotypes.tmp > ${tissue}.${chrom}.phenotypes.bed
    gzip ${tissue}.${chrom}.phenotypes.bed
    # filter phenotype groups to the chromosome of interest
    cat phenotype-groups.txt | grep "^${chrom}:" > ${tissue}.${chrom}.phenotype-groups.txt
    #zcat ${tissue}.${chrom}.phenotypes.bed.gz | awk 'NR>1' | cut -f4 > all-phenotypes.txt
    #cat all-phenotypes.txt | perl -pe 's/:/\\t/g' | cut -f5 > phenotype-groups.tmp
    #paste all-phenotypes.txt phenotype-groups.tmp > ${tissue}.${chrom}.phenotype-groups.txt
    # filter genotypes to the individuals of interest
    zcat ${tissue}.${chrom}.phenotypes.bed.gz | awk 'NR==1' | cut -f5- | perl -pe 's/\\t/\\n/g' | perl -pe 's/^/0\\t/' > keep.txt
    plink -bfile $chrom --output-chr M --keep-allele-order --make-bed --maf $MAF_THRESHOLD --out ${tissue}.${chrom} --keep keep.txt --exclude $exclude_variants --indiv-sort file keep.txt
    """

}


workflow ciseqtl {
    exclude_variants = Channel.fromPath(EXCLUDE_VARIANTS)
    PHENOTYPE_FILES = Channel.fromPath(PHENOTYPE_GLOB).map({it -> [it.getName().tokenize('.')[0], it]}) // tissue, phenotypes
    GENOTYPE_FILES = Channel.fromPath(PLINK_GLOB).map({it -> [it.getName().tokenize('.')[0], it]}).groupTuple() // chrom, plink files
    GENOTYPE_FILES.combine(PHENOTYPE_FILES).combine(exclude_variants) | cis_eqtl_get_genotypes_and_phenotypes
}

workflow cissqtl {
    exclude_variants = Channel.fromPath(EXCLUDE_VARIANTS)
    PHENOTYPE_GROUP_GLOB = params.phenotype_group_glob
    PHENOTYPE_GROUP_FILES = Channel.fromPath(PHENOTYPE_GROUP_GLOB).map({it -> [it.getName().tokenize('.')[0], it]}) // tissue, phenotype groups
    PHENOTYPE_FILES = Channel.fromPath(PHENOTYPE_GLOB).map({it -> [it.getName().tokenize('.')[0], it]}) // tissue, phenotypes
    GENOTYPE_FILES = Channel.fromPath(PLINK_GLOB).map({it -> [it.getName().tokenize('.')[0], it]}).groupTuple() // chrom, plink files
    GENOTYPE_FILES.combine(PHENOTYPE_FILES.combine(PHENOTYPE_GROUP_FILES, by: 0)).combine(exclude_variants) | cis_sqtl_get_genotypes_and_phenotypes_and_groups
}