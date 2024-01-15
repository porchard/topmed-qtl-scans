# TOPMed e/sQTL scans

## Setup

### Dependencies

* Singularity (v. 3)
* NextFlow (v. >= 21.04.0)

### Configuration

Set up a nextflow config as appropriate for your compute environment. Some general resource parameters are also set in the processes of the individual *.nf files.

Update the paths at the top of the Makefile (first 6 lines) to point to wherever you'll be running scans (ROOT is the current repo directory; DATA is the top directory where input data files (described below) are stored; WORK is where analyses will be run).

### Data formats

* Per-chromosome plink genotype files for autosomes and chrX should be placed in DATA/genotypes, named according to the scheme `{chrom}.{bim,bed,fam}`.
* Tab-separated phenotype files (see tensorQTL README) are named `{PREFIX}.tensorqtl-in.phenotypes.bed.gz`. 
  * For cross-ancestry scans, phenotype files should be in DATA/cross-ancestry/cis-{e,s}qtl and PREFIX should be the tissue name.
  * For per-ancestry scans, phenotype files should be in DATA/ancestry-specific/{e,s}qtl and PREFIX should be `{TISSUE}___{ANCESTRY}`, e.g. Whole_blood___AFR
  * For saturation analysis, phenotype files shoudl be in DATA/saturation/{e,s}qtl and PREFIX should be `{TISSUE}_{NUMBER_SAMPLES}`, e.g. `Whole_blood_2500`
* Tab-separated covariate files (see tensorQTL README) are named `{PREFIX}.tensorqtl-in.{NUMBER_PHENOTYPE_PCS_INCLUDED}.covariates.tsv`. Prefixes follow the above conventions, and files should be in the same file as the corresponding phenotype file.
* Phenotype group files (for sQTL scans, mapping introns --> genes) should be named `{PREFIX}.leafcutter.phenotype_groups.txt`. Prefixes follow the above conventions, and files should be in the same file as the corresponding phenotype file.

A list of variants to exclude from the cis scans should be at DATA/variants-exclude.txt

### Cross mappability files
The following GENCODE v30 cross mappability files (generated using https://github.com/porchard/crossmap-nextflow) should be placed in DATA/crossmap/gencode-v30:
* crossmap-gencode-v30/results/snp-mappability/snp_mappability_100mer_2mismatch.bed.gz .
* crossmap-gencode-v30/results/snp-mappability/snp_mappability_100mer_2mismatch.bed.gz.tbi .
* crossmap-gencode-v30/results/gene-mappability/gene_mappability/gene_mappability.txt .
* crossmap-gencode-v30/results/cross_mappability/crossmap.txt .

The GTF used for the crossmapping pipeline (named gencode.v30.annotation.gtf) should be in DATA/crossmap/gencode-v30 as well.


## Cross-ancestry scans
### cis-eQTLs

1. Prep genotypes and phenotypes (subset to samples of interest, etc): `make cross-ancestry-ciseqtl-genotypes-and-phenotypes-maf001`
2. Run permutation scan: `make cross-ancestry-ciseqtl-permutations-maf001`
3. Run nominal scan: `make cross-ancestry-ciseqtl-nominal-maf001`
4. Run susie: `make cross-ancestry-ciseqtl-nominal-maf001`

Wait for each pipeline to complete before running the next one. Replacing 'maf001' with 'maf0001' runs the whole blood MAF 0.1% scans instead.

Results will be in $(WORK)/cross-ancestry.

### cis-sQTLs

1. Prep genotypes and phenotypes (subset to samples of interest, etc): `make cross-ancestry-cissqtl-genotypes-and-phenotypes-maf001`
2. Run permutation scan: `make cross-ancestry-cissqtl-permutations-maf001`
3. Run nominal scan: `make cross-ancestry-cissqtl-nominal-maf001`
4. Run susie: `make cross-ancestry-cissqtl-nominal-maf001`

Wait for each pipeline to complete before running the next one. Replacing 'maf001' with 'maf0001' runs the whole blood MAF 0.1% scans instead.

Results will be in $(WORK)/cross-ancestry.

### trans-eQTL

1. Run trans-eQTL permutations and fine-mapping: `make cross-ancestry-trans-eqtl-maf005`

Results will be in $(WORK)/cross-ancestry.

### trans-sQTL
1. Run trans-sQTL permutations and fine-mapping: `make cross-ancestry-trans-sqtl-maf005`

Results will be in $(WORK)/cross-ancestry.

### GWAS on gene expression and splicing phenotype PCs
1. `make cross-ancestry-gwas-on-cis-eqtl-phenotype-pcs``
2. `make cross-ancestry-gwas-on-cis-sqtl-phenotype-pcs``

## Ancestry specific scans

Follow the above instructions for cross-ancestry, replacing `cross-ancestry` with `ancestry-specific`.

## Saturation analysis

### cis-eQTLs
1. Prep genotypes and phenotypes (subset to samples of interest, etc): `make saturation-ciseqtl-genotypes-and-phenotypes-maf001`
2. Run permutation scan: `make saturation-ciseqtl-permutations-maf001`
3. Run susie: `make saturation-ciseqtl-susie-maf001`

Wait for each pipeline to complete before running the next one.

Results will be in $(WORK)/saturation.

### trans-eQTLs

1. Run trans-eQTL permutations and fine-mapping: `make saturation-trans-eqtl-maf005`

Results will be in $(WORK)/saturation.
