ROOT=/home/porchard
BIN=$(ROOT)/bin
SCRATCH=/scratch/dcmb_dept_root/dcmb_dept/porchard
DATA=$(SCRATCH)/data
WORK=$(SCRATCH)/work
ANALYSIS=$(WORK)/$@


define NL


endef




# cis-eQTL
cross-ancestry-ciseqtl-genotypes-and-phenotypes-maf0001: ANALYSIS=$(WORK)/cross-ancestry/cis-eqtl/genotypes-and-phenotypes/maf0001
cross-ancestry-ciseqtl-genotypes-and-phenotypes-maf0001:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --phenotype_glob '$(DATA)/cross-ancestry/cis-eqtl/Whole*.phenotypes.bed.gz' --maf 0.001 --plink_glob '$(DATA)/genotypes/chr*' --exclude_variants $(DATA)/variants-exclude.txt -entry ciseqtl $(ROOT)/prep-genotypes.nf &

cross-ancestry-ciseqtl-genotypes-and-phenotypes-maf001: ANALYSIS=$(WORK)/cross-ancestry/cis-eqtl/genotypes-and-phenotypes/maf001
cross-ancestry-ciseqtl-genotypes-and-phenotypes-maf001:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --phenotype_glob '$(DATA)/cross-ancestry/cis-eqtl/*.phenotypes.bed.gz' --maf 0.01 --plink_glob '$(DATA)/genotypes/chr*' --exclude_variants $(DATA)/variants-exclude.txt -entry ciseqtl $(ROOT)/prep-genotypes.nf &

cross-ancestry-ciseqtl-permutations-maf001: GP=$(WORK)/cross-ancestry/cis-eqtl/genotypes-and-phenotypes/maf001/results
cross-ancestry-ciseqtl-permutations-maf001: ANALYSIS=$(WORK)/cross-ancestry/cis-eqtl/permutations/maf001
cross-ancestry-ciseqtl-permutations-maf001:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --covariate_glob '$(DATA)/cross-ancestry/cis-eqtl/*.covariates.tsv' --phenotype_glob '$(GP)/phenotypes/*.phenotypes.bed.gz' --maf 0.01 --plink_glob '$(GP)/genotypes/*' -entry permutations $(ROOT)/tensorqtl-cis-eqtl.nf &

cross-ancestry-ciseqtl-permutations-maf0001: GP=$(WORK)/cross-ancestry/cis-eqtl/genotypes-and-phenotypes/maf0001/results
cross-ancestry-ciseqtl-permutations-maf0001: ANALYSIS=$(WORK)/cross-ancestry/cis-eqtl/permutations/maf0001
cross-ancestry-ciseqtl-permutations-maf0001:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --covariate_glob '$(DATA)/cross-ancestry/cis-eqtl/*100.covariates.tsv' --phenotype_glob '$(GP)/phenotypes/*.phenotypes.bed.gz' --maf 0.001 --plink_glob '$(GP)/genotypes/*' -entry permutations $(ROOT)/tensorqtl-cis-eqtl.nf &

cross-ancestry-ciseqtl-nominal-maf0001: GP=$(WORK)/cross-ancestry/cis-eqtl/genotypes-and-phenotypes/maf0001/results
cross-ancestry-ciseqtl-nominal-maf0001: ANALYSIS=$(WORK)/cross-ancestry/cis-eqtl/nominal/maf0001
cross-ancestry-ciseqtl-nominal-maf0001:
	mkdir -p $(ANALYSIS)/covariates
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --covariate_glob '$(DATA)/cross-ancestry/cis-eqtl/*.covariates.tsv' --phenotype_glob '$(GP)/phenotypes/*.phenotypes.bed.gz' --maf 0.001 --plink_glob '$(GP)/genotypes/*' -entry nominal $(ROOT)/tensorqtl-cis-eqtl.nf &

cross-ancestry-ciseqtl-nominal-maf001: GP=$(WORK)/cross-ancestry/cis-eqtl/genotypes-and-phenotypes/maf001/results
cross-ancestry-ciseqtl-nominal-maf001: ANALYSIS=$(WORK)/cross-ancestry/cis-eqtl/nominal/maf001
cross-ancestry-ciseqtl-nominal-maf001:
	mkdir -p $(ANALYSIS)/covariates
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --covariate_glob '$(DATA)/cross-ancestry/cis-eqtl/*.covariates.tsv' --phenotype_glob '$(GP)/phenotypes/*.phenotypes.bed.gz' --maf 0.01 --plink_glob '$(GP)/genotypes/*' -entry nominal $(ROOT)/tensorqtl-cis-eqtl.nf &

cross-ancestry-ciseqtl-susie-maf001: GP=$(WORK)/cross-ancestry/cis-eqtl/genotypes-and-phenotypes/maf001/results
cross-ancestry-ciseqtl-susie-maf001: ANALYSIS=$(WORK)/cross-ancestry/cis-eqtl/susie/maf001
cross-ancestry-ciseqtl-susie-maf001:
	mkdir -p $(ANALYSIS)/permutations
	cp $(WORK)/cross-ancestry/cis-eqtl/permutations/maf001/results/permutations/merge-chroms/Whole_blood.100.cis_qtl.txt.gz $(ANALYSIS)/permutations/
	cp $(WORK)/cross-ancestry/cis-eqtl/permutations/maf001/results/permutations/merge-chroms/Lung.75.cis_qtl.txt.gz $(ANALYSIS)/permutations/
	$(foreach t,T_cell Monocyte PBMC Nasal_epithelial,cp $(WORK)/cross-ancestry/cis-eqtl/permutations/maf001/results/permutations/merge-chroms/$(t).30.cis_qtl.txt.gz $(ANALYSIS)/permutations/$(NL))
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --permutation_glob '$(ANALYSIS)/permutations/*' --covariate_glob '$(DATA)/cross-ancestry/cis-eqtl/*.covariates.tsv' --phenotype_glob '$(GP)/phenotypes/*.phenotypes.bed.gz' --maf 0.01 --plink_glob '$(GP)/genotypes/*' -entry susie $(ROOT)/tensorqtl-cis-eqtl.nf &

cross-ancestry-ciseqtl-susie-maf0001: GP=$(WORK)/cross-ancestry/cis-eqtl/genotypes-and-phenotypes/maf0001/results
cross-ancestry-ciseqtl-susie-maf0001: ANALYSIS=$(WORK)/cross-ancestry/cis-eqtl/susie/maf0001
cross-ancestry-ciseqtl-susie-maf0001:
	mkdir -p $(ANALYSIS)/permutations
	cp $(WORK)/cross-ancestry/cis-eqtl/permutations/maf0001/results/permutations/merge-chroms/Whole_blood.100.cis_qtl.txt.gz $(ANALYSIS)/permutations/
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --permutation_glob '$(ANALYSIS)/permutations/*' --covariate_glob '$(DATA)/cross-ancestry/cis-eqtl/*.covariates.tsv' --phenotype_glob '$(GP)/phenotypes/*.phenotypes.bed.gz' --maf 0.001 --plink_glob '$(GP)/genotypes/*' -entry susie $(ROOT)/tensorqtl-cis-eqtl.nf &

cross-ancestry-ciseqtl-conditional-maf001: GP=$(WORK)/cross-ancestry/cis-eqtl/genotypes-and-phenotypes/maf001/results
cross-ancestry-ciseqtl-conditional-maf001: ANALYSIS=$(WORK)/cross-ancestry/cis-eqtl/conditional/maf001
cross-ancestry-ciseqtl-conditional-maf001:
	mkdir -p $(ANALYSIS)/permutations
	cp $(WORK)/cross-ancestry/cis-eqtl/permutations/maf001/results/permutations/merge-chroms/Whole_blood.100.cis_qtl.txt.gz $(ANALYSIS)/permutations/
	cp $(WORK)/cross-ancestry/cis-eqtl/permutations/maf001/results/permutations/merge-chroms/Lung.75.cis_qtl.txt.gz $(ANALYSIS)/permutations/
	$(foreach t,T_cell Monocyte PBMC Nasal_epithelial,cp $(WORK)/cross-ancestry/cis-eqtl/permutations/maf001/results/permutations/merge-chroms/$(t).30.cis_qtl.txt.gz $(ANALYSIS)/permutations/$(NL))
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --permutation_glob '$(ANALYSIS)/permutations/*' --covariate_glob '$(DATA)/cross-ancestry/cis-eqtl/*.covariates.tsv' --phenotype_glob '$(GP)/phenotypes/*.phenotypes.bed.gz' --maf 0.01 --plink_glob '$(GP)/genotypes/*' -entry conditional $(ROOT)/tensorqtl-cis-eqtl.nf &

cross-ancestry-ciseqtl-conditional-maf0001: GP=$(WORK)/cross-ancestry/cis-eqtl/genotypes-and-phenotypes/maf0001/results
cross-ancestry-ciseqtl-conditional-maf0001: ANALYSIS=$(WORK)/cross-ancestry/cis-eqtl/conditional/maf0001
cross-ancestry-ciseqtl-conditional-maf0001:
	mkdir -p $(ANALYSIS)/permutations
	cp $(WORK)/cross-ancestry/cis-eqtl/permutations/maf0001/results/permutations/merge-chroms/Whole_blood.100.cis_qtl.txt.gz $(ANALYSIS)/permutations/
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --permutation_glob '$(ANALYSIS)/permutations/*' --covariate_glob '$(DATA)/cross-ancestry/cis-eqtl/*.covariates.tsv' --phenotype_glob '$(GP)/phenotypes/*.phenotypes.bed.gz' --maf 0.001 --plink_glob '$(GP)/genotypes/*' -entry conditional $(ROOT)/tensorqtl-cis-eqtl.nf &


# cis-sQTL
cross-ancestry-cissqtl-genotypes-and-phenotypes-maf001: ANALYSIS=$(WORK)/cross-ancestry/cis-sqtl/genotypes-and-phenotypes/maf001
cross-ancestry-cissqtl-genotypes-and-phenotypes-maf001:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --phenotype_glob '$(DATA)/cross-ancestry/cis-sqtl/*.phenotypes.bed.gz' --phenotype_group_glob '$(DATA)/cross-ancestry/cis-sqtl/*.phenotype_groups.txt' --maf 0.01 --plink_glob '$(DATA)/genotypes/chr*' --exclude_variants $(DATA)/variants-exclude.txt -entry cissqtl $(ROOT)/prep-genotypes.nf &

cross-ancestry-cissqtl-genotypes-and-phenotypes-maf0001: ANALYSIS=$(WORK)/cross-ancestry/cis-sqtl/genotypes-and-phenotypes/maf0001
cross-ancestry-cissqtl-genotypes-and-phenotypes-maf0001:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --phenotype_glob '$(DATA)/cross-ancestry/cis-sqtl/Whole_blood*.phenotypes.bed.gz' --phenotype_group_glob '$(DATA)/cross-ancestry/cis-sqtl/Whole_blood*.phenotype_groups.txt' --maf 0.001 --plink_glob '$(DATA)/genotypes/chr*' --exclude_variants $(DATA)/variants-exclude.txt -entry cissqtl $(ROOT)/prep-genotypes.nf &

cross-ancestry-cissqtl-permutations-maf001: GP=$(WORK)/cross-ancestry/cis-sqtl/genotypes-and-phenotypes/maf001/results
cross-ancestry-cissqtl-permutations-maf001: ANALYSIS=$(WORK)/cross-ancestry/cis-sqtl/permutations/maf001
cross-ancestry-cissqtl-permutations-maf001:
	mkdir -p $(ANALYSIS)/data
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --covariate_glob '$(DATA)/cross-ancestry/cis-sqtl/*.covariates.tsv' --phenotype_glob '$(GP)/phenotypes/*.phenotypes.bed.gz' --phenotype_group_glob '$(GP)/phenotype-groups/*' --maf 0.01 --plink_glob '$(GP)/genotypes/*' -entry permutations $(ROOT)/tensorqtl-cis-sqtl.nf &

cross-ancestry-cissqtl-permutations-maf0001: GP=$(WORK)/cross-ancestry/cis-sqtl/genotypes-and-phenotypes/maf0001/results
cross-ancestry-cissqtl-permutations-maf0001: ANALYSIS=$(WORK)/cross-ancestry/cis-sqtl/permutations/maf0001
cross-ancestry-cissqtl-permutations-maf0001:
	mkdir -p $(ANALYSIS)/data
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --covariate_glob '$(DATA)/cross-ancestry/cis-sqtl/Whole*.covariates.tsv' --phenotype_glob '$(GP)/phenotypes/*.phenotypes.bed.gz' --phenotype_group_glob '$(GP)/phenotype-groups/*' --maf 0.001 --plink_glob '$(GP)/genotypes/*' -entry permutations $(ROOT)/tensorqtl-cis-sqtl.nf &

cross-ancestry-cissqtl-nominal-maf001: GP=$(WORK)/cross-ancestry/cis-sqtl/genotypes-and-phenotypes/maf001/results
cross-ancestry-cissqtl-nominal-maf001: ANALYSIS=$(WORK)/cross-ancestry/cis-sqtl/nominal/maf001
cross-ancestry-cissqtl-nominal-maf001:
	mkdir -p $(ANALYSIS)/data
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --covariate_glob '$(DATA)/cross-ancestry/cis-sqtl/*.10.covariates.tsv' --phenotype_glob '$(GP)/phenotypes/*.phenotypes.bed.gz' --phenotype_group_glob '$(GP)/phenotype-groups/*' --permutation_glob '$(WORK)/cross-ancestry/cis-sqtl/permutations/maf001/results/cis/merged-permutations/*' --maf 0.01 --plink_glob '$(GP)/genotypes/*' -entry nominal $(ROOT)/tensorqtl-cis-sqtl.nf &

cross-ancestry-cissqtl-nominal-maf0001: GP=$(WORK)/cross-ancestry/cis-sqtl/genotypes-and-phenotypes/maf0001/results
cross-ancestry-cissqtl-nominal-maf0001:  ANALYSIS=$(WORK)/cross-ancestry/cis-sqtl/nominal/maf0001
cross-ancestry-cissqtl-nominal-maf0001:
	mkdir -p $(ANALYSIS)/data
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --covariate_glob '$(DATA)/cross-ancestry/cis-sqtl/Whole*.10.covariates.tsv' --phenotype_glob '$(GP)/phenotypes/*.phenotypes.bed.gz' --phenotype_group_glob '$(GP)/phenotype-groups/*' --permutation_glob '$(WORK)/cross-ancestry/cis-sqtl/permutations/maf0001/results/cis/merged-permutations/*' --maf 0.001 --plink_glob '$(GP)/genotypes/*' -entry nominal $(ROOT)/tensorqtl-cis-sqtl.nf &

cross-ancestry-cissqtl-susie-maf001: GP=$(WORK)/cross-ancestry/cis-sqtl/genotypes-and-phenotypes/maf001/results
cross-ancestry-cissqtl-susie-maf001: ANALYSIS=$(WORK)/cross-ancestry/cis-sqtl/susie/maf001
cross-ancestry-cissqtl-susie-maf001:
	#mkdir -p $(ANALYSIS)/data
	$(foreach t,T_cell Monocyte PBMC Nasal_epithelial Whole_blood Lung,#ln -s $(WORK)/cross-ancestry/cis-sqtl/nominal/maf001/results/cis/nominal/$(t).cis_qtl.signif_pairs.parquet $(ANALYSIS)/data/$(NL))
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --signif_pairs_glob '$(ANALYSIS)/data/*.parquet' --covariate_glob '$(DATA)/cross-ancestry/cis-sqtl/*.10.covariates.tsv' --phenotype_glob '$(GP)/phenotypes/*.phenotypes.bed.gz' --maf 0.01 --plink_glob '$(GP)/genotypes/*' -entry susie $(ROOT)/tensorqtl-cis-sqtl.nf &

cross-ancestry-cissqtl-susie-maf0001: GP=$(WORK)/cross-ancestry/cis-sqtl/genotypes-and-phenotypes/maf0001/results
cross-ancestry-cissqtl-susie-maf0001: ANALYSIS=$(WORK)/cross-ancestry/cis-sqtl/susie/maf0001
cross-ancestry-cissqtl-susie-maf0001:
	mkdir -p $(ANALYSIS)/data
	$(foreach t,Whole_blood,ln -s $(WORK)/cross-ancestry/cis-sqtl/nominal/maf0001/results/cis/nominal/$(t).cis_qtl.signif_pairs.parquet $(ANALYSIS)/data/$(NL))
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --signif_pairs_glob '$(ANALYSIS)/data/*signif_pairs.parquet' --covariate_glob '$(DATA)/cross-ancestry/cis-sqtl/*.10.covariates.tsv' --phenotype_glob '$(GP)/phenotypes/*.phenotypes.bed.gz' --maf 0.001 --plink_glob '$(GP)/genotypes/*' -entry susie $(ROOT)/tensorqtl-cis-sqtl.nf &

cross-ancestry-cissqtl-conditional-maf001: GP=$(WORK)/cross-ancestry/cis-sqtl/genotypes-and-phenotypes/maf001/results
cross-ancestry-cissqtl-conditional-maf001: ANALYSIS=$(WORK)/cross-ancestry/cis-sqtl/conditional/maf001
cross-ancestry-cissqtl-conditional-maf001:
	mkdir -p $(ANALYSIS)/data
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --covariate_glob '$(DATA)/cross-ancestry/cis-sqtl/*.10.covariates.tsv' --phenotype_glob '$(GP)/phenotypes/*.phenotypes.bed.gz' --phenotype_group_glob '$(GP)/phenotype-groups/*' --permutation_glob '$(WORK)/cross-ancestry/cis-sqtl/permutations/maf001/results/cis/merged-permutations/*' --maf 0.01 --plink_glob '$(GP)/genotypes/*' -entry conditional $(ROOT)/tensorqtl-cis-sqtl.nf &

cross-ancestry-cissqtl-conditional-maf0001: GP=$(WORK)/cross-ancestry/cis-sqtl/genotypes-and-phenotypes/maf0001/results
cross-ancestry-cissqtl-conditional-maf0001: ANALYSIS=$(WORK)/cross-ancestry/cis-sqtl/conditional/maf0001
cross-ancestry-cissqtl-conditional-maf0001:
	mkdir -p $(ANALYSIS)/data
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --covariate_glob '$(DATA)/cross-ancestry/cis-sqtl/*.10.covariates.tsv' --phenotype_glob '$(GP)/phenotypes/*.phenotypes.bed.gz' --phenotype_group_glob '$(GP)/phenotype-groups/*' --permutation_glob '$(WORK)/cross-ancestry/cis-sqtl/permutations/maf0001/results/cis/merged-permutations/*' --maf 0.001 --plink_glob '$(GP)/genotypes/*' -entry conditional $(ROOT)/tensorqtl-cis-sqtl.nf &

# GWAS on cis-eQTL phenotype PCs 
cross-ancestry-gwas-on-cis-eqtl-phenotype-pcs: ANALYSIS=$(WORK)/cross-ancestry/gwas-on-cis-eqtl-phenotype-pcs
cross-ancestry-gwas-on-cis-eqtl-phenotype-pcs:
	mkdir -p $(ANALYSIS)/data
	ln -s $(DATA)/cross-ancestry/cis-eqtl/Lung.tensorqtl-in.75.covariates.tsv $(ANALYSIS)/data/
	ln -s $(DATA)/cross-ancestry/cis-eqtl/Monocyte.tensorqtl-in.30.covariates.tsv $(ANALYSIS)/data/
	ln -s $(DATA)/cross-ancestry/cis-eqtl/Nasal_epithelial.tensorqtl-in.30.covariates.tsv $(ANALYSIS)/data/
	ln -s $(DATA)/cross-ancestry/cis-eqtl/PBMC.tensorqtl-in.30.covariates.tsv $(ANALYSIS)/data/
	ln -s $(DATA)/cross-ancestry/cis-eqtl/T_cell.tensorqtl-in.30.covariates.tsv $(ANALYSIS)/data/
	ln -s $(DATA)/cross-ancestry/cis-eqtl/Whole_blood.tensorqtl-in.100.covariates.tsv $(ANALYSIS)/data/
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --covariate_glob '$(ANALYSIS)/data/*.covariates.tsv' --mac 5 --plink_glob '$(DATA)/genotypes/chr*' -entry gwas $(ROOT)/tensorqtl-gwas-on-phenotype-pcs.nf &

# GWAS on cis-sQTL phenotype PCs 
cross-ancestry-gwas-on-cis-sqtl-phenotype-pcs: ANALYSIS=$(WORK)/cross-ancestry/gwas-on-cis-sqtl-phenotype-pcs
cross-ancestry-gwas-on-cis-sqtl-phenotype-pcs:
	mkdir -p $(ANALYSIS)/data
	ln -s $(DATA)/cross-ancestry/cis-sqtl/*.tensorqtl-in.10.covariates.tsv $(ANALYSIS)/data/
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --covariate_glob '$(ANALYSIS)/data/*.covariates.tsv' --mac 5 --plink_glob '$(DATA)/genotypes/chr*' -entry gwas $(ROOT)/tensorqtl-gwas-on-phenotype-pcs.nf &

cross-ancestry-trans-eqtl-maf005: ANALYSIS=$(WORK)/cross-ancestry/trans-eqtl/maf005
cross-ancestry-trans-eqtl-maf005:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results \
					--gencode_gtf $(DATA)/crossmap/gencode-v30/gencode.v30.annotation.gtf \
					--cross_mappability $(DATA)/crossmap/gencode-v30/crossmap.txt \
					--gene_mappability $(DATA)/crossmap/gencode-v30/gene_mappability.txt \
					--snp_mappability $(DATA)/crossmap/gencode-v30/snp_mappability_100mer_2mismatch.bed.gz \
					--snp_mappability_index $(DATA)/crossmap/gencode-v30/snp_mappability_100mer_2mismatch.bed.gz.tbi \
					--covariate_glob '$(DATA)/cross-ancestry/cis-eqtl/*.covariates.tsv' --phenotype_glob '$(DATA)/cross-ancestry/cis-eqtl/*.phenotypes.bed.gz' \
					--maf 0.05 --maf_finemapping 0.05 --plink_glob '$(DATA)/genotypes/chr*' -entry trans $(ROOT)/tensorqtl-trans-eqtl.nf &


cross-ancestry-trans-sqtl-maf005: ANALYSIS=$(WORK)/cross-ancestry/trans-sqtl/maf005
cross-ancestry-trans-sqtl-maf005:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results \
					--gencode_gtf $(DATA)/crossmap/gencode-v30/gencode.v30.annotation.gtf \
					--cross_mappability $(DATA)/crossmap/gencode-v30/crossmap.txt \
					--gene_mappability $(DATA)/crossmap/gencode-v30/gene_mappability.txt \
					--snp_mappability $(DATA)/crossmap/gencode-v30/snp_mappability_100mer_2mismatch.bed.gz \
					--snp_mappability_index $(DATA)/crossmap/gencode-v30/snp_mappability_100mer_2mismatch.bed.gz.tbi \
					--phenotype_groups_glob '$(DATA)/cross-ancestry/cis-sqtl/*.leafcutter.phenotype_groups.txt' \
					--covariate_glob '$(DATA)/cross-ancestry/cis-sqtl/*.covariates.tsv' --phenotype_glob '$(DATA)/cross-ancestry/cis-sqtl/*.phenotypes.bed.gz' \
					--maf 0.05 --plink_glob '$(DATA)/genotypes/chr*' -entry trans $(ROOT)/tensorqtl-trans-sqtl.nf &


# saturation analysis
saturation-ciseqtl-genotypes-and-phenotypes-maf001: ANALYSIS=$(WORK)/saturation/cis-eqtl/genotypes-and-phenotypes/maf001
saturation-ciseqtl-genotypes-and-phenotypes-maf001:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --phenotype_glob '$(DATA)/saturation/eqtl/*.phenotypes.bed.gz' --maf 0.01 --plink_glob '$(DATA)/genotypes/chr*' --exclude_variants $(DATA)/variants-exclude.txt -entry ciseqtl $(ROOT)/prep-genotypes.nf &

saturation-ciseqtl-permutations-maf001: GP=$(WORK)/saturation/cis-eqtl/genotypes-and-phenotypes/maf001/results
saturation-ciseqtl-permutations-maf001: ANALYSIS=$(WORK)/saturation/cis-eqtl/permutations/maf001
saturation-ciseqtl-permutations-maf001:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --covariate_glob '$(DATA)/saturation/eqtl/*.covariates.tsv' --phenotype_glob '$(GP)/phenotypes/*.phenotypes.bed.gz' --maf 0.01 --plink_glob '$(GP)/genotypes/*' -entry permutations $(ROOT)/saturation-cis-eqtl.nf &

saturation-ciseqtl-susie-maf001: GP=$(WORK)/saturation/cis-eqtl/genotypes-and-phenotypes/maf001/results
saturation-ciseqtl-susie-maf001: ANALYSIS=$(WORK)/saturation/cis-eqtl/susie/maf001
saturation-ciseqtl-susie-maf001:
	mkdir -p $(ANALYSIS)/permutations
	cp $(WORK)/saturation/cis-eqtl/permutations/maf001/results/permutations/merge-chroms/*.cis_qtl.txt.gz $(ANALYSIS)/permutations/
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --permutation_glob '$(ANALYSIS)/permutations/*' --covariate_glob '$(DATA)/saturation/eqtl/*.covariates.tsv' --phenotype_glob '$(GP)/phenotypes/*.phenotypes.bed.gz' --maf 0.01 --plink_glob '$(GP)/genotypes/*' -entry susie $(ROOT)/saturation-cis-eqtl.nf &

saturation-gwas-on-cis-eqtl-phenotype-pcs: ANALYSIS=$(WORK)/saturation/gwas-on-cis-eqtl-phenotype-pcs
saturation-gwas-on-cis-eqtl-phenotype-pcs:
	mkdir -p $(ANALYSIS)/data
	ln -s $(DATA)/saturation/eqtl/*.tensorqtl-in.100.covariates.tsv $(ANALYSIS)/data/
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --covariate_glob '$(ANALYSIS)/data/*.covariates.tsv' --mac 5 --plink_glob '$(DATA)/genotypes/chr*' -entry gwas $(ROOT)/tensorqtl-gwas-on-phenotype-pcs.nf &

saturation-trans-eqtl-maf005: ANALYSIS=$(WORK)/saturation/trans-eqtl/maf005
saturation-trans-eqtl-maf005:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results \
					--gencode_gtf $(DATA)/crossmap/gencode-v30/gencode.v30.annotation.gtf \
					--cross_mappability $(DATA)/crossmap/gencode-v30/crossmap.txt \
					--gene_mappability $(DATA)/crossmap/gencode-v30/gene_mappability.txt \
					--snp_mappability $(DATA)/crossmap/gencode-v30/snp_mappability_100mer_2mismatch.bed.gz \
					--snp_mappability_index $(DATA)/crossmap/gencode-v30/snp_mappability_100mer_2mismatch.bed.gz.tbi \
					--covariate_glob '$(DATA)/saturation/eqtl/*.covariates.tsv' --phenotype_glob '$(DATA)/saturation/eqtl/*.phenotypes.bed.gz' \
					--maf 0.05 --maf_finemapping 0.05 --plink_glob '$(DATA)/genotypes/chr*' -entry trans $(ROOT)/saturation-trans-eqtl.nf &



saturation-cissqtl-genotypes-and-phenotypes-maf001: ANALYSIS=$(WORK)/saturation/cis-sqtl/genotypes-and-phenotypes/maf001
saturation-cissqtl-genotypes-and-phenotypes-maf001:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --phenotype_glob '$(DATA)/saturation/sqtl/*.phenotypes.bed.gz' --phenotype_group_glob '$(DATA)/saturation/sqtl/*.phenotype_groups.txt' --maf 0.01 --plink_glob '$(DATA)/genotypes/chr*' --exclude_variants $(DATA)/variants-exclude.txt -entry cissqtl $(ROOT)/prep-genotypes.nf &

saturation-cissqtl-permutations-maf001: GP=$(WORK)/saturation/cis-sqtl/genotypes-and-phenotypes/maf001/results
saturation-cissqtl-permutations-maf001: ANALYSIS=$(WORK)/saturation/cis-sqtl/permutations/maf001
saturation-cissqtl-permutations-maf001:
	mkdir -p $(ANALYSIS)/data
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --covariate_glob '$(DATA)/saturation/sqtl/*.covariates.tsv' --phenotype_glob '$(GP)/phenotypes/*.phenotypes.bed.gz' --phenotype_group_glob '$(GP)/phenotype-groups/*' --maf 0.01 --plink_glob '$(GP)/genotypes/*' -entry permutations $(ROOT)/saturation-cis-sqtl.nf &

saturation-cissqtl-nominal-maf001: GP=$(WORK)/saturation/cis-sqtl/genotypes-and-phenotypes/maf001/results
saturation-cissqtl-nominal-maf001: ANALYSIS=$(WORK)/saturation/cis-sqtl/nominal/maf001
saturation-cissqtl-nominal-maf001:
	mkdir -p $(ANALYSIS)/data
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --covariate_glob '$(DATA)/saturation/sqtl/*.10.covariates.tsv' --phenotype_glob '$(GP)/phenotypes/*.phenotypes.bed.gz' --phenotype_group_glob '$(GP)/phenotype-groups/*' --permutation_glob '$(WORK)/saturation/sqtl/permutations/maf001/results/cis/merged-permutations/*' --maf 0.01 --plink_glob '$(GP)/genotypes/*' -entry nominal $(ROOT)/saturation-cis-sqtl.nf &



# ancestry-specific analyses
ancestry-specific-ciseqtl-genotypes-and-phenotypes: ANALYSIS=$(WORK)/ancestry-specific/cis-eqtl/genotypes-and-phenotypes
ancestry-specific-ciseqtl-genotypes-and-phenotypes:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --phenotype_glob '$(DATA)/ancestry-specific/eqtl/*.phenotypes.bed.gz' --maf 0.01 --plink_glob '$(DATA)/genotypes/chr*' --exclude_variants $(DATA)/variants-exclude.txt -entry ciseqtl $(ROOT)/prep-genotypes.nf &

ancestry-specific-ciseqtl-permutations: GP=$(WORK)/ancestry-specific/cis-eqtl/genotypes-and-phenotypes/results
ancestry-specific-ciseqtl-permutations: ANALYSIS=$(WORK)/ancestry-specific/cis-eqtl/permutations
ancestry-specific-ciseqtl-permutations:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --covariate_glob '$(DATA)/ancestry-specific/eqtl/*.covariates.tsv' --phenotype_glob '$(GP)/phenotypes/*.phenotypes.bed.gz' --maf 0.01 --plink_glob '$(GP)/genotypes/*' -entry permutations $(ROOT)/tensorqtl-cis-eqtl.nf &

ancestry-specific-ciseqtl-nominal: GP=$(WORK)/ancestry-specific/cis-eqtl/genotypes-and-phenotypes/results
ancestry-specific-ciseqtl-nominal: ANALYSIS=$(WORK)/ancestry-specific/cis-eqtl/nominal
ancestry-specific-ciseqtl-nominal:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --covariate_glob '$(DATA)/ancestry-specific/eqtl/*.covariates.tsv' --phenotype_glob '$(GP)/phenotypes/*.phenotypes.bed.gz' --maf 0.01 --plink_glob '$(GP)/genotypes/*' -entry nominal $(ROOT)/tensorqtl-cis-eqtl.nf &

ancestry-specific-ciseqtl-susie: GP=$(WORK)/ancestry-specific/cis-eqtl/genotypes-and-phenotypes/results
ancestry-specific-ciseqtl-susie: PERMUTATIONS=$(WORK)/ancestry-specific/cis-eqtl/permutations/results/permutations/merge-chroms
ancestry-specific-ciseqtl-susie: ANALYSIS=$(WORK)/ancestry-specific/cis-eqtl/susie
ancestry-specific-ciseqtl-susie:
	mkdir -p $(ANALYSIS)/permutations
	$(foreach x,Lung___EUR.20 Monocyte___EUR.15 Nasal_epithelial___EUR.20 PBMC___AFR.20 PBMC___EAS.10 PBMC___EUR.20 T_cell___EUR.15 Whole_blood___AFR.40 Whole_blood___EUR.40,cp $(PERMUTATIONS)/$(x).cis_qtl.txt.gz $(ANALYSIS)/permutations/$(NL))
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --permutation_glob '$(ANALYSIS)/permutations/*' --covariate_glob '$(DATA)/ancestry-specific/eqtl/*.covariates.tsv' --phenotype_glob '$(GP)/phenotypes/*.phenotypes.bed.gz' --maf 0.01 --plink_glob '$(GP)/genotypes/*' -entry susie $(ROOT)/tensorqtl-cis-eqtl.nf &



ancestry-specific-cissqtl-genotypes-and-phenotypes: ANALYSIS=$(WORK)/ancestry-specific/cis-sqtl/genotypes-and-phenotypes
ancestry-specific-cissqtl-genotypes-and-phenotypes:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --phenotype_glob '$(DATA)/ancestry-specific/sqtl/*.phenotypes.bed.gz' --phenotype_group_glob '$(DATA)/ancestry-specific/sqtl/*.phenotype_groups.txt' --maf 0.01 --plink_glob '$(DATA)/genotypes/chr*' --exclude_variants $(DATA)/variants-exclude.txt -entry cissqtl $(ROOT)/prep-genotypes.nf &

ancestry-specific-cissqtl-permutations: GP=$(WORK)/ancestry-specific/cis-sqtl/genotypes-and-phenotypes/results
ancestry-specific-cissqtl-permutations: ANALYSIS=$(WORK)/ancestry-specific/cis-sqtl/permutations
ancestry-specific-cissqtl-permutations:
	mkdir -p $(ANALYSIS)/data
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --covariate_glob '$(DATA)/ancestry-specific/sqtl/*.covariates.tsv' --phenotype_glob '$(GP)/phenotypes/*.phenotypes.bed.gz' --phenotype_group_glob '$(GP)/phenotype-groups/*' --maf 0.01 --plink_glob '$(GP)/genotypes/*' -entry permutations $(ROOT)/tensorqtl-cis-sqtl.nf &

ancestry-specific-cissqtl-nominal: GP=$(WORK)/ancestry-specific/cis-sqtl/genotypes-and-phenotypes/results
ancestry-specific-cissqtl-nominal: ANALYSIS=$(WORK)/ancestry-specific/cis-sqtl/nominal
ancestry-specific-cissqtl-nominal:
	mkdir -p $(ANALYSIS)/data
	$(foreach t,Lung___EUR Nasal_epithelial___EUR PBMC___EUR Whole_blood___AFR Whole_blood___EUR,ln -s $(DATA)/ancestry-specific/sqtl/$(t).tensorqtl-in.10.covariates.tsv $(ANALYSIS)/data/$(NL))
	$(foreach t,Monocyte___EUR PBMC___AFR PBMC___EAS T_cell___EUR,ln -s $(DATA)/ancestry-specific/sqtl/$(t).tensorqtl-in.5.covariates.tsv $(ANALYSIS)/data/$(NL))
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --covariate_glob '$(ANALYSIS)/data/*.covariates.tsv' --phenotype_glob '$(GP)/phenotypes/*.phenotypes.bed.gz' --phenotype_group_glob '$(GP)/phenotype-groups/*' --permutation_glob '$(WORK)/ancestry-specific/cis-sqtl/permutations/results/cis/merged-permutations/*' --maf 0.01 --plink_glob '$(GP)/genotypes/*' -entry nominal $(ROOT)/tensorqtl-cis-sqtl.nf &

ancestry-specific-cissqtl-susie: GP=$(WORK)/ancestry-specific/cis-sqtl/genotypes-and-phenotypes/results
ancestry-specific-cissqtl-susie: ANALYSIS=$(WORK)/ancestry-specific/cis-sqtl/susie
ancestry-specific-cissqtl-susie:
	mkdir -p $(ANALYSIS)/data/sig-pairs
	mkdir -p $(ANALYSIS)/data/covariates
	$(foreach t,Lung___EUR Nasal_epithelial___EUR PBMC___EUR Whole_blood___AFR Whole_blood___EUR,ln -s $(DATA)/ancestry-specific/sqtl/$(t).tensorqtl-in.10.covariates.tsv $(ANALYSIS)/data/covariates/$(NL))
	$(foreach t,Monocyte___EUR PBMC___AFR PBMC___EAS T_cell___EUR,ln -s $(DATA)/ancestry-specific/sqtl/$(t).tensorqtl-in.5.covariates.tsv $(ANALYSIS)/data/covariates/$(NL))
	$(foreach t,Lung___EUR Nasal_epithelial___EUR PBMC___EUR Whole_blood___AFR Whole_blood___EUR Monocyte___EUR PBMC___AFR PBMC___EAS T_cell___EUR,ln -s $(WORK)/ancestry-specific/cis-sqtl/nominal/results/cis/nominal/$(t).cis_qtl.signif_pairs.parquet $(ANALYSIS)/data/sig-pairs/$(NL))
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --signif_pairs_glob '$(ANALYSIS)/data/sig-pairs/*.parquet' --covariate_glob '$(ANALYSIS)/data/covariates/*.covariates.tsv' --phenotype_glob '$(GP)/phenotypes/*.phenotypes.bed.gz' --maf 0.01 --plink_glob '$(GP)/genotypes/*' -entry susie $(ROOT)/tensorqtl-cis-sqtl.nf &

ancestry-specific-gwas-on-cis-eqtl-phenotype-pcs: ANALYSIS=$(WORK)/ancestry-specific/gwas-on-cis-eqtl-phenotype-pcs
ancestry-specific-gwas-on-cis-eqtl-phenotype-pcs:
	mkdir -p $(ANALYSIS)/data
	ln -s $(DATA)/ancestry-specific/eqtl/Lung___EUR.tensorqtl-in.20.covariates.tsv $(ANALYSIS)/data/
	ln -s $(DATA)/ancestry-specific/eqtl/Monocyte___EUR.tensorqtl-in.15.covariates.tsv $(ANALYSIS)/data/
	ln -s $(DATA)/ancestry-specific/eqtl/Nasal_epithelial___EUR.tensorqtl-in.20.covariates.tsv $(ANALYSIS)/data/
	ln -s $(DATA)/ancestry-specific/eqtl/PBMC___AFR.tensorqtl-in.20.covariates.tsv $(ANALYSIS)/data/
	ln -s $(DATA)/ancestry-specific/eqtl/PBMC___EAS.tensorqtl-in.10.covariates.tsv $(ANALYSIS)/data/
	ln -s $(DATA)/ancestry-specific/eqtl/PBMC___EUR.tensorqtl-in.20.covariates.tsv $(ANALYSIS)/data/
	ln -s $(DATA)/ancestry-specific/eqtl/T_cell___EUR.tensorqtl-in.15.covariates.tsv $(ANALYSIS)/data/
	ln -s $(DATA)/ancestry-specific/eqtl/Whole_blood___AFR.tensorqtl-in.40.covariates.tsv $(ANALYSIS)/data/
	ln -s $(DATA)/ancestry-specific/eqtl/Whole_blood___EUR.tensorqtl-in.40.covariates.tsv $(ANALYSIS)/data/
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --covariate_glob '$(ANALYSIS)/data/*.covariates.tsv' --mac 5 --plink_glob '$(DATA)/genotypes/chr*' -entry gwas $(ROOT)/tensorqtl-gwas-on-phenotype-pcs.nf &

# GWAS on cis-sQTL phenotype PCs 
ancestry-specific-gwas-on-cis-sqtl-phenotype-pcs: ANALYSIS=$(WORK)/ancestry-specific/gwas-on-cis-sqtl-phenotype-pcs
ancestry-specific-gwas-on-cis-sqtl-phenotype-pcs:
	mkdir -p $(ANALYSIS)/data
	$(foreach t,Lung___EUR Nasal_epithelial___EUR PBMC___EUR Whole_blood___AFR Whole_blood___EUR,ln -s $(DATA)/ancestry-specific/sqtl/$(t).tensorqtl-in.10.covariates.tsv $(ANALYSIS)/data/$(NL))
	$(foreach t,Monocyte___EUR PBMC___AFR PBMC___EAS T_cell___EUR,ln -s $(DATA)/ancestry-specific/sqtl/$(t).tensorqtl-in.5.covariates.tsv $(ANALYSIS)/data/$(NL))
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --covariate_glob '$(ANALYSIS)/data/*.covariates.tsv' --mac 5 --plink_glob '$(DATA)/genotypes/chr*' -entry gwas $(ROOT)/tensorqtl-gwas-on-phenotype-pcs.nf &


ancestry-specific-trans-eqtl-maf005: ANALYSIS=$(WORK)/ancestry-specific/trans-eqtl/maf005
ancestry-specific-trans-eqtl-maf005:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results \
					--gencode_gtf $(DATA)/crossmap/gencode-v30/gencode.v30.annotation.gtf \
					--cross_mappability $(DATA)/crossmap/gencode-v30/crossmap.txt \
					--gene_mappability $(DATA)/crossmap/gencode-v30/gene_mappability.txt \
					--snp_mappability $(DATA)/crossmap/gencode-v30/snp_mappability_100mer_2mismatch.bed.gz \
					--snp_mappability_index $(DATA)/crossmap/gencode-v30/snp_mappability_100mer_2mismatch.bed.gz.tbi \
					--covariate_glob '$(DATA)/ancestry-specific/eqtl/*.covariates.tsv' --phenotype_glob '$(DATA)/ancestry-specific/eqtl/*.phenotypes.bed.gz' \
					--maf 0.05 --maf_finemapping 0.05 --plink_glob '$(DATA)/genotypes/chr*' -entry trans $(ROOT)/tensorqtl-trans-eqtl.nf &


ancestry-specific-trans-sqtl-maf005: ANALYSIS=$(WORK)/ancestry-specific/trans-sqtl/maf005
ancestry-specific-trans-sqtl-maf005:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results \
					--gencode_gtf $(DATA)/crossmap/gencode-v30/gencode.v30.annotation.gtf \
					--cross_mappability $(DATA)/crossmap/gencode-v30/crossmap.txt \
					--gene_mappability $(DATA)/crossmap/gencode-v30/gene_mappability.txt \
					--snp_mappability $(DATA)/crossmap/gencode-v30/snp_mappability_100mer_2mismatch.bed.gz \
					--snp_mappability_index $(DATA)/crossmap/gencode-v30/snp_mappability_100mer_2mismatch.bed.gz.tbi \
					--phenotype_groups_glob '$(DATA)/ancestry-specific/sqtl/*.leafcutter.phenotype_groups.txt' \
					--covariate_glob '$(DATA)/ancestry-specific/sqtl/*.covariates.tsv' --phenotype_glob '$(DATA)/ancestry-specific/sqtl/*.phenotypes.bed.gz' \
					--maf 0.05 --plink_glob '$(DATA)/genotypes/chr*' -entry trans $(ROOT)/tensorqtl-trans-sqtl.nf &





