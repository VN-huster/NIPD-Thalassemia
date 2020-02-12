mkdir -p workdir/output/haplotype_consistency_plot
cd workdir/output/haplotype_consistency_plot/
Rscript ../../../script/beta_fetus_hap_consistency.v6.R 
Rscript ../../../script/alpha_fetus_hap_consistency.v6.R 
