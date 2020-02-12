mkdir -p workdir/output/nipt_consistency_plot
cd workdir/output/nipt_consistency_plot/
ln -sf ../F*_nipt.OR*out .
Rscript alpha_nipt_consistency.v2.R 
Rscript beta_nipt_consistency.v2.R 
