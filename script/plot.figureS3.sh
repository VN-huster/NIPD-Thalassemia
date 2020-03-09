mkdir -p workdir/output/nipt_consistency_plot
cd workdir/output/nipt_consistency_plot/
ln -sf ../F*_nipt.OR*out .
ln -sf ../../../script/*relation.list.sort ./
Rscript ../../../script/alpha_nipt_consistency.R 
Rscript ../../../script/beta_nipt_consistency.R 
