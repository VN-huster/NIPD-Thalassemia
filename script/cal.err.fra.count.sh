for i in {01..59};do perl script/err_fra_count_20200208.pl workdir/input/F$i.vcf.gz workdir/output/F$i;done > tee workdir/output/err.fra.count.txt
