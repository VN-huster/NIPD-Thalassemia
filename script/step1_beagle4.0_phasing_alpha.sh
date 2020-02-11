#!/bin/bash
date
script=`basename $0`
if [ $# -lt 4 ]
then
    echo -e "\033[1;32m$script\033[0;0m \033[34m<fam_vcf> <population_vcf> <prefix> <cnv_info>\033[0m"
    exit 1
fi
		
fam_vcf=$1
popu_vcf=$2
prefix=$3
cnv=$4
#perl /THL4/home/bgi_guofengyu/work/haplotyping/script/aratio_filter.pl /THL4/home/bgi_guofengyu/work/haplotyping/freq/selected_0.5M.chr11.freq.txt $fam_vcf ${prefix}.filter.vcf.gz
perl script/aratio_filter.pl /THL4/home/bgi_guofengyu/work/haplotyping/freq/selected_1M.chr16.freq.txt $fam_vcf ${prefix}.filter.vcf.gz
echo perl /THL4/home/bgi_guofengyu/work/haplotyping/script/noninvasive_analysis/get_parent_vcf.pl ${prefix}.filter.vcf.gz ${prefix}.fam.filter.vcf.gz
perl /THL4/home/bgi_guofengyu/work/haplotyping/script/noninvasive_analysis/get_parent_vcf.pl ${prefix}.filter.vcf.gz ${prefix}.fam.filter.vcf.gz
echo perl /vol6/home/bgi_zhuyaping/project/hap_analysis/50fam/NIPD_analysis/alpha.pl $cnv ${prefix}.fam.filter.vcf.gz ${prefix}.fam.filter.alpha.headers ${prefix}.fam.filter.alpha.vcf
perl /vol6/home/bgi_zhuyaping/project/hap_analysis/50fam/NIPD_analysis/alpha.pl $cnv ${prefix}.fam.filter.vcf.gz ${prefix}.fam.filter.alpha.headers ${prefix}.fam.filter.alpha.vcf
cat ${prefix}.fam.filter.alpha.headers ${prefix}.fam.filter.alpha.vcf |bgzip -c >${prefix}.fam.filter.alpha.vcf.gz
tabix -f -p vcf ${prefix}.fam.filter.alpha.vcf.gz
vcf-merge $popu_vcf ${prefix}.fam.filter.alpha.vcf.gz 2>${prefix}.fam.merge.log |bgzip -c > ${prefix}.fam.merge.vcf.gz
if [ $? -ne 0 ]; then
    echo "vcf-merge failed"
    exit 1
fi
perl /THL4/home/bgi_guofengyu/work/haplotyping/script/nomissing.pl ${prefix}.fam.merge.vcf.gz ${prefix}.fam.merge.nomissing.vcf.gz
#/vol6/home/bgi_thmed/TJ_CLUSTER/jdk1.8.0_121/bin/java -Xmx9000m -jar /THL4/home/bgi_guofengyu/work/hap_software/beagle.21Jan17.6cc.jar gt=${prefix}.fam.merge.nomissing.vcf.gz out=${prefix}.phased
java -Xmx9000m -jar /THL4/home/bgi_guofengyu/work/hap_software/beagle.r1399.jar gt=${prefix}.fam.merge.nomissing.vcf.gz usephase=true out=${prefix}.phased
if [ $? -ne 0 ]; then
    echo "beagle failed"
    exit 1
fi
date
exit 0
