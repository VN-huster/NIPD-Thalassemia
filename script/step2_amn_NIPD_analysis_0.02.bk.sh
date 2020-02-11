#!/bin/bash
date
script=`basename $0`
if [ $# -lt 5 ]
then
    echo -e "\033[1;32m$script\033[0;0m \033[34m<phased_vcf> <fam_vcf> <mother_mutation> <father_mutation> <prefix>\033[0m"
    exit 1
fi

phased_vcf=$1
fam_vcf=$2
mother_mutation=$3
father_mutation=$4
prefix=$5

perl /THL4/home/bgi_guofengyu/work/haplotyping/script/noninvasive_analysis/get_parent_hap_v2.pl $phased_vcf $mother_mutation $father_mutation $prefix
perl /THL4/home/bgi_guofengyu/work/haplotyping/script/amn_noninvasive_analysis/amn_famvcf_analysis.pl $fam_vcf $prefix
Rscript script/HMM_dbi_pat_0.02.R $prefix
if [ $? -ne 0 ]; then
    echo "HMM_pat failed"
    exit 1
fi
perl /THL4/home/bgi_guofengyu/work/haplotyping/script/noninvasive_analysis/pat_origin_hap.pl $prefix
#Rscript /THL4/home/bgi_guofengyu/work/haplotyping/script/noninvasive_analysis/HMM_dbi_mat.R $prefix
Rscript /vol6/home/bgi_zhuyaping/project/hap_analysis/20181201_paper/beagle_refhalf_previous/HMM_dbi_mat.R $prefix
if [ $? -ne 0 ]; then
    echo "HMM_mat failed"
    exit 1
fi
perl /THL4/home/bgi_guofengyu/work/haplotyping/script/amn_noninvasive_analysis/amn_split_v2_pat.pl $prefix
perl /THL4/home/bgi_guofengyu/work/haplotyping/script/amn_noninvasive_analysis/amn_split_v2_mat.pl $prefix
date
exit 0
