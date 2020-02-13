#!/bin/bash
date
script=`basename $0`
if [ $# -lt 3 ]
then
    echo -e "\033[1;32m$script\033[0;0m \033[34m<fam_vcf> <population_vcf> <prefix>\033[0m"
    exit 1
fi
		
fam_vcf=$1
popu_vcf=$2
prefix=$3
perl script/aratio_filter.pl db/selected_0.5M.chr11.freq.txt $fam_vcf ${prefix}.filter.vcf.gz
perl script/get_parent_vcf.pl ${prefix}.filter.vcf.gz ${prefix}.fam.filter.vcf.gz
tabix -f -p vcf ${prefix}.fam.filter.vcf.gz
vcf-merge $popu_vcf ${prefix}.fam.filter.vcf.gz 2>${prefix}.fam.merge.log |bgzip -c > ${prefix}.fam.merge.vcf.gz
if [ $? -ne 0 ]; then
    echo "vcf-merge failed"
    exit 1
fi
perl script/nomissing.pl ${prefix}.fam.merge.vcf.gz ${prefix}.fam.merge.nomissing.vcf.gz
java -Xmx9000m -jar bin/beagle.r1399.jar gt=${prefix}.fam.merge.nomissing.vcf.gz usephase=true out=${prefix}.phased
if [ $? -ne 0 ]; then
    echo "beagle failed"
    exit 1
fi
date
exit 0
