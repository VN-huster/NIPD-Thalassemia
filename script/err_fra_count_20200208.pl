#!/usr/bin/perl
use strict;
die "perl $0 <vcf> <prefix>\n" if @ARGV<2;
my $vcf = shift;
my $prefix  = shift;
open VCF, "zcat -f $vcf|"||die $!;
my ($mat_id, $pat_id);
my $plasma_id;
my (@count, @count_het);
while(<VCF>){
    chomp;
    next if /^##/;
    my ($chr, $pos, $id , $ref, $alt ,$qual, $filter, $info, $format,@sample_ID)=split;
	for (my $i=0;$i<=$#sample_ID;$i++){
		$mat_id=$i if $sample_ID[$i]=~m/mat/i;
		$pat_id=$i if $sample_ID[$i]=~m/pat/i;
		$plasma_id=$i if $sample_ID[$i]=~m/plasma/i;
	}
    last;
}
close VCF;
my $snpa;
my $snpb;
my $ff_total;
my ($err, $n_for_err, $fra, $n_for_fra);
open FRA, "> $prefix.fra" or die $!;

open VCF, "tabix -B $vcf /THL4/home/bgi_guofengyu/work/NIPT/HPC/config/npitv3-2P1_capture_targets.bed|"||die $!;
while(<VCF>){
    chomp;
    s/\|/\//g;
    s/1\/0/0\/1/g;
    my ($chr, $pos, $id , $ref, $alt ,$qual, $filter, $info, $format, @sample)=split;
	$chr=~s/^chr//;
	next unless grep {$_ eq $chr} 1..22;
	$filter eq "PASS" or $filter eq "." or next;
	my $mat=$sample[$mat_id];
	my $pat=$sample[$pat_id];
	my $plasma=$sample[$plasma_id];
	next if $mat eq "./." or $pat eq "./." or $plasma eq "./.";
	next if $mat eq "." or $pat eq "." or $plasma eq ".";
	my $mat_gt=(split /:/, $mat)[0];
	my $pat_gt=(split /:/, $pat)[0];
	my @f=split /:/, $format;
	my %hash_format;
	@hash_format{@f}=0..$#f;
	if(!exists $hash_format{'AD'}){
		next;
	}
	my ($plasma_gt, $plasma_AD)=(split /:/, $plasma)[$hash_format{'GT'}, $hash_format{'AD'}];
	my ($ad0, $ad1)=split /,/, $plasma_AD;
	my $dp=$ad0+$ad1;
	next if $dp<30;
	if ($pat_gt eq "0/0" && $mat_gt eq "0/0"){
		$snpa++;
		$err+=$ad1;
		$n_for_err+=$ad0+$ad1;
	}elsif ($pat_gt eq "1/1" && $mat_gt eq "1/1"){
		$snpa++;
		$err+=$ad0;
		$n_for_err+=$ad0+$ad1;
	}elsif($pat_gt eq "0/0" && $mat_gt eq "1/1"){
		$snpb++;
		$fra+=$ad0;
		$n_for_fra+=$ad0+$ad1;
		$ff_total+=2*$ad0/$dp;
		print FRA "$chr\t$pos\t$ref\t$alt\t$ad0\t$dp\t" . (2*$ad0/$dp)."\n";
	}elsif($pat_gt eq "1/1" && $mat_gt eq "0/0"){
		$snpb++;
		$fra+=$ad1;
		$n_for_fra+=$ad0+$ad1;
		$ff_total+=2*$ad1/$dp;
		print FRA "$chr\t$pos\t$ref\t$alt\t$ad1\t$dp\t" . (2*$ad1/$dp)."\n";
	}
}
my $err_rate=$err/$n_for_err;
my $fetal_fra1=$fra*2/$n_for_fra;
my $fetal_fra2=$ff_total/$snpb;
print "$prefix\t$snpa\t$snpb\t$err_rate\t$fetal_fra1\t$fetal_fra2\n";
#print "sample\tsnp\thet_snp\n";
