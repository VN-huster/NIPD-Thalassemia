#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
die "perl $0 <hapfile> <fam VCF> <mat mutation> <pat mutation><prefix>\n" if @ARGV<5;
my $hapfile=shift;
my $vcf = shift;
my $mat_mut=shift;
my $pat_mut=shift;
my $prefix = shift;
my %hash_proband = (
    "0/1:0/0:0/0" => 0,
    "0/1:0/0:0/1" => 1,
    "0/1:1/1:0/1" => 0,
    "0/1:1/1:1/1" => 1,
    "0/1:0/1:0/0" => 0,
    "0/1:0/1:1/1" => 1,
    "0/0:1/1:0/1" => 0,
    "1/1:0/1:0/1" => 1,
);
my %hash_G = ( #Mat:Pat:Proband
    "0/1:0/0" =>0,
    "0/1:1/1" =>1
);
my %hash_fra = ( #Mat:Pat
    "0/0:1/1" => 1,
    "1/1:0/0" => 0,
);
open HAP, $hapfile or die $!;
my %hash_hap;
my $type;
while(<HAP>){
    chomp;
    s/chrX/23/;
    s/chrY/24/;
    my ($chr, $pos, $ref, $alt, $Pathap1, $Pathap2, $Mathap1, $Mathap2)=split;
    $chr =~ s/chr//;
    if($chr==16){
	$type=1;
    }elsif($chr==11){
	$type=2;
    }
    $hash_hap{"$chr\t$pos\t$ref\t$alt"}="$Pathap1\t$Pathap2\t$Mathap1\t$Mathap2";
    s/chrX/23/;
    s/chrY/24/;
    $chr =~ s/chr//;
}
open VCF, "zcat -f $vcf|" or die $!;
my $outfile;
if($type==1){
	$outfile="$prefix.alpha.OUT";
}elsif($type==2){
	$outfile="$prefix.beta.OUT";
}else{
	die "Not alpha or beta\n";
}
open OUT, "> $outfile" or die $!;
print OUT "chr\tpos\tref\talt\tM0\tM1\tP0\tP1\tMathap1\tMathap2\tPathap1\tPathap2\n";
my ($mat_id, $pat_id, $proband_id);
while(<VCF>){
    chomp;
    /^##/ and next;
    if (/^#/){
	my @t = split /\t/;
	for my $i (0..$#t){
	    $mat_id = $i if $t[$i] =~ m/mat/i;
	    $pat_id = $i if $t[$i] =~ m/pat/i;
	    $proband_id = $i if $t[$i] =~ m/proband|fetus/i;
	}
    }
    die "Can't find ID\n" unless  ($mat_id && $pat_id && $proband_id);
    s/\|/\//g;
    s/chrX/23/;
    s/chrY/24/;
    my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,$mat,$pat,$proband) = (split)[0,1,2,3,4,5,6,7,8,$mat_id,$pat_id,$proband_id];
    $chr =~ s/chr//;
    my @GT = ($ref,$alt);
    next if ($mat eq "./." or $pat eq "./." or $proband eq "./.");
    next unless exists $hash_hap{"$chr\t$pos\t$ref\t$alt"};
    my ($Pathap1, $Pathap2, $Mathap1, $Mathap2) = split /\t/, $hash_hap{"$chr\t$pos\t$ref\t$alt"};
    my $gtm = (split /:/, $mat)[0]; 
    my $gtp = (split /:/, $pat)[0]; 
    my $gtpr = (split /:/, $proband)[0]; 
    my @mat_hap = split /\//, $gtm;
    my @pat_hap = split /\//, $gtp;
    my @proband_hap = split /\//, $gtpr;
    $gtm =~ s/1\/0/0\/1/;
    $gtp =~ s/1\/0/0\/1/;
    $gtpr =~ s/1\/0/0\/1/;
    if (exists $hash_proband{"$gtm:$gtp:$gtpr"}){
	my @GT_mat = split /\//, $gtm;
	my @GT_pat = split /\//, $gtp;
	my @GT_proband = split /\//, $gtpr;
	my $fam_key = $hash_proband{"$gtm:$gtp:$gtpr"};
	my $M0 = $fam_key;
	my $P0 = ($M0==$GT_proband[0])?$GT_proband[1]:$GT_proband[0];
	my $M1 = ($M0==$GT_mat[0])?$GT_mat[1]:$GT_mat[0];
	my $P1 = ($P0==$GT_pat[0])?$GT_pat[1]:$GT_pat[0];
	print OUT "$chr\t$pos\t$ref\t$alt\t$M0\t$M1\t$P0\t$P1\t$Mathap1\t$Mathap2\t$Pathap1\t$Pathap2\n";
    }elsif (exists $hash_proband{"$gtp:$gtm:$gtpr"}){
	my @GT_mat = split /\//, $gtm;
	my @GT_pat = split /\//, $gtp;
	my @GT_proband = split /\//, $gtpr;
	my $fam_key = $hash_proband{"$gtp:$gtm:$gtpr"};
	my $P0 = $fam_key;
	my $M0 = ($P0==$GT_proband[0])?$GT_proband[1]:$GT_proband[0];
	my $M1 = ($M0==$GT_mat[0])?$GT_mat[1]:$GT_mat[0];
	my $P1 = ($P0==$GT_pat[0])?$GT_pat[1]:$GT_pat[0];
	print OUT "$chr\t$pos\t$ref\t$alt\t$M0\t$M1\t$P0\t$P1\t$Mathap1\t$Mathap2\t$Pathap1\t$Pathap2\n";
    }
}
close VCF;
close OUT;
system "perl $Bin/noninvasive_Haplotype_err_compute.pl $hapfile $outfile $outfile.stat $mat_mut $pat_mut";
print "perl $Bin/noninvasive_Haplotype_err_compute.pl $hapfile $outfile $outfile.stat $mat_mut $pat_mut\n";
if($type==1){
    system "Rscript $Bin/noninvasive_fam_hap_plot_alpha.R $outfile $prefix.alpha";
}elsif($type==2){
    system "Rscript $Bin/noninvasive_fam_hap_plot_beta.R $outfile $prefix.beta";
}
