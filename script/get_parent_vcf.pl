use strict;
use List::Util qw/max min/;
die "perl $0 <vcf> <outvcf> [floor-depth] [het_threshold] [hom_threshold]\n" if @ARGV<2;
my $vcf = shift;
my $outvcf = shift;
my $floor_depth = shift;
my $het_threshold = shift;
my $hom_threshold = shift;
open VCF, "zcat -f $vcf|" or die $!;
open OUT,"|bgzip -c > $outvcf" or die $!;
my ($mat_id, $pat_id);
my ($GT_id, $AD_id);
my %format;
while(<VCF>){
    if (/^##/){
	print OUT;
	next;
    }else{
	chomp;
	my @t = split;
	if(/^#/){
	    for my $i (9..$#t){
		$mat_id = $i if $t[$i] =~ m/mat/i;
		$pat_id = $i if $t[$i] =~ m/pat/i;
	    }
	    die "Cat find sample's ID\n" unless ($mat_id && $pat_id);
		#print OUT "$t[0]\t$t[1]\t$t[2]\t$t[3]\t$t[4]\t$t[5]\t$t[6]\t$t[7]\t$t[8]\t$t[$mat_id]\t$t[$pat_id]\n";
	}
	next if $t[$mat_id] eq "." or $t[$mat_id] eq "./.";
	next if $t[$pat_id] eq "." or $t[$pat_id] eq "./.";
	my @f = split /:/, $t[8];
	@format{@f} = 0..$#f;
	($GT_id, $AD_id) = ($format{"GT"}, $format{"AD"});
	my ($gt_mat, $ad_mat) = (split /:/, $t[$mat_id])[$GT_id, $AD_id];
	my ($gt_pat, $ad_pat) = (split /:/, $t[$pat_id])[$GT_id, $AD_id];
	my ($ad_mat0, $ad_mat1) = split /,/, $ad_mat;
	my ($ad_pat0, $ad_pat1) = split /,/, $ad_pat;
	my $dp_mat = $ad_mat0 + $ad_mat1;
	my $dp_pat = $ad_pat0 + $ad_pat1;
	print OUT "$t[0]\t$t[1]\t$t[2]\t$t[3]\t$t[4]\t$t[5]\t$t[6]\t$t[7]\t$t[8]\t$t[$mat_id]\t$t[$pat_id]\n";
    }
}
close VCF;
close OUT;
