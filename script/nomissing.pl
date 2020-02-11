use strict;
die "perl $0 <vcf> <outvcf>\n" if @ARGV<2;
my $vcf = shift;
my $outfile = shift;
open VCF, "zcat -f $vcf|" or die $!;
open OUT, "|bgzip -c > $outfile" or die $!;
while(<VCF>){
    chomp;
    if (/^#/){
	print OUT $_ . "\n";
	next;
    }else{
	my @t = split;
	next if $t[4] =~ m/,/;
	my $tag = 0;
	my $missing;
	for my $i (9..$#t){
	    if ($t[$i] eq "."){
		$tag = 1;
		next;
	    }
	    if ($t[$i] eq "./."){
		$missing++;
	    }
	}
	$tag++ if $missing>5;
	if ($tag){
	    next;
	}
	print OUT $_ . "\n";
    }
}
close VCF;
close OUT;
