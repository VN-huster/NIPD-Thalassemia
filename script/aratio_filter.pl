use strict;
die "perl $0 <freqfile> <vcf> <outfile>\n" if @ARGV<3;

my $positionlist = shift;
my $vcffile = shift;
my $outfile = shift;
my %hash;
open POSLIST, $positionlist or die $!;
while(<POSLIST>){
    chomp;
    my ($chr, $pos, $ref, $alt, $a_num, $a_freq) = split;
    $hash{"$chr\t$pos\t$ref\t$alt"}++; 
}
close POSLIST;
open VCFFILE, "zcat -f $vcffile|" or die $!;
open OUT, "|bgzip -c > $outfile" or die $!;
while(<VCFFILE>){
    if (/^#/){
	print OUT;
	next;
    }
    my ($chr, $pos, $id, $ref, $alt) = split (/\t/, $_, 6);
    next if $alt =~ /,/;
    print OUT if $hash{"$chr\t$pos\t$ref\t$alt"};
}
close VCFFILE;
close OUT;
