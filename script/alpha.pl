use strict;
my $samIDlist = shift;
my $vcf=shift;
my $outfile1=shift;
my $outfile=shift;
my %hash;
open LIST, $samIDlist or die "can not open $samIDlist:$!\n";
while(<LIST>){
    chomp;
    #A-THAL1-mat    .	SEA/N	SEA 0/1	del37	0/0 del42   0/0
    #A-THAL1-pat    .	SEA/N	SEA 0/1	del37	0/0 del42   0/0
    #A-THAL1-proband	.   SEA/N   SEA	0/1 del37   0/0	del42	0/0
    my @t = split;
    $hash{$t[0]} = \@t;
}
open IN, "zcat -f $vcf |" or die $!;
open OUT1, "> $outfile1" or die $!;
open OUT, "|sort -V > $outfile" or die $!;
my @ID;
my %headers;
my @h;
my ($tag1, $tag2)=(0,0);
while(<IN>){
    chomp;
    if(/^##/){
	print OUT1 "$_\n";
	next;
    }elsif(/^#/){
	my @t = split;
	@h = @t;
        print OUT1 "$_\n";
    }else{
	my @t = split;
	@headers{@t} = 0..$#t;
	if ($tag1==0){
	    my (@g_sea, @g_del37, @g_del42);
	    for my $i (9..$#t){
		die "$h[$i]\n" if !exists $hash{$h[$i]};
		push @g_sea, $hash{$h[$i]}->[4];
		push @g_del37, $hash{$h[$i]}->[6];
		push @g_del42, $hash{$h[$i]}->[8];
	    }    
	    my $str1 = join "\t", @g_sea;
	    my $str2 = join "\t", @g_del37;
	    my $str3 = join "\t", @g_del42;
	    print OUT "chr16\t215400\t.\tNN\tN\t.\t.\t.\tGT\t$str1\n";
	    print OUT "chr16\t223300\t.\tNN\tN\t.\t.\t.\tGT\t$str2\n";
	    print OUT "chr16\t219817\t.\tNN\tN\t.\t.\t.\tGT\t$str3\n";
	    $tag1++;
	    #next;
	}
        print OUT "$_\n";
    }
}
