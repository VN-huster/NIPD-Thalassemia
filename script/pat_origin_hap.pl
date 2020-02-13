use strict;
die "perl $0 <prefix>\n" if @ARGV<1;
my $prefix = shift;
my $hapfile = $prefix . ".OUT";
my $fout = $prefix . ".fout";
open HAPFILE, $hapfile or die $!;
open FOUT, $fout or die $!;
my @FHAP = <FOUT>;
my %header;
my (@F0_s, @F0_e);
my (@F1_s, @F1_e);
my $prior = 0;
my $prior_pos;
for my $i (0..$#FHAP){
    chomp $FHAP[$i];
    my @t = split /\t/, $FHAP[$i];
    @header{@t} = 0..$#t;
    if ($t[$header{"fetal_F"}]!=$prior && $t[$header{"fetal_F"}]==1){
	push @F0_s, $t[1];
	if ($prior){
	   push @F1_e, $prior_pos; 
	}
    }
    if ($t[$header{"fetal_F"}]!=$prior && $t[$header{"fetal_F"}]==2){
	push @F1_s, $t[1];
	if ($prior){
	   push @F0_e, $prior_pos;
	}
    }
    if ($i == $#FHAP){
	push @F0_e, $t[1] if $t[$header{"fetal_F"}]==1;
	push @F1_e, $t[1] if $t[$header{"fetal_F"}]==2;
    }
    $prior = $t[$header{"fetal_F"}];
    $prior_pos = $t[1];
}
close FOUT;
print "$#F0_s\t$#F0_e\n$#F1_s\t$#F1_e\n";
print "@F0_s\n";
print "@F0_e\n";
print "@F1_s\n";
print "@F1_e\n";

my %h;
open OUT, "> $prefix.pat_origin.hap" or die $!;
while(<HAPFILE>){
    chomp;
    my @t = split;
    if (/pos/i){
	@h{@t} = 0..$#t;
	print OUT "$_\tF_origin\tfetal_F\n";
	next;
    }
    my $F_origin;
    my $fetal_F;
    if ($t[4]==$t[5]){
	$F_origin = $t[4];
	$fetal_F = 3;
    }
    for my $i (0..$#F0_s){
	if ($t[1]>=$F0_s[$i] && $t[1]<=$F0_e[$i]){
	    $fetal_F = 1;
	    $F_origin = $t[$h{"F0"}];
	    last;
	}
    }
    for my $i (0..$#F1_s){
	if ($t[1]>=$F1_s[$i] && $t[1]<=$F1_e[$i]){
	    $fetal_F = 2;
	    $F_origin = $t[$h{"F1"}];
	    last;
	}
    }
    print OUT $_ . "\t" . $F_origin . "\t" . $fetal_F . "\n" if $fetal_F;
}
close OUT;
