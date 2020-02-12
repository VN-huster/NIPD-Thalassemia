use strict;
die "perl $0 <hapfile> <hap_info> <stat_result> <mat_mut> <pat_mut>\n" if @ARGV<5;
my $hapfile=shift;
my $file = shift;
my $outfile = shift;
my $mat_mut=shift;
my $pat_mut=shift;
#my $famID = $file=~s/.*\///r;
#$famID =~ s/\.OUT//;
my $famID=(split(/\_/,$file))[0];
my $p_mat=0;
my $p_pat=0;
open HP,$hapfile or die $!;
while(<HP>){
    chomp;
    my @ar=split;
    #$chr\t$pos\t$ref\t$alt\t$F0\t$F1\t$M0\t$M1
    if ($ar[6]!=$ar[7]){
	$p_mat++;
    }
    if ($ar[4]!=$ar[5]){
	$p_pat++;
    }
}
close HP;

open IN, $file or die $!;
my $n;
my ($n_mat, $n_mat0, $n_mat1);
my ($n_pat, $n_pat0, $n_pat1);
my ($m1,$m2,$p1,$p2);
my $mat_pos=(split(/\:/,$mat_mut))[1];
my $pat_pos=(split(/\:/,$pat_mut))[1];
while(<IN>){
    chomp;
    /pos/i and next;
    $n++;
    my ($chr, $pos, $ref, $alt, $M0, $M1, $P0, $P1, $mat0, $mat1, $pat0, $pat1) = split;
    if ($M0 != $M1){
	$n_mat++;
	$M0 == $mat0 and $n_mat0++;
	$M0 == $mat1 and $n_mat1++;
    }
    if ($P0 != $P1){
	$n_pat++;
	$P0 == $pat0 and $n_pat0++;
	$P0 == $pat1 and $n_pat1++;
    }
    if ($pos==$mat_pos){
	$m1=($M0==1)?'4':'5';
	$m2=($mat0==1)?'8':'9';
    }
    if($pos==$pat_pos){
	$p1=($P0==1)?'6':'7';
	$p2=($pat0==1)?'10':'11';
    }
}
close IN;
#print "$mat_pos\t$pat_pos\t$m1\t$m2\t$p1\t$p2\n";
my $mat_r = ($n_mat0>$n_mat1)?$n_mat0:$n_mat1;
my $pat_r = ($n_pat0>$n_pat1)?$n_pat0:$n_pat1;
my $tag_mat="raw_methord";
if ($m1 ne ""){
    $tag_mat="compare_mat_mut_hap";
    if (($m1==4 && $m2==8)||($m1==5 && $m2==9)){
	$mat_r=$n_mat0;
    }else{
	$mat_r=$n_mat1;
    }
}
my $tag_pat="raw_methord";
if ($p1 ne ""){
    $tag_pat="compare_pat_mut_hap";
    if (($p1==6 && $p2==10)||($p1==7 && $p2==11)){
	$pat_r=$n_pat0;
    }else{
	$pat_r=$n_pat1;
    }
}
#my $mat_r = ($n_mat0>$n_mat1)?$n_mat0:$n_mat1;
#my $pat_r = ($n_pat0>$n_pat1)?$n_pat0:$n_pat1;
open OUT, "> $outfile" or die $!;
#print OUT "$famID\t$n\t$n_mat\t$mat_r\t" . $mat_r/$n_mat . "\t$n_pat\t$pat_r\t" . $pat_r/$n_pat . "\n";
print OUT "$famID\t$p_mat\t$n_mat\t$mat_r\t" . $mat_r/$n_mat . "\t$p_pat\t$n_pat\t$pat_r\t" . $pat_r/$n_pat . "\t$tag_mat\t$tag_pat\n";
close OUT;

open M,">$outfile.mout.plot" or die $!;
open F,">$outfile.fout.plot" or die $!;
print M "chr\tpos\tparent\tfetus\n";
print F "chr\tpos\tparent\tfetus\n";
open IN, $file or die $!;
while(<IN>){
    chomp;
    /pos/i and next;
    my @ar=split;
    my ($chr, $pos, $ref, $alt, $M0, $M1, $P0, $P1, $mat0, $mat1, $pat0, $pat1) = split;
    if ($M0 != $M1){
	if($m1 eq ''){
	    if ($n_mat0>$n_mat1){
		print M "$chr\t$pos\t$M0\t$mat0\n";
	    }else{
		print M "$chr\t$pos\t$M0\t$mat1\n";
	    }
	}else{
	    print M "$chr\t$pos\t$ar[$m1]\t$ar[$m2]\n";
	}
    }
    if ($P0 != $P1){
	if ($p1 eq ''){
	    if ($n_pat0>$n_pat1){
		print F "$chr\t$pos\t$P0\t$pat0\n";
	    }else{
		print F "$chr\t$pos\t$P0\t$pat1\n";
	    }
	}else{
	    print F "$chr\t$pos\t$ar[$p1]\t$ar[$p2]\n";
	}
    }
}
close IN;
close M;
close F;
