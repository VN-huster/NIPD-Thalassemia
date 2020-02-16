use strict;
die "perl $0 <VCF> <prefix> [plasmavcf]\n" if @ARGV<2;
my $vcf = shift;
my $prefix = shift;
my $plasmavcf = shift;

my $floor_depth = 30;
my ($het, $hom) = (0.40, 0.1);
#my ($het, $hom) = (0.20, 0.1);
my $hap = $prefix . ".parent.hap";
my %hash_fra = ( #Mat:Pat
    "0/0:1/1" => 1,
    "1/1:0/0" => 0,
);
my %hash_amn=(
	"0/0:0/1:0/0" => 0,
	"0/0:0/1:0/1" => 1,
	"1/1:0/1:0/1" => 0,
	"1/1:0/1:1/1" => 1,
	"0/1:0/1:0/0" => 0,
	"0/1:0/1:1/1" => 1,
	#"0/1:0/1:0/1" => 0.5,
);
my %hash_mp;
if($plasmavcf){
    open MP,"< $plasmavcf" or die$!;
    while(<MP>){
        chomp;
        my($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,$sample)=split;
        $hash_mp{"$chr\t$pos\t$ref\t$alt"}=$_;
    }
    close MP;
}
open HAP, $hap or die $!;
my %hash;
while(<HAP>){
    chomp;
    my ($chr, $pos, $ref, $alt, $F0, $F1, $M0, $M1) = split;
    $chr =~ s/chr//;
    $hash{"$chr\t$pos\t$ref\t$alt"} = "$F0\t$F1\t$M0\t$M1";
}
open VCF, "zcat -f $vcf|" or die $!;
open OUT, "> $prefix.OUT" or die $!;
print OUT "chr\tpos\tref\talt\tF0\tF1\tM0\tM1\tMp_ref\tMp_alt\tamn_F\tamn_M\n";
open FRA, "> $prefix.fra" or die $!;
print FRA "chr\tpos\tref\talt\tn\tfra\n";
open MAD, "> $prefix.mad" or die $!;
print MAD "chr\tpos\tMgt\tMref\tMalt\n";
my ($mat_id, $pat_id, $plasma_id, $amn_id);
while(<VCF>){
    chomp;
    /^##/ and next;
    s/\|/\//g;
    s/1\/0/0\/1/g;
    s/chrX/23/g;
    s/chrY/24/g;
    my @t = split;
    if(/^#/){
	for my $i (9..$#t){
	    $mat_id = $i if $t[$i] =~ m/mat/i;
	    $pat_id = $i if $t[$i] =~ m/pat/i;
	    $plasma_id = $i if $t[$i] =~ m/plasma/i;
	    $amn_id = $i if $t[$i] =~ m/fetus|amn|chor/i;
	}
	die "Cat find sample's ID\n" unless ($mat_id && $pat_id && $plasma_id && $amn_id);
	next;
    }
    my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,$mat,$pat,$mp,$amn) = @t[0,1,2,3,4,5,6,7,8,$mat_id,$pat_id,$plasma_id,$amn_id];
    $chr =~ s/chr//;
    my @GT = ($ref,$alt);
    next if ($mat eq "./." or $pat eq "./." or $mp eq "./.");
    my %format;
    my @f = split /:/, $t[8];
    @format{@f} = 0..$#f;
    my ($GT_id, $AD_id) = ($format{"GT"}, $format{"AD"});
    my ($gt_mat, $ad_mat) = (split /:/, $t[$mat_id])[$GT_id, $AD_id];
    my ($gt_pat, $ad_pat) = (split /:/, $t[$pat_id])[$GT_id, $AD_id];
    my ($gt_plasma, $ad_plasma) = (split /:/, $t[$plasma_id])[$GT_id, $AD_id];
    my ($gt_amn, $ad_amn) = (split /:/, $t[$amn_id])[$GT_id, $AD_id];
    exists$hash_mp{"$chr\t$pos\t$ref\t$alt"} 
        and $mp=(split /\t/,$hash_mp{"$chr\t$pos"})[-1]
        and ($gt_plasma,$ad_plasma)=split /:/,$mp;
    my @ad_mat = split /,/, $ad_mat;
    my @ad_pat = split /,/, $ad_pat;
    my @ad_plasma = split /,/, $ad_plasma;
    my $dp_mat = $ad_mat[0] + $ad_mat[1];
    my $dp_pat = $ad_pat[0] + $ad_pat[1];
    my $dp_plasma = $ad_plasma[0] + $ad_plasma[1];
    next if $dp_mat <$floor_depth or $dp_pat <$floor_depth or $dp_plasma <$floor_depth;
    my $mat_min = (sort{$a<=>$b}@ad_mat)[0];
    my $pat_min = (sort{$a<=>$b}@ad_pat)[0];
    my $plasma_min = (sort{$a<=>$b}@ad_plasma)[0];
    if ($gt_mat eq '0/1'){
	next unless $mat_min/$dp_mat > $het;
    }else{
	next unless $mat_min/$dp_mat < $hom;
    }
    if ($gt_pat eq '0/1'){
	next unless $pat_min/$dp_pat > $het;
    }else{
	next unless $pat_min/$dp_pat < $hom;
    }
    if($chr <23 and exists $hash_fra{"$gt_mat:$gt_pat"}){
	my $k = $hash_fra{"$gt_mat:$gt_pat"};
	my $fra = 2*$ad_plasma[$k]/($ad_plasma[0]+$ad_plasma[1]);
	print FRA "$chr\t$pos\t$ref\t$alt\t$k\t$fra\n";
    }
    next unless exists $hash{"$chr\t$pos\t$ref\t$alt"};
    my ($F0, $F1, $M0, $M1) = split /\t/, $hash{"$chr\t$pos\t$ref\t$alt"};
	my ($amn_F, $amn_M);
	if(exists $hash_amn{"$gt_mat:$gt_pat:$gt_amn"}){
		$amn_F=($hash_amn{"$gt_mat:$gt_pat:$gt_amn"} eq $F0) ? 1 : 2;
	}elsif($F0==$F1){
		$amn_F=1.5;
	}else{
		$amn_F=1.5;
	}
	if(exists $hash_amn{"$gt_pat:$gt_mat:$gt_amn"}){
		$amn_M=($hash_amn{"$gt_pat:$gt_mat:$gt_amn"} eq $M0) ? 1 : 2;
	}elsif($M0==$M1){
		$amn_M=1.5;
	}else{
		$amn_M=1.5;
	}
    print OUT "$chr\t$pos\t$ref\t$alt\t$F0\t$F1\t$M0\t$M1\t$ad_plasma[0]\t$ad_plasma[1]\t$amn_F\t$amn_M\n";
    $chr<=23 and print MAD "$chr\t$pos\t$gt_mat\t$ad_mat[0]\t$ad_mat[1]\n";
}
close FRA;
close OUT;
close MAD;
