#!/usr/bin/perl -w 
use strict;
die "perl $0 <prefix>\n" if @ARGV<1;
my $prefix  = shift;

my %hashM;my %hashF;
my ($m0,$m1,$f0,$f1)=(0,0,0,0);
open MOUT,"zcat -f $prefix\_nipt.OR.mout|"or die $!;
while(<MOUT>){
    chomp;
    my @ar=split(/\t/,$_);
    if (/^chr/){
	@hashM{@ar}=0..$#ar;
	next;
    }
    if ($ar[$hashM{"fetal_M"}]==1){
	$m0++;
    }elsif($ar[$hashM{"fetal_M"}]==2){
	$m1++;
    }elsif($ar[$hashM{"fetal_F"}]==0){
	if ($ar[$hashM{"P0"}]>0.5){
	    $m0++;
	}elsif($ar[$hashM{"P0"}]<0.5){
	    $m1++;
	}
    }
}
close MOUT;
open FOUT,"zcat -f $prefix\_nipt.OR.fout|"or die $!;
while(<FOUT>){
    chomp;
    my @ar=split(/\t/,$_);
    if (/^chr/){
	@hashF{@ar}=0..$#ar;
	next;
    }
    if ($ar[$hashF{"fetal_F"}]==1){
	$f0++;
    }elsif($ar[$hashF{"fetal_F"}]==2){
	$f1++;
    }elsif($ar[$hashF{"OddRatio_F"}]==0){
	if ($ar[$hashF{"P0"}]>0.5){
	    $f0++;
	}elsif($ar[$hashF{"P0"}]<0.5){
	    $f1++;
	}
    }
}
close MOUT;

print "$prefix\t$m0\t$m1\t$f0\t$f1\n";
