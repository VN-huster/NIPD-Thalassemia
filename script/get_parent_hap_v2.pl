use strict;
die "perl $0 <vcf> <mother_mutataion_pos> <father_mutation_pos> <prefix>\n" if @ARGV<2;
my $vcf = shift;
my $mother_mutation = shift;
my $father_mutation = shift;
my $prefix = shift;
my ($mat_chr, $mat_pos) = split /:/, $mother_mutation;
my ($pat_chr, $pat_pos) = split /:/, $father_mutation;
my $outfile = $prefix . ".parent.hap";
open VCF, "zcat -f $vcf|" or die $!;
open OUT, "> $outfile" or die $!;
my ($mat_id, $pat_id);
my ($h_mat, $h_pat) = (-1, -1);
while(<VCF>){
    chomp;
    /^##/ and next;
    my @t = split;
    if (/^#/){
	for my $i (9..$#t){
	    $mat_id = $i if $t[$i] =~ m/mat/i;
	    $pat_id = $i if $t[$i] =~ m/pat/i;
	}
	die "Can't find mat/pat\n" unless $mat_id && $pat_id;
	next;
    }
    my ($chr, $pos, $ref, $alt, $mat, $pat) = @t[0,1,3,4,$mat_id,$pat_id];
    my $mat_gt = (split /:/, $mat)[0];
    my $pat_gt = (split /:/, $pat)[0];
    my ($M0, $M1) = split /\|/, $mat_gt;
    my ($F0, $F1) = split /\|/, $pat_gt;
    if ($chr eq $pat_chr and $pos == $pat_pos){
	$h_pat = 1 if $F0==0; 
    }
    if ($chr eq $mat_chr and $pos == $mat_pos){
	$h_mat = 1 if $M0==0; 
    }
    print OUT "$chr\t$pos\t$ref\t$alt\t$F0\t$F1\t$M0\t$M1\n";
}
close OUT;
close VCF;

if ($h_mat==1 and $h_pat==1){
    `awk -v OFS='\t' '{print \$1,\$2,\$3,\$4,\$6,\$5,\$8,\$7}' $outfile >$outfile.new`;
    `rm -rf $outfile`;
    rename "$outfile.new", $outfile;
}elsif($h_mat==1 and $h_pat==-1){
    `awk -v OFS='\t' '{print \$1,\$2,\$3,\$4,\$5,\$6,\$8,\$7}' $outfile >$outfile.new`;
    `rm -rf $outfile`;
    rename "$outfile.new", $outfile;
}elsif($h_mat==-1 and $h_pat==1){
    `awk -v OFS='\t' '{print \$1,\$2,\$3,\$4,\$6,\$5,\$7,\$8}' $outfile >$outfile.new`;
    `rm -rf $outfile`;
    rename "$outfile.new", $outfile;
}
