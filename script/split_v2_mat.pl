
use 5.010;
#use feature qw(say);

$#ARGV<0 and die"$0 prefix [gender]\n";
$md=shift;
$gender=shift;
$Gene=q[db/NIPT_V2_gene.txt];
$Plot_R=q[script/fetal_analysis/plot_v2.R];

open GENE,"< $Gene" or die$!;
my(%start,%end);
while(<GENE>){
    chomp;
    my($chr,$start,$end,$gene)=(split /\t/,$_)[0,1,2,5];
    $chr=~s/chrX/23/g;
    $chr=~s/chrY/24/g;
    $chr=~s/chr//g;
    $start{$gene}=$chr*1e9+$start;
    $end{$gene}  =$chr*1e9+$end;
}


open IN,"< $md.mout" or die$!;
open MOUT,"> $md.mout.plot" or die$!;
while(<IN>){
    chomp;
    @ln=split /\t/,$_;
    $ln[0]=~s/chr//g;
    $Chr=$ln[0];
    $Pos=$ln[1];
    $P0=$ln[-3];
    $P1=$ln[-2];
    $fetal_M=$ln[-1];

    my$pos=$Chr*1e9+$Pos;
    my$Region=&getGENE($pos);
    $Region or next;
    $Region=~s/;$//g or die"$pos\t$Region\n";
    my@Region=split /;/,$Region;
    for my$region(@Region){
        $fetal_M and print MOUT "$Chr\t$Pos\t$P0\t$P1\t$fetal_M\t$region\n";
    }
}
close IN;
close MOUT;
say "Rscript $Plot_R $md.mout.plot $region $Gene";
system"Rscript $Plot_R $md.mout.plot $region $Gene";
sub getGENE(){
    my$pos=$_[0];
    my$region="";
    for my$gene(keys%start){
        my$start=$start{$gene};
        my$end=$end{$gene};
        given($pos){
            when($_ <=  $start-1e6  ){                              
            }
            when($_ >   $end+1e6    ){                              
            }
=cut
            when($_ <=  $start-1e7  ){                              
            }
            when($_ >   $end+1e7    ){                              
            }
            when($_ <=  $start-3e6  ){
                $region=$region.$gene."_US_3M-10M\t$gene;";
            }
            when($_ <=  $start-2e6  ){
                $region=$region.$gene."_US_2M-3M\t$gene;";
            }
            when($_ <=  $start-1e6  ){
                $region=$region.$gene."_US_1M-2M\t$gene;";
            }
=cut
            when($_ <=  $start-5e5  ){
                $region=$region.$gene."_US_500K-1M\t$gene;";
            }
            when($_ <=  $start      ){
                $region=$region.$gene."_US_500K\t$gene;";
            }
            when($_ <=  $end        ){
                $region=$region.$gene."\t$gene;";            
            }
            when($_ <=  $end+5e5    ){
                $region=$region.$gene."_DS_500K\t$gene;";
            }
            when($_ <=  $end+1e6    ){
                $region=$region.$gene."_DS_500K-1M\t$gene;";
            }
=cut
            when($_ <=  $end+2e6    ){
                $region=$region.$gene."_DS_1M-2M\t$gene;";
            }
            when($_ <=  $end+3e6    ){
                $region=$region.$gene."_DS_2M-3M\t$gene;";
            }
            when($_ <=  $end+1e7    ){
                $region=$region.$gene."_DS_3M-10M\t$gene;";
            }
=cut
            default                  {
                say"$pos error";
            }
        }
    }
    return $region;
}




=cut
cut -f 1,2,25,26,29,30,31,32 test3.fmout.xls |head
Chr Pos AMN_F   AMN_M   P0  P1  fetal_F fetal_M
1   5001488 1   0   0.95    0.05    1   0
1   15009937    2   0   0.05    0.95    1   0
1   16001588    1   0   0.95    0.05    1   0
1   22000026    1   0   0.95    0.05    1   0
1   23000592    1   0   0.95    0.05    1   0
1   31000290    1   0   0.95    0.05    1   0
1   36021738    1   0   0.95    0.05    1   0
1   40007013    1   0   0.95    0.05    1   0
1   41006213    1   0   0.95    0.05    1   0
chr1	35246790	35251965	" 5,176 "	+	GJB3	NM_024009.2		
chr1	45965856	45976739	" 10,884 "	+	MMACHC	NM_015506.2		
chr11	5246696	5248301	" 1,606 "	-	HBB	NM_000518.4		
chr11	112097088	112104696	" 7,609 "	+	PTS	NM_000317.2		
chr12	103232104	103311381	" 79,278 "	-	PAH	NM_000277.1		
chr13	20761602	20767114	" 5,513 "	-	GJB2	NM_004004.5		
chr13	52506805	52585630	" 78,826 "	-	ATP7B	NM_000053.3		
chr16	222846	223709	 864 	+	HBA2	NM_000517.4		
chr16	226679	227520	 842 	+	HBA1	NM_000558.3		
chr17	78075355	78093679	" 18,325 "	+	GAA	NM_000152.3		
