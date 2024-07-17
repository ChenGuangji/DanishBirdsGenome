### outputs: Contig_N50	Contig_L50	Contig_Num	Scaffold_N50	Scaffold_L50	Scaffold_Num	Total_Len(bp)	Len_Scaffold_longer100(bp)	Gap_Len(bp)	Gap_Ratio(%)
fa=$1
Bin=`dirname $0`
perl ${Bin}/get_scaftig.pl $fa > $fa.scaftig
${Bin}/seq_n50 $fa.scaftig | awk '{if($1~/N50/){printf $2"\t"$3"\t"}else if($1~/100bp/){printf $2"\t"}}' > $fa.stat
${Bin}/seq_n50 $fa | awk '{if($1~/N50/){printf $2"\t"$3"\t"}else if($1~/Total/){printf $NF"\t"}else if($1~/100bp/){printf $2"\t"}}' | awk '{printf $1"\t"$2"\t"$4"\t"$3"\t"}' >> $fa.stat
perl ${Bin}/stat_N.pl $fa >> $fa.stat
rm $fa.scaftig
