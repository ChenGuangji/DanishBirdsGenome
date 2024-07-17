#!/usr/bin/perl -w
use strict;
use File::Basename qw(basename dirname);

die "Usage:perl $0 GFF.file\n" if(@ARGV!=1);

my $GFF_file=shift;
my $name = basename($GFF_file);
my %Gene;
open (IN,$GFF_file) || die "$!\n";
while(<IN>){
	chomp;
	if(/^#/){next;}
	my @c=split /\t/;
	@c[3,4]=@c[4,3] if($c[3]>$c[4]);
	if ($c[2] eq "mRNA"){
		my $id=$1 if($c[8]=~/ID=(\S+?);/ || $c[8]=~/GenePrediction\s+(\S+)/);
		@{$Gene{$id}{mRNA}}=@c;
	}
	if($c[2] eq "CDS"){
		my $id=$1 if($c[8]=~/Parent=(\S+?);/ || $c[8]=~/GenePrediction\s+(\S+)/);
		push @{$Gene{$1}{CDS}},[@c];
	}
}
close IN;

#print "Gene\tGene_Len\tCDS_Len\tStrain\tScore\tExon_Num\tExon_Len\tIntron_Len\n";

my ($gene_n,$geneL_n,$cdsL_n,$exon_n,$exonL_n,$intro_n,$introL_n,$single_exon);
foreach my $id (keys %Gene){
	my $mRNAlen= $Gene{$id}{mRNA}[4]-$Gene{$id}{mRNA}[3]+1;
	my $score= $Gene{$id}{mRNA}[5];
	my $string= $Gene{$id}{mRNA}[6];
	my $cdslen=0;
	my @exonlen;
	my @intronlen;
	@{$Gene{$id}{CDS}}= sort {$a->[3] <=> $b->[3]} @{$Gene{$id}{CDS}};
	if (@{$Gene{$id}{CDS}} == 1){$single_exon++};
	for(my $i=0; $i< @{$Gene{$id}{CDS}}; $i++){
		my $elen= $Gene{$id}{CDS}[$i][4]-$Gene{$id}{CDS}[$i][3]+1;
		$cdslen= $cdslen+$elen;
		push @exonlen, $elen;
		$exon_n++;
		$exonL_n+=$elen;
		next if($i==0);
		my $ilen= $Gene{$id}{CDS}[$i][3]-$Gene{$id}{CDS}[$i-1][4]-1;
		push @intronlen, $ilen;
		$intro_n++;
		$introL_n+=$ilen;
	}
	#print "$id\t$mRNAlen\t$cdslen\t$string\t$score\t", scalar @exonlen,"\t";
#	print "$id\t$mRNAlen\t$cdslen\t", scalar @exonlen,"\n";
	#print join (",",@exonlen), "\t";
	#print join (",",@intronlen), "\n";
	$gene_n++;
	$geneL_n+=$mRNAlen;
	$cdsL_n+=$cdslen;
}
my $avg_geneL=int($geneL_n/$gene_n);
my $avg_cdsL=int($cdsL_n/$gene_n);
my $avg_exon=int($exon_n/$gene_n);
my $avg_exonL=int($exonL_n/$exon_n);
my $avg_intro=int($intro_n/$gene_n);
my $avg_introL=int($introL_n/$intro_n);
print "$name\t$gene_n\t$avg_geneL\t$avg_cdsL\t$avg_exon\t$avg_exonL\t$avg_intro\t$avg_introL\t$single_exon\n";
