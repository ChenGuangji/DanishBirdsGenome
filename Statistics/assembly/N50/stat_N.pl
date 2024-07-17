#!/usr/bin/perl -w
use strict;


my $t_num=0;
my $noN_num=0;
my $longer_100_num=0;

open (IN,"$ARGV[0]") or die $!;
while (<IN>){
	chomp;
	next if $_ =~ m/>/;
	my $seq_len = length($_);
	$t_num += $seq_len;
	(my $subseq = $_)=~ s/N//ig;
	$noN_num += length($subseq);
}

close IN;

my $gap_num = $t_num-$noN_num;
my $N_rate = 100*$gap_num/$t_num;
printf "%3s\t%3s\t%.2f\n",$t_num,$gap_num,$N_rate;

