#!/usr/bin/env perl
#### Author: Shengdi Li ####
#### 2022 April. 19 ####
#### Extraction ATAC-seq transposon insertion frequency metric from external data (Schep et al) 
#### the bedgraph used was generated from 'pyatac ins --bam' using the mapped reads from all lin samples in the cited publication, indicating Tn5 insertion frequencies spaning the genome. 
## Author: Shengdi Li
use strict;
use Getopt::Long;
my $window_size = 50;
#### Load external data ####
GetOptions('h' => \my $display_help, 'in=s' => \my $input_meta_file, 'out=s' => \my $outfile, 'w:i'=>\ $window_size, 'bed=s' => \my $bedfile);

if($display_help)
{
	print "Command:\n\tperl extract_ATACseq_ins_freq.pl [-h] -in INPUT_META_LIST -out OUTFILE -bed BEDFILE [-w WINDOW_SIZE]\n";
	print "Function: \n\tExtract mean insertion frequencies at given list of target coordinates (mean of bases in the window defined by a window size)...\n\n";
	
	print "Usage:\n";
	print "\t-h\tPrint help info.\n";
	print "\t-in\tInput sample meta sheet.\n";
	print "\t-out\tOutput file path. \n";
	print "\t-bed\tBedgraph of ATACseq insertion rates outputted by pyatac ins.\n";
	print "\t-w {50}\tLength of window to extract mean insertion frequency.\n";
	
	exit 1;
}
if((!$input_meta_file)||(!$outfile)||(!$window_size)||(!$bedfile))
{
	print "Missing arguments.\n";
	print "Use -h for looking at usage...\n";
	exit 1;
}


open(IN,$input_meta_file);
my $header =<IN>;
open(OUT,">".$outfile);
open(INS,$bedfile);
print OUT "#id\tATACseq_ins_freq\n";
my %ins_freq;
while(my $read_line=<INS>)
{
	chomp $read_line;
	if($read_line=~/^chr/)
	{
		my ($chr,$start,$end,$value)=split("\t",$read_line);
		$ins_freq{$chr}{$start}=$value;
	}
}
print "Finish Loading ATAC seq dataset...\n";
while(my $read_line=<IN>)
{
	chomp $read_line;
	my ($id,$chr,$strand,$coord,$guide)=split("\t",$read_line);
	#TEST IF VARIANT IS IN OPEN CHROMATIN REGIONS
	my $sum_ins=0;
	my $i=0;
	for(my $j=$coord-$window_size/2;$j<$coord+$window_size/2;$j++)
	{
		$sum_ins+=$ins_freq{$chr}{$j};
		$i++;
	}
	print OUT $id,"\t",sprintf("%.3f",log($sum_ins+1)/log(2)-log($i)/log(2)),"\n";
}