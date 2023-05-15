#!/usr/bin/env perl
#### Author: Shengdi Li ####
#### 2023 Feb. 16 ####
#### extract reads that overlap the target region for remapping ####

use Getopt::Long;
use Bio::SeqIO;
use strict;

GetOptions('h' => \my $display_help, 'bam=s' => \my $bam, 'region=s' => \my $region, 'out=s' =>\my $output);
if($display_help)
{
	print "Command: \n\tperl extract_reads.pl [-h] -bam INPUT_BAM -region REGION -out OUTPUT\n\n";
	print "Function: \n\tExtract reads that overlap the target region for remapping...\n\n";
	print "Usage: \n";
	print "\t-h\tPrint help info.\n";
	print "\t-bam\tThe input bam files (separated with ; if there are multiple).\n";
	print "\t-region\tThe region to extract (e.g chr3:249310-249355).\n";
	print "\t-out\tOutput file.\n";
	exit 1;
}

if((!$bam)||(!$region)||(!$output))
{
	print "Missing arguments.\n";
	print "Use -h for looking at usage...\n";
	exit 1;
}
sub min(){
	if($_[0]<$_[1])
	{
		return $_[0]
	}
	else
	{
		return $_[1]
	}
}
my ($chr,$start,$end) = split(/[\:\-]/,$region);
my $duplicate_flag = 1024;
open(OUT,">".$output.".tmp.sam");
foreach my $bamfile (split(";",$bam))
{
	open(IN,"samtools view -h ".$bamfile." | ");
	while(my $read_line=<IN>)
	{
		if($read_line =~ /^@/)
		{
			print OUT $read_line;
		}
		my ($id,$flag,$chr_read,$pos,$MAPQ,$CIGAR,$chr_mate,$pos_mate,$frag_length,$seq,$qual) = split("\t",$read_line);
		if(($chr_read eq $chr)&&($chr_mate eq "=")&&(&min($pos,$pos_mate)<=$start)&&((&min($pos,$pos_mate) + abs($frag_length))>=$end)&&($flag < $duplicate_flag)&&(abs($frag_length) < 1000))
		{
			print OUT $read_line;
		}
	}
}
#system("samtools bam2fq -1 ".$output.".1.fq -2 ".$output.".2.fq ".$output.".tmp.sam");
system("samtools bam2fq ".$output.".tmp.sam > ".$output.".fq");
system("rm ".$output.".tmp.sam")
