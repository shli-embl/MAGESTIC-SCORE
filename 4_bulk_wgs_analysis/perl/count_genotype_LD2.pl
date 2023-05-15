#!/usr/bin/env perl
#### Author: Shengdi Li ####
#### 2023 Feb. 16 ####
#### count reads after remapping ####

use Getopt::Long;
use Bio::SeqIO;
use strict;

GetOptions('h' => \my $display_help, 'dir=s' => \my $dir,'bam=s' => \my $bam, 'out=s' =>\my $output);
if($display_help)
{
	print "Command: \n\tperl count_genotype_LD2.pl [-h] -bam INPUT_BAM -out OUTPUT\n\n";
	print "Function: \n\tCount reads after remapping...\n\n";
	print "Usage: \n";
	print "\t-h\tPrint help info.\n";
	print "\t-bam\tThe input bam files (separated with ; if there are multiple).\n";
	print "\t-out\tOutput file.\n";
	exit 1;
}

if((!$dir)||(!$bam)||(!$output))
{
	print "Missing arguments.\n";
	print "Use -h for looking at usage...\n";
	exit 1;
}

open(OUT,">".$output);
foreach my $bamfile (split(";",$bam))
{
#	print $bamfile,"\n";
	open(IN,"samtools view -h ".$dir."/".$bamfile." | ");
	my $mut = 0;
	my $wt = 0;
	my $deletion = 0;
	while(my $read_line=<IN>)
	{
		my ($id,$flag,$chr_read,$pos,$MAPQ,$CIGAR,$chr_mate,$pos_mate,$frag_length,$seq,$qual) = split("\t",$read_line);
		if($MAPQ > 0)
		{
			if($chr_read eq "wt")
			{
				$wt++;
			}
			elsif($chr_read eq "mut")
			{
				$mut++;
			}
			elsif(($chr_read eq "deletion1")||($chr_read eq "deletion2"))
			{
				$deletion++;
			}
		}
	}
	print OUT $bamfile,"\t",$wt,"\t",$mut,"\t",$deletion,"\n";
}

