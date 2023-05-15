#!/usr/bin/env perl
#### Author: Shengdi Li ####
#### 2023 Feb. 16 ####
#### count read number containing the edited or unedited genotype at a genomic range ####

use Getopt::Long;
use Bio::SeqIO;
use strict;

GetOptions('h' => \my $display_help, 'bam=s' => \my $bam, 'variant=s' => \my $variant, 'up=s' => \my $up_seq, 'down=s' => \my $down_seq, 'out=s' =>\my $output);
if($display_help)
{
	print "Command: \n\tperl count_genotype.pl [-h] -bam INPUT_BAM -chr CHR -start START -end END -gt1 GENOTYPE1 -gt2 GENOTYPE2 -out OUTPUT\n\n";
	print "Function: \n\tCount the number of reads in a given genomic range that contain a substring matching the given genotypes...\n\n";
	print "Usage: \n";
	print "\t-h\tPrint help info.\n";
	print "\t-bam\tThe input bam files (separated with ; if there are multiple).\n";
	print "\t-variant\tThe variant genotype to match (e.g chr8_287671_AGGTGG_GAATTC).\n";
	print "\t-up\tThe upstream sequence to the scaned region.\n";
	print "\t-down\tThe downstream sequence to the scaned region.\n";
	print "\t-out\tOutput file.\n";
	exit 1;
}

my %complement = ("A"=>"T","G"=>"C","C"=>"G","T"=>"A","N"=>"N");


if((!$bam)||(!$variant)||(!$output)||(!$up_seq)||(!$down_seq))
{
	print "Missing arguments.\n";
	print "Use -h for looking at usage...\n";
	exit 1;
}

sub reverseComplement(){
	my $seq = $_[0];
	$seq = reverse $seq;
	for(my $i=0;$i<length($seq);$i++)
	{
		if($complement{substr($seq,$i,1)})
		{
			substr($seq,$i,1)=$complement{substr($seq,$i,1)}
		}
		else
		{
			print "found unknown character in the alignment ",substr($seq,$i,1),"\n";
			exit 0
		}
	}
	return $seq;
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
my ($chr,$start,$genotype1,$genotype2) = split("_",$variant);
my $end = $start + length($genotype1) - 1;
my $duplicate_flag = 1024;
my $genotype1_up = $up_seq.$genotype1;
my $genotype1_up_rc = &reverseComplement($genotype1_up);
my $genotype1_down = $genotype1.$down_seq;
my $genotype1_down_rc = &reverseComplement($genotype1_down);
my $genotype2_up = $up_seq.$genotype2;
my $genotype2_up_rc = &reverseComplement($genotype2_up);
my $genotype2_down = $genotype2.$down_seq;
my $genotype2_down_rc = &reverseComplement($genotype2_down);
my $count_gt1 = 0;
my $count_gt2 = 0;
foreach my $bamfile (split(";",$bam))
{
	open(IN,"samtools view ".$bamfile." | ");
	while(my $read_line=<IN>)
	{
		my ($id,$flag,$chr_read,$pos,$MAPQ,$CIGAR,$chr_mate,$pos_mate,$frag_length,$seq,$qual) = split("\t",$read_line);
	#	print $chr_read,"\t",$chr_mate,"\t",&min($pos,$pos_mate),"\t",(&min($pos,$pos_mate) + abs($frag_length)),"\t",$flag,"\n";
		if(($chr_read eq $chr)&&($chr_mate eq "=")&&(&min($pos,$pos_mate)<=$start)&&((&min($pos,$pos_mate) + abs($frag_length))>=$end)&&($flag < $duplicate_flag)&&(abs($frag_length) < 1000))
		{
	#		print $read_line;
	#		print $seq,"\n";
			if(($seq =~ /$genotype1_up/)||($seq =~ /$genotype1_up_rc/)||($seq =~ /$genotype1_down/)||($seq =~ /$genotype1_down_rc/))
			{
				$count_gt1++;
			}
			elsif(($seq =~ /$genotype2_up/)||($seq =~ /$genotype2_up_rc/)||($seq =~ /$genotype2_down/)||($seq =~ /$genotype2_down_rc/))
			{
				$count_gt2++;
			}
		}
	}
}	

open(OUT,">".$output);
print OUT "#id\t","variant","\t",$genotype1,"\t",$genotype2,"\n";
print OUT $bam,"\t",$variant,"\t",$count_gt1,"\t",$count_gt2,"\n";
