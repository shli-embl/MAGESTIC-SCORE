#!/usr/bin/env perl
#### Author: Shengdi Li ####
#### 2022 April. 19 ####
#### Calculate global sequence repetitiveness of a given sequence in a genome ####
#### This is the step 1 of the algorithm, which creates a bed file containing the found shortest unique substring (SUS) length at each genomic location ####
#### A masking window size can be defined, so that homologou sequence started within such window (center at test site) is not considered in search SUS. ####

use Getopt::Long;
use Bio::SeqIO;
use strict;

my $mask_offset=5000;
my $word_size=11;


GetOptions('h' => \my $display_help, 'out=s' =>\my $output_file, 'in=s'=>\my $genome_file, 'min_word_size:i'=>\$word_size, ,'mask_dist_offset:i'=>\$mask_offset);

if($display_help)
{
	print "Command:\n\tperl calculate_Igr_construct_sus_library.pl [-h] -in INPUT_GENOME_FASTA -out OUTPUT_FILE [-min_word_size MINIMUM_WORD_SIZE -mask_window_size MASKED_WINDOW_SIZE]\n\n";
	print "Function: \n\tCalculating the Global sequence repetativeness for a given list of target sites.\n\n";
	
	print "Usage:\n";
	print "\t-h\tPrint help info.\n";
	print "\t-in\tReference fasta sequence file.\n";
	print "\t-out\tOutput file.\n";
	print "\t-min_word_size {11}\tMinimum word size when searching for homology seed sequences. A too small word size will lead to significant computational time.\n";
	print "\t-mask_dist_offset {5000}\tLength of window centered at test location for masking. Any homology started within the masked window will not be considered as SUS.\n";
	
	exit 1;
}

if((!$output_file)||(!$genome_file)||(!$word_size)||(!$mask_offset))
{
	print "Missing arguments.\n";
	print "Use -h for looking at usage...\n";
	exit 1;
}

sub reverseComplement(){
	my $string = $_[0];
	my %nucleotide = ("A"=>"T","T"=>"A","C"=>"G","G"=>"C","N"=>"N");
	$string = reverse $string;
	for(my $i=0;$i<length($string);$i++)
	{
		substr($string,$i,1)=$nucleotide{substr($string,$i,1)}
	}
	return $string;
}

my %word_library;
## Loading genome sequences into hash table ##
my $genome_seqio=Bio::SeqIO->new(-file=>$genome_file, -format=>"fasta");
my %chromosome_sequence;
my @chrs;
while(my $chromosome_seq=$genome_seqio->next_seq)
{
	$chromosome_sequence{$chromosome_seq->display_id}=$chromosome_seq->seq;
	push @chrs,$chromosome_seq->display_id;
}

## step 1-1 ##
#### create N-mer library using genome sequence ####

foreach my $chromosome (@chrs)
{
		my $tmp_chr_seq = $chromosome_sequence{$chromosome};
		my $length = length($tmp_chr_seq);
		my $i=0;
		while($i < ($length-$word_size+1))
		{
			$word_library{substr($tmp_chr_seq,$i,$word_size)}{$chromosome}{$i}=1;
			$i++;			
		}		
}
print "finished creating N-mer library...\n";
open(OUT,">".$output_file);
## step 1-2 ##
#### for each position, calculate length of SUS ####
#### this is done based on the indexed N-mer poistions ####
foreach my $chromosome (@chrs)
{
	print "processing ".$chromosome." ...\n";
	my $tmp_chr_seq = $chromosome_sequence{$chromosome};
	my $i=0;
	while($i < (length($tmp_chr_seq)-$word_size+1))
	{
		my $current_word_1 = substr($tmp_chr_seq,$i,$word_size);
		my $current_word_2 = &reverseComplement(substr($tmp_chr_seq,$i,$word_size));
		my $SUS =0;
		my $SUS_trim=0;	
		my $best_hit=0;
		## forward strand ##
		foreach my $chr(keys %{$word_library{$current_word_1}})
		{
			foreach my $pos(keys %{$word_library{$current_word_1}{$chr}})
			{
				## if the string is within defined offset distance to the tested location, then leave out the string ##
				if(($chr eq $chromosome)&&(abs($i - $pos) < $mask_offset))
				{
				}
				else
				{
					my $tmp_str1=$current_word_1;
					my $tmp_str2=$current_word_1;
					my $tmp_chr_seq2=$chromosome_sequence{$chr};
					my $last_base1=$i+$word_size;
					my $last_base2=$pos+$word_size;
					while(($tmp_str1 eq $tmp_str2)&&($last_base1<length($tmp_chr_seq))&&($last_base2<length($tmp_chr_seq2)))
					{
						$tmp_str1=$tmp_str1.substr($tmp_chr_seq,$last_base1,1);
						$tmp_str2=$tmp_str2.substr($tmp_chr_seq2,$last_base2,1);
						$last_base1++;
						$last_base2++;
					}
					if(length($tmp_str1)>$SUS)
					{
						$SUS = length($tmp_str1);
						$SUS_trim = $SUS - $word_size;
					}
				}
			}				
		}
		## backward strand ##
		foreach my $chr(keys %{$word_library{$current_word_2}})
		{
			foreach my $pos(keys %{$word_library{$current_word_2}{$chr}})
			{
				## if the string is within defined offset distance to the tested location, then leave out the string ##
				if(($chr eq $chromosome)&&(abs($i - $pos) < $mask_offset))
				{
				}
				else
				{
					my $tmp_str1=$current_word_2;
					my $tmp_str2=$current_word_2;
					my $tmp_chr_seq2=$chromosome_sequence{$chr};
					my $last_base1=$i+$word_size;
					my $last_base2=$pos-1;
					while(($tmp_str1 eq $tmp_str2)&&($last_base1<length($tmp_chr_seq))&&($last_base2>=0))
					{
						$tmp_str1=&reverseComplement(substr($tmp_chr_seq,$last_base1,1)).$tmp_str1;
						$tmp_str2=substr($tmp_chr_seq2,$last_base2,1).$tmp_str2;
						$last_base1++;
						$last_base2--;
					}
					if(length($tmp_str1)>$SUS)
					{
						$SUS = length($tmp_str1);
						$SUS_trim = $SUS - $word_size;
					}
				}
			}
		}
		
		## output SUS ##
		print OUT $chromosome,"\t",$i+1,"\t",$i+1,"\t",$SUS,"\t",$SUS_trim,"\n";
		$i++;
	}
}