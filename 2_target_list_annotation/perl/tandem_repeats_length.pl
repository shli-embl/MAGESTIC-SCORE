#!/usr/bin/env perl
#### Author: Shengdi Li ####
#### 2022 April. 19 ####
#### Calculate length metrics of tandem repeats present at target window ####
## Author: Shengdi Li
use strict;
use Bio::SeqIO;
use Getopt::Long;

my $min_match_bp=3;
my $max_mismatch=1;
my $threshold_match_filter=3;

sub min()
{
	if($_[0]>$_[1])
	{
		return $_[1];
	}
	else
	{
		return $_[0];
	}
}
sub max()
{
	if($_[0]<$_[1])
	{
		return $_[1];
	}
	else
	{
		return $_[0];
	}
}

GetOptions('h' => \my $display_help, 'chr=s' => \my $chr, 'coord=s' => \my $coord, 'out=s' =>\my $out, 'genome=s' =>\my $genome, 'window=s' =>\my $window, 'minM:i'=>\$min_match_bp, 'maxMM:i'=>\$max_mismatch, 'filM:i'=>\$threshold_match_filter);

if($display_help)
{
	print "Command:\n\tperl tandem_repeats.pl [-h] -in INPUT_META_LIST -out OUTFILE -w WINDOW_SIZE [-minM MIN_MATCH_BP -maxMM MAX_MISMATCH_BP -filM MIN_MATCH_BP_OUTPUT]";
	print "Function: \n\tSearch for tandem repeats present at the (window_size X 2) target region. A minimum of 7bp perfect match is considered as seed region and 1bp is allowed by default to extend the matched homology...\n\n";
	
	print "Usage:\n";
	print "\t-h\tPrint help info.\n";
	print "\t-chr\tFrom which chromosome id the sequence should be extracted.\n";
	print "\t-coord\tFrom which coordinate on the chromosome id the sequence should be treated as target.\n";
	print "\t-window\tLength of the sequence up- and down-stream the coordinate to search for TRLs.\n";
	print "\t-out\tOutput file.\n";
	print "\t-genome\tGenome sequence file.\n";
	print "\t-minM {3}\tMinimal number of perfect matched bases for a pair of sequences to be regarded as potential homology pair (seed).\n";
	print "\t-maxMM {1}\tMax number of mismatches (when extending the seeds) when concatenating adjacent perfect matched pairs\n";
	print "\t-filM {3}\tMinimal number of matched bases (mismatch allowed)of homology pairs to be output.\n";
	
	exit 1;
}
if((!$chr)||(!$coord)||(!$out)||(!$genome)||(!$window))
{
	print "Missing arguments.\n";
	print "Use -h for looking at usage...\n";
	exit 1;
}

#Load yeast_genome;
my $chr_seq;
my $seqio=Bio::SeqIO->new(-file=>$genome,-format=>"fasta");
while(my $chromosome_seq=$seqio->next_seq)
{
	if($chromosome_seq->display_id eq $chr)
	{
		$chr_seq=$chromosome_seq->seq;
		last;
	}
}

my %vocabulary;
my $left_seq=substr($chr_seq,&max(0,$coord-1-$window),$coord-1-&max(0,$coord-1-$window));
my $right_seq=substr($chr_seq,$coord-1,&min(length($chr_seq)+1-$coord,$window));

#Build n-mer vocabulary
for(my $i=1;$i<(length($left_seq)-$min_match_bp+1);$i++)
{
	if($vocabulary{substr($left_seq,$i-1,$min_match_bp)}{"L"})
	{	
		$vocabulary{substr($left_seq,$i-1,$min_match_bp)}{"L"}=$vocabulary{substr($left_seq,$i-1,$min_match_bp)}{"L"}.",".$i;
	}
	else
	{
		$vocabulary{substr($left_seq,$i-1,$min_match_bp)}{"L"}=$i;
	}
}
for(my $i=1;$i<(length($right_seq)-$min_match_bp+1);$i++)
{
	if($vocabulary{substr($right_seq,$i-1,$min_match_bp)}{"R"})
	{	
		$vocabulary{substr($right_seq,$i-1,$min_match_bp)}{"R"}=$vocabulary{substr($right_seq,$i-1,$min_match_bp)}{"R"}.",".$i;
	}
	else
	{
		$vocabulary{substr($right_seq,$i-1,$min_match_bp)}{"R"}=$i;
	}
}

#Merge matched pairs, get a list of perfectly matched intervals
my %position_matched;
foreach my $word(keys %vocabulary)
{
	if(($vocabulary{$word}{"L"})&&($vocabulary{$word}{"R"}))
	{
		my @word_left=split(",",$vocabulary{$word}{"L"});
		my @word_right=split(",",$vocabulary{$word}{"R"});
		foreach my $coord_left(@word_left)
		{
			foreach my $coord_right(@word_right)
			{
				$position_matched{$coord_left}{$coord_right}=$word
			}
		}
	}
}
##scan from first position and merge adjacent matches
for(my $i=1; $i<length($left_seq);$i++)
{
	if($position_matched{$i})
	{
		foreach my $right(%{$position_matched{$i}})
		{
			my $extension=1;
			while($position_matched{$i+$extension}{$right+$extension})
			{
				$position_matched{$i}{$right}=substr($left_seq,$i-1,length($position_matched{$i}{$right})+1);
				delete $position_matched{$i+$extension}{$right+$extension};
				$extension++;
			}
		}
	}
}

#concatenate flanking homolog elements if <= min_misatch_allowed
for(my $i=1; $i<length($left_seq);$i++)
{
	if($position_matched{$i})
	{
		foreach my $right(keys%{$position_matched{$i}})
		{
			my $length=length($position_matched{$i}{$right});
			my $left_start=$i;
			my $right_start=$right;
			my $left_end=$i+$length-1;
			my $right_end=$right+$length-1;
			for(my $j=$left_end+1;$j<=($left_end+$max_mismatch+1);$j++)
			{
				if($position_matched{$j})
				{
					foreach my $right_j(keys %{$position_matched{$j}})
					{
						if((($right_j-$right_end)<=($max_mismatch+1))&&(($right_j-$right_end)>0))
						{
							$position_matched{$i}{$right}=$position_matched{$i}{$right}."(";
							if($left_end==($j-1))
							{
								$position_matched{$i}{$right}=$position_matched{$i}{$right}."-/";
							}
							else
							{
								$position_matched{$i}{$right}=$position_matched{$i}{$right}.lc(substr($left_seq,$left_end,$j-$left_end-1))."/";
							}
							if($right_end==($right_j)-1)
							{
								$position_matched{$i}{$right}=$position_matched{$i}{$right}."-";
							}
							else
							{
								$position_matched{$i}{$right}=$position_matched{$i}{$right}.lc(substr($right_seq,$right_end,$right_j-$right_end-1))
							}
							$position_matched{$i}{$right}=$position_matched{$i}{$right}.")".$position_matched{$j}{$right_j};
							delete $position_matched{$j}{$right_j};
						}
						
					}
				}
			}
		}
	}
}
#get homology length for ranking
my %n_matched_bases;
foreach my $left_key(keys %position_matched)
{
	foreach my $right_key(keys %{$position_matched{$left_key}})
	{
		my $tmp=$position_matched{$left_key}{$right_key};
		$tmp=~s/[a-z]//g;
		$tmp=~s/[\(\)\/]//g;
		$n_matched_bases{$left_key."_".$right_key}=length($tmp);
	}
}

#calculate the longest shared substring, allowing 1 mismatch, between N bp left to pam site vs. N bp right to pam site
my $metricD;
my $homolog_seq;
foreach my $key(sort {$n_matched_bases{$b}<=>$n_matched_bases{$a}} keys %n_matched_bases)
{
	my @split=split("_",$key);
	my $right_key=$split[1];
	my $left_key=$split[0];
	my $offset=$right_key+length($left_seq)-$left_key;
	if($n_matched_bases{$key}>$threshold_match_filter)
	{
		if($metricD < $n_matched_bases{$key})
		{
			$metricD=$n_matched_bases{$key};
			$homolog_seq=$position_matched{$left_key}{$right_key};
		}
	}
}

## output 
open(OUT,">".$out);
print OUT "#id\tTRL".$window."x2\tbest_match_sequence\n";
if($metricD)
{
	print OUT "X\t",$metricD;
}
else
{
	print OUT "X\t0";
}
print OUT "\t",$homolog_seq,"\n";


