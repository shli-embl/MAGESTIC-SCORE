#!/usr/bin/env perl
#### Author: Shengdi Li ####
#### 2022 April. 19 ####
#### Calculate global sequence repetitiveness of a given sequence in a genome ####
#### This is the step 2 of the algorithm, which creates a bed file containing the found shortest unique substring (SUS) length at each genomic location ####
#### A masking window size can be defined, so that homologou sequence started within such window (center at test site) is not considered in searching SUS. ####

use Getopt::Long;
use Bio::SeqIO;
use strict;

GetOptions('h' => \my $display_help, 'chr=s' =>\my $chr, 'coord=s' =>\my $coord, 'out=s' =>\my $out, 'w=i'=>\my $window_size,'bed_obs=s'=>\my $bed_observe, 'bed_exp=s'=>\my $bed_expect);

if($display_help)
{
	print "Command:\n\tperl calculate_igr_get_score.pl [-h] -chr CHROMOSOME -coord COORDINATE -out OUTPUT_FILE -w WINDOW_SIZE -bed_obs BED_FILE_SUS -bed_exp BED_FILE_SUS_SHUFFLED \n\n";
	print "Function: \n\tCalculating the Global sequence repetativeness for a given list of target sites.\n\n";
	
	print "Usage:\n";
	print "\t-h\tPrint help info.\n";
	print "\t-chr\tSite chromosome id.\n";
	print "\t-coord\tSite coordinate.\n";
	print "\t-out\tOutput file.\n";
	print "\t-w\tWindow_length for calculating AUC.\n";
	print "\t-bed_obs\tBed file which contains five columns: 1)chr 2)position 3)position 4)SUS 5)trimmed SUS.\n";
	print "\t-bed_exp\tBed file in same format as -bedobs input using shuffled genome.\n";
	
	exit 1;
}

if((!$chr)||(!$coord)||(!$out)||(!$window_size)||(!$bed_observe)||(!$bed_expect))
{
	print "Missing arguments.\n";
	print "Use -h for looking at usage...\n";
	exit 1;
}

#### get mean SUS and word size from bed expect ####
open(IN,$bed_expect);
my $word_size=0;
my $sum=0;
my $i=0;
while(my $read_line=<IN>)
{
	chomp $read_line;
	my ($chr,$start,$end,$SUS,$SUS_trim)=split("\t",$read_line);
	if(!$word_size)
	{
		$word_size=$SUS-$SUS_trim;
	}
	$sum = $sum + $word_size + $SUS_trim;
	$i++;
}
my $mean_SUS = $sum/$i;
close IN;
#### load sample sheet and calculate AUC ####
open(BED,$bed_observe);
my %SUS_obs;
while(my $read_line=<BED>)
{
	chomp $read_line;
	my ($chr,$start,$end,$SUS,$SUS_trim)=split("\t",$read_line);
	$SUS_obs{$chr}{$start}=$SUS_trim + $word_size;
}
## calculate log(AUC_OBS/AUC_EXP)
open(OUT,">".$out);
print OUT "#id\tIGR",$window_size,"\n";
my $AUC=0;
my $n=0;
for(my $i=$coord-$window_size/2+1;$i<=($coord+$window_size/2);$i++)
{
	if($SUS_obs{$chr}{$i}>0)
	{
		$AUC+=$SUS_obs{$chr}{$i};
		$n++;
	}
}
if($n)
{
	print OUT "X","\t",sprintf("%.3f",log($AUC)/log(10)-log($mean_SUS*$n)/log(10)),"\n";
}
else
{
	print OUT "X","\t","\n";
}

