#!/usr/bin/env perl
#### Author: Shengdi Li ####
#### 2022 April. 19 ####
#### Extract strand-specific RNA-seq read density at a target window of 50bp ####
## Author: Shengdi Li
use strict;
use Getopt::Long;
my $window_size = 50;
#### Load external data ####
GetOptions('h' => \my $display_help, 'in=s' => \my $input_meta_file, 'out=s' => \my $outfile, 'w:i'=>\ $window_size, 'forward=s' => \my $forward_depth, 'reverse=s' => \my $reverse_depth);

if($display_help)
{
	print "Command:\n\tperl extract_transcription_level.pl [-h] -in INPUT_META_LIST -out OUTFILE -forward FORWARD_STRAND_READ_DEPTH -reverse REVERSE_STRAND_READ_DEPTH [-w WINDOW_SIZE]\n";
	print "Function: \n\tExtract strand-specific RNA-seq read density at a list of target sites...\n\n";
	
	print "Usage:\n";
	print "\t-h\tPrint help info.\n";
	print "\t-in\tInput sample meta sheet.\n";
	print "\t-out\tOutput file path. \n";
	print "\t-forward\tBed files of read depth for forward reads.\n";
	print "\t-reverse\tBed files of read depth for reverse reads.\n";
	print "\t-w {50}\tLength of window to extract mean read depth.\n";
	
	exit 1;
}
if((!$input_meta_file)||(!$outfile)||(!$window_size)||(!$forward_depth)||(!$reverse_depth))
{
	print "Missing arguments.\n";
	print "Use -h for looking at usage...\n";
	exit 1;
}

my %forward;
my %reverse;
open(FORWARD,$forward_depth);
open(REVERSE,$reverse_depth);
while(my $read_line=<FORWARD>)
{
	chomp $read_line;
	my ($chr,$pos,$value) = split("\t",$read_line);
	$forward{$chr}{$pos}=$value;
}
while(my $read_line=<REVERSE>)
{
	chomp $read_line;
	my ($chr,$pos,$value) = split("\t",$read_line);
	$reverse{$chr}{$pos}=$value;
}
print "Finish Loading RNA-seq dataset...\n";
open(IN,$input_meta_file);
my $header = <IN>;
open(OUT,">".$outfile);
print OUT "#id\tTranscription_Cas9_bind_template_strand\tTranscription_Cas9_bind_non_template_strand\n";
while(my $read_line=<IN>)
{
	chomp $read_line;
	my ($id,$chr,$strand,$coord,$guide)=split("\t",$read_line);

	my $RD_forward=0;
	my $RD_reverse=0;
	for(my $i=($coord-$window_size/2+1);$i<=($coord+$window_size/2);$i++)
	{
		$RD_forward+=$forward{$chr}{$i};
		$RD_reverse+=$reverse{$chr}{$i};
	}
	
	$RD_forward=sprintf("%.3f",log($RD_forward+1)/log(2)-log($window_size)/log(2));
	$RD_reverse=sprintf("%.3f",log($RD_reverse+1)/log(2)-log($window_size)/log(2));
	
	print OUT $id,"\t";
	### First output expression level when Cas9 binds template strand; next when Cas9 binds non-template strand ##
	if($strand eq "+")
	{
		print OUT $RD_forward,"\t",$RD_reverse,"\n";
	}
	else
	{
		print OUT $RD_reverse,"\t",$RD_forward,"\n";
	}
	
}