#!/usr/bin/env perl
#### Author: Shengdi Li ####
#### 2022 Nov. 2 ####
#### Summarizing TRL scores from the TRL program output ####
use Getopt::Long;
use strict;

GetOptions('h' => \my $display_help, 'in=s' => \my $input_meta_file, 'workdir=s' => \my $data_dir, 'w=s' => \my $w, 'out=s' =>\my $output_file_summary);
if($display_help)
{
	print "Command:\n\tperl summary_trl.pl [-h] -in INPUT_TARGET_LIST -workdir PATH_TO_TRL_OUTPUT -out OUTPUT_SUMMARY -w WINDOW_SIZE\n\n";
	print "Function: \n\tSummarizing output into a single file from TRL...\n\n";
	
	print "Usage:\n";
	print "\t-h\tPrint help info.\n";
	print "\t-in\tInput list of targets (col1: id, col2: chr, col3: strand, col4: PAM coord, col5: guide sequence).\n";
	print "\t-workdir\tDirectory containing the outputs of TRL.\n";
	print "\t-w\tWindow size (only used for searching the naming of ir outputs; e.g. 250bp window size has the suffix \"xxx.trl2500x2.txt\").\n";
	print "\t-out\tOutput file.\n";

	exit 1;
}

if((!$input_meta_file)||(!$data_dir)||(!$output_file_summary)||(!$w))
{
	print "Missing arguments.\n";
	print "Use -h for looking at usage...\n";
	exit 1;
}


#start processing samples
open(INDEX,$input_meta_file);
my $header =<INDEX>;
open(OUT,">".$output_file_summary);
print OUT "#id\tTRL",$w,"x2\n";
while(my $read_line=<INDEX>)
{
	chomp $read_line;
	my ($id,$chr,$strand,$coord,$guide)=split("\t",$read_line);
	open(TMP,$data_dir."/".$id.".trl".$w."x2.out");
	my $line=<TMP>;
	$line=<TMP>;
	chomp $line;
	my @read_column=split("\t",$line);
	print OUT $id,"\t",$read_column[1],"\n";
}