#!/usr/bin/env perl
#### Author: Shengdi Li ####
#### 2022 Nov. 2 ####
#### Summarizing off target prediction from Casoffinder ####
use Getopt::Long;
use strict;

GetOptions('h' => \my $display_help, 'in=s' => \my $input_meta_file, 'workdir=s' => \my $data_dir, 'out=s' =>\my $output_file_summary, 'm=s' =>\my $m, 'out2=s' =>\my $output_list);
if($display_help)
{
	print "Command:\n\tperl summary_offtarget.pl [-h] -in INPUT_TARGET_LIST -workdir PATH_TO_CASOFFINDER_OUTPUT -m MISMATCHES -out OUTPUT_SUMMARY [-out2 OUTPUT_OFFTARGET_LIST] \n\n";
	print "Function: \n\tSummarizing output into a single file from Cas-offinder...\n\n";
	
	print "Usage:\n";
	print "\t-h\tPrint help info.\n";
	print "\t-in\tInput list of targets (col1: id, col2: chr, col3: strand, col4: PAM coord, col5: guide sequence).\n";
	print "\t-workdir\tDirectory containing the outputs of Casoffinder.\n";
	print "\t-m\tMaximum number of mismatches to look for.\n";
	print "\t-out\tOutput file.\n";

	exit 1;
}

if((!$input_meta_file)||(!$data_dir)||(!$output_file_summary)||(!$m))
{
	print "Missing arguments.\n";
	print "Use -h for looking at usage...\n";
	exit 1;
}


#start processing samples
open(INDEX,$input_meta_file);
my $header =<INDEX>;
open(OUT,">".$output_file_summary);
print OUT "#id";
if($output_list){
	open(OUT2,">".$output_list);
	print OUT2 "#id\tguide\tchr\tcoord\tmatch\tstrand\tmismatch\n";
}
for(my $i = 0;$i<=$m;$i++)
{
	print OUT "\tm",$i
}
print OUT "\n";
while(my $read_line=<INDEX>)
{
	chomp $read_line;
	my ($id,$chr,$strand,$coord,$guide)=split("\t",$read_line);
	open(TMP,$data_dir."/".$id.".ot.txt");
	my %ms;
	for(my $i = 0;$i<=$m;$i++)
	{
		$ms{$i}=0
	}
	while(my $read_line1 = <TMP>)
	{
		chomp $read_line1;
		my ($guide,$chr,$coord,$match,$strand,$mm) = split("\t",$read_line1);
		$ms{$mm}++;
		if($output_list){
			print OUT2 $id,"\t",$guide,"\t",$chr,"\t",$coord,"\t",$match,"\t",$strand,"\t",$mm,"\n";
		}
	}
	print OUT $id;
	for(my $i = 0;$i<=$m;$i++)
	{
		print OUT "\t",$ms{$i}
	}
	print OUT "\n";
}