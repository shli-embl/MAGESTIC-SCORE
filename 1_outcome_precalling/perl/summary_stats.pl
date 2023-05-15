#!/usr/bin/env perl
#### Author: Shengdi Li ####
#### 2022 Mar. 29 ####
#### Summarizing statistics from mapped reads (.bam) ####
#### Get read numbers, mapping rate and coverage depth per WGS dataset ####
use Getopt::Long;
use Bio::SeqIO;
use strict;

GetOptions('h' => \my $display_help, 'in=s' => \my $input_meta_file, 'workdir=s' => \my $data_dir, 'out=s' =>\my $output_file_summary);
if($display_help)
{
	print "Command:\n\tperl summary_stats.pl [-h] -in INPUT_LIST -workdir PATH_TO_DATA_FOLDER -out OUTPUT_SUMMARY \n\n";
	print "Function: \n\tSummarizing statistics for a given list of samples (with pre-processing done with read_process_batch_final.pl)...\n\n";
	
	print "Usage:\n";
	print "\t-h\tPrint help info.\n";
	print "\t-in\tInput list file containing the meta data of samples.\n";
	print "\t-workdir\tProvide the output folder of read_process_batch_final.pl.\n";
	print "\t-out\tOutput file containing the summary of per-sample statistics.\n";

	exit 1;
}

if((!$input_meta_file)||(!$data_dir)||(!$output_file_summary))
{
	print "Missing arguments.\n";
	print "Use -h for looking at usage...\n";
	exit 1;
}


#start processing samples
open(INDEX,$input_meta_file);
my $header=<INDEX>;
open(OUT,">".$output_file_summary);
print OUT "#sample_name\ttotal_reads\tmapping_rate\tave_coverage\n";
while(my $read_line=<INDEX>)
{
	chomp $read_line;
	my ($sample_name,$path1,$path2,$chr,$strand,$pam_coord,$v_coord,$d_start,$d_end,$ref,$alt,$donor_WGS,$target_vars,$syn_vars)=split("\t",$read_line);
	if(-e $data_dir."/bam/".$sample_name.".clean.bam")
	{
		## calculate per-base coverage from samtools depth output ##
		system("samtools depth -a ".$data_dir."/bam/".$sample_name.".clean.bam > tmp.".$sample_name.".depth");
		open(INPUT,"tmp.".$sample_name.".depth");
		my $sum_mapped_bases = 0;
		my $n_sites = 0;
		while(my $read_line1=<INPUT>)
		{
			chomp $read_line1;
			my @read_column1=split("\t",$read_line1);
			if($read_column1[0] ne "chrM")
			{
				$sum_mapped_bases+=$read_column1[2];
				$n_sites++;
			}
		
		}
		my $average_depth_per_base;
		if($n_sites!=0)
		{
			$average_depth_per_base = sprintf("%.2f",$sum_mapped_bases/$n_sites);
		}
		else
		{
			$average_depth_per_base = "non_zero";
		}
		system("rm tmp.".$sample_name.".depth");
		## read mapping statistics from .flagstats ##
		my ($total,$map_rate);
		open(FLAGSTAT,$data_dir."/stats/".$sample_name.".flagstat");
		while(my $read_line=<FLAGSTAT>)
		{
			if($read_line=~/^(\d+) \+ \d+ in total/)
			{
				$total = $1;
			}
			if($read_line=~/mapped \((\d+\.\d+\%) \:/)
			{
				$map_rate = $1;
			}
		}
		print OUT $sample_name,"\t",$total,"\t",$map_rate,"\t",$average_depth_per_base,"\n";
	}
	else
	{
		print OUT $sample_name,"\tNA\tNA\tNA\n";
	}
}