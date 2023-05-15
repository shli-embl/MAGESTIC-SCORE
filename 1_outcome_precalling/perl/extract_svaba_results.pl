#!/usr/bin/env perl
#### Author: Shengdi Li ####
#### 2022 April. 11 ####
#### batch running SVABA Structural Variant calling ####

use Getopt::Long;
use Bio::SeqIO;
use strict;

my $window=5000;

GetOptions('h' => \my $display_help, 'in=s' => \my $input_meta_file, 'dir=s' => \my $vcf_dir, 'out=s' =>\my $output_file_summary, 'w=i' =>\$window);

if($display_help)
{
	print "Command:\n\tperl extract_svaba_result.pl [-h] -in INPUT_LIST -workdir PATH_TO_VCF_FILES -out OUTPUT_SUMMARY\n\n";
	print "Function: \n\tSummarize svaba results...\n\n";
	
	print "Usage:\n";
	print "\t-h\tPrint help info.\n";
	print "\t-in\tInput list file containing the meta data of samples.\n";
	print "\t-dir\tThe directory containing filtered and unfiltered vcfs.\n";
	print "\t-out\tOutput file containing the summary of svaba results.\n";
	print "\t-w\tWindow size within which an SV to be regarded as on target.\n";
	
	exit 1;
}

if((!$input_meta_file)||(!$vcf_dir)||(!$output_file_summary))
{
	print "Missing arguments.\n";
	print "Use -h for looking at usage...\n";
	exit 1;
}

open(INDEX,$input_meta_file);
my $header = <IN>;
open(OUT,">".$output_file_summary);
print OUT "sample_name\tOT_SV\tOT_INDEL\n";
while(my $read_line=<INDEX>)
{
	chomp $read_line;
	my ($sample_name,$path1,$path2,$chr,$strand,$pam_pos,$v_coord,$d_start,$d_end,$ref,$alt,$donor_wgs,$designed_v,$syn_v)=split("\t",$read_line);
	my @OT_SV;
	my @OT_INDEL;
	## read the unfiltered.sv.vcf file ##
	open(VCF,"$vcf_dir/unfil_vcf/$sample_name.svaba.unfiltered.sv.vcf");
	while(my $read_line1=<VCF>)
	{
		if($read_line1!~/^#/)
		{
			chomp $read_line1;
			my @read_column1=split("\t",$read_line1);
			if($read_column1[0] eq $chr)
			{
				if(($read_column1[1] > $v_coord - $window) && ($read_column1[1] < $v_coord + $window))
				{
					if($read_column1[7]=~/IMPRECISE/)
					{
						push @OT_SV,$read_column1[0].":".$read_column1[1]."_".$read_column1[3]."_to_".$read_column1[4].":".$read_column1[6].":IMPRECISE";
					}
					else
					{
						push @OT_SV,$read_column1[0].":".$read_column1[1]."_".$read_column1[3]."_to_".$read_column1[4].":".$read_column1[6];
					}
				}
			}
		}
	}
	## read the unfiltered.indel.vcf file ##
	open(VCF,"$vcf_dir/unfil_vcf/$sample_name.svaba.unfiltered.indel.vcf");
	while(my $read_line1=<VCF>)
	{
		if($read_line1!~/^#/)
		{
			chomp $read_line1;
			my @read_column1=split("\t",$read_line1);
			if($read_column1[0] eq $chr)
			{
				if(($read_column1[1] > $v_coord - $window) && ($read_column1[1] < $v_coord + $window))
				{
					if($read_column1[7]=~/IMPRECISE/)
					{
						push @OT_INDEL,$read_column1[0].":".$read_column1[1]."_".$read_column1[3]."_to_".$read_column1[4].":".$read_column1[6].":IMPRECISE";
					}
					else
					{
						push @OT_INDEL,$read_column1[0].":".$read_column1[1]."_".$read_column1[3]."_to_".$read_column1[4].":".$read_column1[6];
					}
				}
			}
		}
	}
	print OUT $sample_name,"\t",join(";",@OT_SV),"\t",join(";",@OT_INDEL),"\n";
}
