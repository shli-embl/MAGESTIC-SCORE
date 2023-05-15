#!/usr/bin/env perl
#### Author: Shengdi Li ####
#### 2022 Nov. 15 ####
#### Count on-target reads ####

use Getopt::Long;
use Bio::SeqIO;
use strict;

#start reading parameters
my $input_meta_file;
my $data_dir;
my $output_file_summary;

GetOptions('h' => \my $display_help, 'in=s' => \my $input_meta_file, 'workdir=s' => \my $data_dir, 'out=s' =>\my $output_file_summary);

if($display_help)
{
	print "Command:\n\tperl summary_count_ontarget_reads.pl [-h] -in INPUT_LIST -workdir PATH_TO_DATA_FOLDER -out OUTPUT_SUMMARY \n\n";
	print "Function: \n\tSummarizing on-target read coverage for a given list of sample (bam files)...\n\n";
	
	print "Usage:\n";
	print "\t-h\tPrint help info.\n";
	print "\t-in\tInput list file containing the meta data of samples.\n";
	print "\t-workdir\tThe path where WGS result sub-folders can be found (bams are required).\n";
	print "\t-out\tOutput file containing the summary of per-sample on-target coverage information.\n";

	exit 1;
}

if((!$input_meta_file)||(!$data_dir)||(!$output_file_summary))
{
	print "Missing arguments.\n";
	print "Use -h for looking at usage...\n";
	exit 1;
}
sub dec2bin {
    my $str = unpack("B32", pack("N", shift));
    $str =~ s/^0+(?=\d)//;   # otherwise you'll get leading zeros
    return $str;
}
sub if_duplicate
{
	my $flag = &dec2bin($_[0]);
	## position of duplicate flag is on 11th
	if(length($flag) < 11)
	{
		return 0;
	}
	else
	{
		return substr($flag,10,1);
		
	}
}
sub length_CIGAR
{
	my $CIGARstring=$_[0];
	my $read_line=$_[1];
	my $length_CIGAR=0;
	while($CIGARstring)
	{
		if($CIGARstring=~/^(\d+)([IDMSH])/)
		{
			my $tmp_length=$1;
			my $tmp_type=$2;
			if(($tmp_type eq "M")||($tmp_type eq "D"))
			{
				$length_CIGAR+=$tmp_length;
			}
			$CIGARstring=~s/^\d+[IDMSH]//;
		}
		elsif($CIGARstring eq "*")
		{
			$CIGARstring=~s/\*//;
		}
		else
		{
			print "CIGAR string unrecognized:",$CIGARstring,"\n";
			print $read_line,"\n";
			exit;
		}
	}
	return $length_CIGAR;
}	

#start processing samples
open(INDEX,$input_meta_file);
my $header = <INDEX>;
open(OUT,">".$output_file_summary);
print OUT "#id\tcoverage_depth\tfrags_overlap_genome\tfrags_overlap_cassete\tfrags_translocation\n";
while(my $read_line=<INDEX>)
{
	chomp $read_line;
	my ($sample_name,$path1,$path2,$chr,$strand,$pam_pos,$v_coord,$d_start,$d_end,$ref,$alt,$donor_wgs,$designed_v,$syn_v)=split("\t",$read_line);	
	
	if(-e $data_dir."/bam/".$sample_name.".clean.bam.target_filtering.anno")
	{
		open(ANNO,$data_dir."/bam/".$sample_name.".clean.bam.target_filtering.anno");
		my $flagstart=0;
		my $variant_cov = 0;
		my $target_read = 0;
		my $excluded_read = 0;
		my $translocation_read = 0;
		my %processed_id;		
		while(my $read_line=<ANNO>)
		{
			chomp $read_line;
			if($flagstart)
			{	
				my ($action,$id,$flag)=split(/[ \t]/,$read_line);
				if((!$processed_id{$id})&&(!&if_duplicate($flag)))
				{
					if($action=~/INCLUDE/)
					{
						$target_read++;
						open(BAM,"samtools view -h ".$data_dir."/bam/".$sample_name.".sort.dedup.recal2.bam | grep \"".$id."\" | ");
						while(my $read_line1=<BAM>)
						{
							chomp $read_line1;
							my ($read_id,$FLAG,$chr2,$read_start,$r4,$r5)=split("\t",$read_line1);
							my $read_end=$read_start + &length_CIGAR($r5)-1;
							if(($chr eq $chr2)&&($v_coord>=$read_start)&&(($v_coord+length($ref)-1)<=$read_end))
							{
								$variant_cov++;
							}
						}
					}
					elsif($action=~/DISCARD/)
					{
						$excluded_read++;
					}
					elsif($action=~/REDITRANS/)
					{
						$translocation_read++;
					}
					$processed_id{$id}++;
				}
				
			}
			elsif($read_line=~/READ\_INFO/)
			{
				$flagstart=1;
			}
		}
		print OUT $sample_name,"\t",$variant_cov,"\t",$target_read,"\t",$excluded_read,"\t",$translocation_read,"\n";
	}
}