#!/usr/bin/env perl
#### Author: Shengdi Li ####
#### 2022 April. 21 ####
#### Extract singletons from joint call vcf ####

use Getopt::Long;
use Bio::SeqIO;
use strict;

GetOptions('h' => \my $display_help, 'in=s' =>\my $input_meta_file, 'out=s' =>\my $output_file, 'vcf=s'=>\my $vcf_file);

if($display_help)
{
	print "Command:\n\tperl generate_singleton_list.pl [-h] -in INPUT_LIST -out OUTPUT_FILE -vcf VCF_FILE \n\n";
	print "Function: \n\tExtract singletons from vcf.\n\n";
	
	print "Usage:\n";
	print "\t-h\tPrint help info.\n";
	print "\t-in\tInput sample meta sheet.\n";
	print "\t-out\tOutput file.\n";
	print "\t-vcf\tVCF file.\n";
	
	exit 1;
}

if((!$output_file)||(!$input_meta_file)||(!$vcf_file))
{
	print "Missing arguments.\n";
	print "Use -h for looking at usage...\n";
	exit 1;
}

open(IN,$input_meta_file);
open(OUT,">".$output_file);
my %PAM_sites;
while(my $read_line=<IN>)
{
	chomp $read_line;
	my ($sample_name,$path1,$path2,$chr,$strand,$pam_pos,$v_coord,$d_start,$d_end,$ref,$alt,$donor_wgs,$designed_v,$syn_v)=split("\t",$read_line);
	$PAM_sites{$sample_name}=$chr."_".$pam_pos;
}
### remake variant file removing individuals with missing data > 50% ###
open(VCF,$vcf_file);
open(TMP1,">".$output_file.".excludedInds.txt");
open(TMP2,">".$output_file.".excludedSites.txt");
my @header;
my %count_missing;
my $line_count=0;
print "Start filtering inds\n";
while(my $read_line=<VCF>)
{
	if($read_line=~/^#CHROM/)
	{
		chomp $read_line;
		@header=split("\t",$read_line);
	}
	elsif($read_line!~/#/)
	{
		chomp $read_line;
		my ($chr,$pos,$col3,$ref,$alt,$QUAL,$filter,$metrics,$fields,@GTs)=split("\t",$read_line);
		my $length=@GTs;
		my $i=0;
		while($i<$length)
		{
			my ($GT,$AD)=split(":",$GTs[$i]);
			if($GT eq ".")
			{
				$count_missing{$header[$i+9]}++;
			}
			
			$i++;
		}
		$line_count++;
	}
}
#### exclude samples with >50 % missing data ####
my $threshold1=$line_count * 0.5;
my %ex_ind;
foreach my $sample_name (keys %count_missing)
{
#	print $sample_name,"\t",$count_missing{$sample_name},"\t",$threshold1,"\n";
	if($count_missing{$sample_name}>$threshold1)
	{
		print TMP1 $sample_name,"\n";
		$ex_ind{$sample_name}++;
	}
}

#### exclude variants with >50 % missing data ####
my $sn=0;
open(VCF,$vcf_file);
my %ex_sites;
print "Start filtering sites\n";
while(my $read_line=<VCF>)
{
	if($read_line=~/^#CHROM/)
	{
		chomp $read_line;
		
		foreach my $x(split("\t",$read_line))
		{
			if($count_missing{$x}<=$threshold1)
			{
				$sn++;
			}
		}
#		print $sn,"\n";
	}
	elsif($read_line!~/#/)
	{
		chomp $read_line;
		my ($chr,$pos,$col3,$ref,$alt,$QUAL,$filter,$metrics,$fields,@GTs)=split("\t",$read_line);
		my $length=@GTs;
		my $i=0;
		#### exclude multi allelic loci ####
		#### exclude lowQual ####
#		print $filter;
		if(($alt!~/\,/)&&($filter eq "PASS"))
		{
			my $missing;
			while($i<$length)
			{
				my ($GT,$AD)=split(":",$GTs[$i]);
				if(($GT eq ".")&&(!$ex_ind{$header[$i+9]}))
				{
					$missing++;
				}
			
				$i++;
			}
#			print $missing,"\t",$sn * 0.5,"\n";
			if($missing > ($sn * 0.5))
			{
				print TMP2 $chr,"\t",$pos,"\t",$ref,"\t",$alt,"\thighMiss","\n";
				$ex_sites{$chr."_".$pos."_".$ref."_".$alt}++;
			}
		}
		else
		{
			print TMP2 $chr,"\t",$pos,"\t",$ref,"\t",$alt,"\tlowQual","\n";
			$ex_sites{$chr."_".$pos."_".$ref."_".$alt}++;
		}
	}
}
print "Start extracting singletons\n";
open(VCF,$vcf_file);
while(my $read_line=<VCF>)
{

	if($read_line!~/#/)
	{
		chomp $read_line;
		my ($chr,$pos,$col3,$ref,$alt,$QUAL,$filter,$metrics,$fields,@GTs)=split("\t",$read_line);
		my $length=@GTs;
		my $singleton;
		my $if_het=0;
		my $out;
		my $i=0;
		#### exclude multi allelic loci ####
		#### exclude lowQual ####
		if(($alt!~/\,/)&&($filter eq "PASS"))
		{
			while(($i<$length)&&($singleton ne "multi"))
			{
				my $sample_name = $header[$i+9];
				my ($GT,$AD)=split(":",$GTs[$i]);
				
				if($GT eq "1")
				{
					my ($AD1,$AD2)=split(",",$AD);
					if($AD1 > 0)
					{
						$if_het = 1;
					}
					if((!$singleton)&&(!$ex_ind{$header[$i+9]})&&(!$ex_sites{$chr."_".$pos."_".$ref."_".$alt}))
					{
						$singleton = $PAM_sites{$sample_name};
						$out=$out.$sample_name."\t".$chr."_".$pos."_".$ref."_".$alt."\n";
						
					}
					elsif(($singleton ne $PAM_sites{$sample_name})&&(!$ex_ind{$header[$i+9]})&&(!$ex_sites{$chr."_".$pos."_".$ref."_".$alt}))
					{
						$singleton="multi";

					}
					elsif((!$ex_ind{$header[$i+9]})&&(!$ex_sites{$chr."_".$pos."_".$ref."_".$alt}))
					{
						$out=$out.$sample_name."\t".$chr."_".$pos."_".$ref."_".$alt."\n";
					}
				}
				$i++;
			}
		}
		if(($singleton ne "multi")&&($singleton)&&(!$if_het))
		{
			#### exclude if overlap target +-200bp ####
			my ($vchr,$vpos)=split("_",$singleton);
			if(($vchr eq $chr)&&($vpos >= ($pos-200))&&($vpos <= ($pos+200)))
			{		
			}
			else
			{
				print OUT $out;
			}
		}
	}
}
