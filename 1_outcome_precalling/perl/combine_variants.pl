#!/usr/bin/env perl
#### Author: Shengdi Li ####
#### 2022 April. 21 ####
#### combine adjacent variants ####

use Getopt::Long;
use Bio::SeqIO;
use strict;

GetOptions('h' => \my $display_help, 'in=s' =>\my $input_file, 'out=s' =>\my $output_file);

if($display_help)
{
	print "Command:\n\tperl combine.pl [-h] -in INPUT_FILE -out OUTPUT_FILE \n\n";
	print "Function: \n\tCombine adjacent variants.\n\n";
	
	print "Usage:\n";
	print "\t-h\tPrint help info.\n";
	print "\t-in\tInput list of variants with sample ids.\n";
	print "\t-out\tOutput file.\n";
	
	exit 1;
}

if((!$output_file)||(!$input_file))
{
	print "Missing arguments.\n";
	print "Use -h for looking at usage...\n";
	exit 1;
}
sub min(){
	if($_[0]>$_[1])
	{
		return $_[1];
	}
	else
	{
		return $_[0];
	}
}
sub max(){
	if($_[0]<$_[1])
	{
		return $_[1];
	}
	else
	{
		return $_[0];
	}
}
sub distance(){
	my ($pos1,$pos2)=split("_",$_[0]);
	my ($pos3,$pos4)=split("_",$_[1]);
	return &max(0,&max($pos1,$pos3)-&min($pos2,$pos4));
}
sub merge(){
	my ($pos1,$pos2)=split("_",$_[0]);
	my ($pos3,$pos4)=split("_",$_[1]);
	return &min($pos1,$pos3)."_".&max($pos2,$pos4);
}

open(IN,$input_file);
open(OUT,">".$output_file);
my %variants;
while(my $read_line=<IN>)
{
	chomp $read_line;
	my ($sample_name,$variant)=split("\t",$read_line);
	my ($chr,$pos,$ref,$alt)=split("_",$variant);
	$variants{$sample_name}{$chr}{$pos."_".($pos + length($ref)-1)}=$variant;
}
foreach my $sample_name (keys %variants)
{
	my $collapse_end=0;
	while(!$collapse_end)
	{
		$collapse_end=1;
		foreach my $chr (keys %{$variants{$sample_name}})
		{
			foreach my $var1 (keys %{$variants{$sample_name}{$chr}})
			{
				foreach my $var2 (keys %{$variants{$sample_name}{$chr}})
				{
					if($var1 ne $var2)
					{
						if(&distance($var1,$var2)<50)
						{
							$variants{$sample_name}{$chr}{&merge($var1,$var2)}=join(",",($variants{$sample_name}{$chr}{$var1},$variants{$sample_name}{$chr}{$var2}));
							$collapse_end=0;
							delete $variants{$sample_name}{$chr}{$var1};
							delete $variants{$sample_name}{$chr}{$var2};
						}
					}
				}
			}
		}
	}
	foreach my $chr (keys %{$variants{$sample_name}})
	{
		foreach my $var1 (keys %{$variants{$sample_name}{$chr}})
		{
			my ($pos1,$pos2)=split("_",$var1);
			print OUT $sample_name,"\t",$chr,"\t",$pos1,"\t",$pos2,"\t",$variants{$sample_name}{$chr}{$var1},"\n";
		}
	}
}