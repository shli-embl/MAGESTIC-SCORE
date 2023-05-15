#!/usr/bin/env perl
#### Author: Shengdi Li ####
#### 2022 April. 19 ####
#### Calculate the distance from pam site to telomere or centromere, in bp or in %

use strict;
use Getopt::Long;
GetOptions('h' => \my $display_help, 'in=s' => \my $input_meta_file, 'out=s' => \my $outfile, 'coord=s' => \my $coordfile);

if($display_help)
{
	print "Command:\n\tperl dist_centromere_telomere.pl [-h] -in INPUT_META_LIST -out OUTFILE -coord COORDINATES_CENTROMERE\n";
	print "Function: \n\tCalculate distances from pam sites to telomere/centromere...\n\n";
	
	print "Usage:\n";
	print "\t-h\tPrint help info.\n";
	print "\t-in\tInput sample meta sheet.\n";
	print "\t-out\tOutput file path. \n";
	print "\t-coord\tFiles containing centromere coordinates and chromosome lengths.\n";
	
	exit 1;
}
if((!$input_meta_file)||(!$outfile)||(!$coordfile))
{
	print "Missing arguments.\n";
	print "Use -h for looking at usage...\n";
	exit 1;
}

sub min
{
	if($_[0]>$_[1])
	{
		return $_[1]
	}
	else
	{
		return $_[0]
	}
}

open(REF,$coordfile);
my %chr_length;
my %cent_coord;
my $header=<REF>;
while(my $read_line=<REF>)
{
	chomp $read_line;
	my ($chr,$length,$cent_left,$cent_right)=split("\t",$read_line);
	$chr_length{$chr}=$length;
	$cent_coord{$chr}{"L"}=$cent_left;
	$cent_coord{$chr}{"R"}=$cent_right;
}

open(IN,$input_meta_file);
my $header=<IN>;
open(OUT,">".$outfile);
print OUT "#id\tbp_to_telomere\tbp_to_centromere\tnorm_distance_to_telomere\n";
while(my $read_line=<IN>)
{
	chomp $read_line;
	my ($id,$chr,$strand,$coord,$guide)=split("\t",$read_line);
	my $dist_to_tel;
	my $dist_to_cen;
	if($coord < $cent_coord{$chr}{"L"})
	{
		$dist_to_tel = $coord;
		$dist_to_cen = $cent_coord{$chr}{"L"} - $coord;
	}
	elsif($coord > $cent_coord{$chr}{"R"})
	{
		$dist_to_tel = $chr_length{$chr}-$coord;
		$dist_to_cen = $coord - $cent_coord{$chr}{"R"};
	}
	else
	{
		$dist_to_tel = $coord;
		$dist_to_cen = 0;
	}
	print OUT $id,"\t",$dist_to_tel,"\t",$dist_to_cen,"\t";
	print OUT sprintf("%.3f",$dist_to_tel/($dist_to_tel + $dist_to_cen)),"\n";

}