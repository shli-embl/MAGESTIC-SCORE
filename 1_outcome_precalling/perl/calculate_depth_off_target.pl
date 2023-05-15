#!/usr/bin/env perl
use strict;
use Getopt::Long;
GetOptions('h' => \my $display_help, 'in=s' => \my $input_meta_file, 'out=s' =>\my $output_file, 'workdir=s' =>\my $workdir);
if($display_help)
{
	print "Command: \n\tperl calculate_depth_off_target.pl [-h] -in INPUT_META_FILE -workdir WORKDIR -out OUTPUT_FILE \n\n";
	print "Function: \n\tGet mapping coverage for the potential off-target positions...\n\n";
	print "Usage: \n";
	print "\t-h\tPrint help info.\n";
	print "\t-in\tInput list file containing the meta data of samples.\n";
	print "\t-out\tOutput file name.\n";
	print "\t-workdir\tWork directory containing BAM alignments.\n";
	exit 1;
}

if((!$input_meta_file)||(!$output_file)||(!$workdir))
{
	print "Missing arguments.\n";
	print "Use -h for looking at usage...\n";
	exit 1;
}

open(IN,$input_meta_file);
my $header=<IN>;
open(OUT,">".$output_file);
print OUT "#id\tchr\tcoord\tvariant\tmin_depth\tmean_depth\n";
while(my $read_line=<IN>)
{
	chomp $read_line;
	my ($sample_name,$guide,$chr,$coord,$match,$strand,$mismatch)=split("\t",$read_line);
	if($strand eq "+")
	{
		system("samtools depth -a -r ".$chr.":".($coord + 1)."-".($coord + 23)." ".$workdir."/bam/".$sample_name.".clean.bam > ".$output_file.".tmp.".$sample_name);
	}
	else
	{
		system("samtools depth -a -r ".$chr.":".($coord - 2)."-".($coord + 20)." ".$workdir."/bam/".$sample_name.".clean.bam > ".$output_file.".tmp.".$sample_name);
	}
	open(TMP,$output_file.".tmp.".$sample_name);
	my $min_depth =1000;
	my $sum_depth =0;
	my $n_site = 0;
	while(my $read_line2=<TMP>)
	{
		chomp $read_line2;
		my ($chr2,$coord2,$depth) = split("\t",$read_line2);
		if($depth < $min_depth)
		{
			$min_depth = $depth;
		}
		$sum_depth = $sum_depth + $depth;
		$n_site++;
	}
	system("rm ".$output_file.".tmp.".$sample_name);
	print OUT $sample_name,"\t",$chr,"\t",$coord,"\t",$min_depth,"\t",sprintf("%.2f",$sum_depth/$n_site),"\n";
}
