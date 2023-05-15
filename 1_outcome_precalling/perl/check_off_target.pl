#!/usr/bin/env perl
use strict;
use Bio::SeqIO;
use Getopt::Long;
GetOptions('h' => \my $display_help, 'in=s' => \my $input_meta_file, 'out=s' =>\my $output_file, 'vcf=s' =>\my $variant_file, 'w=s' =>\my $seq_length);
if($display_help)
{
	print "Command: \n\tperl check_off_target.pl [-h] -in INPUT_META_FILE -var VARIANT_FILE -w SEQUENCE_LENGTH -out OUTPUT_FILE \n\n";
	print "Function: \n\tFind best match from alignment and compute edit distance...\n\n";
	print "Usage: \n";
	print "\t-h\tPrint help info.\n";
	print "\t-in\tInput list file containing the meta data of samples.\n";
	print "\t-out\tOutput file name.\n";
	print "\t-vcf\tThe vcf file.\n";
	print "\t-w\tWindow size (length of region centered at PAM to be checked).\n";
	exit 1;
}

if((!$input_meta_file)||(!$variant_file)||(!$output_file)||(!$seq_length))
{
	print "Missing arguments.\n";
	print "Use -h for looking at usage...\n";
	exit 1;
}

open(VAR,$variant_file);
my %variants;
my @ids;
while(my $read_line=<VAR>)
{
	if($read_line =~ /^#CHROM/)
	{
		chomp $read_line;
		@ids = split("\t",$read_line);
#		print $ids[9],"\n";
	}
	elsif($read_line !~ /^#/)
	{
		my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@gts)=split("\t",$read_line);
		my @alts = split(",",$alt);
		for(my $i = 0; $i < (my $l = @gts);$i++)
		{
#			print "!!$i\t$l\n";
			my $gt = substr($gts[$i],0,1);
			if(($gt ne ".")&&($gt ne "0"))
			{
				$gt = $alts[$gt-1];
				$variants{$chr}{$pos}{$ids[$i+9]}=$chr."_".$pos."_".$ref."_to_".$gt;
#				print $i,"\t",$ids[1],"\n";
#				print $chr,"\t",$pos,"\t",$ids[$i+9],"\t",$gt,"\n";
			}
		}
	}
}
open(IN,$input_meta_file);
open(OUT,">".$output_file);
print OUT "#id\tchr\tcoord\tvariant\n";
my %guides;
while(my $read_line=<IN>)
{
	chomp $read_line;
	my ($sample_name,$guide,$chr,$coord,$match,$strand,$mismatch)=split("\t",$read_line);
	my $var;
	for(my $i = $coord - $seq_length/2; ($i<= $coord + $seq_length/2);$i++)
	{
		if($variants{$chr}{$i}{$sample_name})
		{
			if(!$var)
			{
				$var = $variants{$chr}{$i}{$sample_name}
			}
			else
			{
				$var = $var.";".$variants{$chr}{$i}{$sample_name}
			}
		}
	}
	print OUT $sample_name,"\t",$chr,"\t",$coord,"\t",$var,"\n";
}

