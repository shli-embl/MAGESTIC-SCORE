#!/usr/bin/env perl
#### Author: Shengdi Li ####

#### 2022 April. 18 ####
#### Extraction of N-mer sequence centered at 20-bp gRNA targeted sequence on a given genome

use strict;
use Getopt::Long;
use Bio::SeqIO;
GetOptions('h' => \my $display_help, 'in=s' => \my $input_meta_file, 'out=s' =>\my $output_file, 'l=s' =>\my $seq_length, 'g=s' =>\my $genome);
if($display_help)
{
	print "Command: \n\tperl extract_nmer.pl [-h] -in INPUT_META_FILE -out OUTPUT_FILE -l LENGTH\n\n";
	print "Function: \n\tExtract n-mer sequence of a given gRNA target site...\n\n";
	print "Usage: \n";
	print "\t-h\tPrint help info.\n";
	print "\t-in\tInput list file containing the meta data of samples.\n";
	print "\t-out\tOutput file name.\n";
	print "\t-l\tLength of sequence to extract.\n";
	print "\t-g\tGenome sequence file.\n";
	exit 1;
}

if((!$input_meta_file)||(!$output_file)||(!$seq_length)||(!$genome))
{
	print "Missing arguments.\n";
	print "Use -h for looking at usage...\n";
	exit 1;
}

sub reversecomplement()
{
	my %complement=("A"=>"T","G"=>"C","T"=>"A","C"=>"G");
	my $in=$_[0];
	my $string = reverse $in;
	for(my $i=0;$i< length($string);$i++)
	{
		substr($string,$i,1)=$complement{substr($string,$i,1)}
	}
	return $string;
}
## input from STD
open(IN,$input_meta_file);
my $header = <IN>;
open(OUT,">".$output_file);
print OUT "#id\ttarget_sequence_".$seq_length."mer\n";
my $seqio=Bio::SeqIO->new(-file=>$genome);
my %chr_reference;
while(my $seq=$seqio->next_seq)
{
	$chr_reference{$seq->display_id}=$seq->seq;
}

while(my $read_line=<IN>)
{
	chomp $read_line;
	my ($id,$chr,$strand,$coord,$guide)=split("\t",$read_line);
	print OUT $id,"\t";
	if($strand eq "+")
	{
		print OUT substr($chr_reference{$chr},$coord-21-($seq_length-20)/2,$seq_length);
	}
	else
	{
		print OUT &reversecomplement(substr($chr_reference{$chr},$coord+2-($seq_length-20)/2,$seq_length));
	}
	print OUT "\n";
}