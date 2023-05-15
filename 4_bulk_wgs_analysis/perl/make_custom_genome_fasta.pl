#!/usr/bin/env perl
#### Author: Shengdi Li ####
#### 2023 Feb. 9 ####
#### make a custom genome from selected chromosome names ####
#### analysis of individual WGS dataset with customized reference genomes ####

use Getopt::Long;
use Bio::SeqIO;
use strict;

GetOptions('h' => \my $display_help, 'in1=s' => \my $fasta1, 'in2=s' => \my $fasta2, 'out=s' =>\my $out_fasta, 'inc_chr=s' => \my $inc_chr_names);
if($display_help)
{
	print "Command: \n\tperl make_custom_genome_fasta.pl [-h] -in1 INPUT_FASTA_1 -in2 INPUT_FASTA_2 -out OUTPUT_FASTA -inc_chr INCLUDE_CHR_NAMES\n\n";
	print "Function: \n\tConcatenating sequences from two fasta files and select a list of chromosomes based on names...\n\n";
	print "Usage: \n";
	print "\t-h\tPrint help info.\n";
	print "\t-in1\tFirst input fasta.\n";
	print "\t-in2\tSecond input fasta.\n";
	print "\t-inc_chr\tNames of chromosomes to include separated by comma.\n";
	print "\t-out\tOutput fasta.\n";
	exit 1;
}

if((!$fasta1)||(!$fasta2)||(!$out_fasta)||(!$inc_chr_names))
{
	print "Missing arguments.\n";
	print "Use -h for looking at usage...\n";
	exit 1;
}

my @chr_names = split(",",$inc_chr_names);
my %chr_name;
foreach my $x(@chr_names)
{
	$chr_name{$x}=1;
}
open(OUT,">".$out_fasta);

my $seqio1=Bio::SeqIO->new(-file=>$fasta1,-format=>"fasta");
while(my $seq1=$seqio1->next_seq)
{
	if($chr_name{$seq1->display_id})
	{
		print OUT ">".$seq1->display_id."\n";
		print OUT $seq1->seq,"\n";
	}
}
my $seqio2=Bio::SeqIO->new(-file=>$fasta2,-format=>"fasta");
while(my $seq2=$seqio2->next_seq)
{
	if($chr_name{$seq2->display_id})
	{
		print OUT ">".$seq2->display_id."\n";
		print OUT $seq2->seq,"\n";
	}
}
