#!/usr/bin/env perl
#### Author: Shengdi Li ####
#### 2022 Nov. 3 ####
#### Randomly shuffle a fasta genome file. CHR length and nucleotide composition remains same. ####
## Author: Shengdi Li
use strict;
use Bio::SeqIO;
use Getopt::Long;
use List::Util qw/shuffle/;

#### Load external data ####
GetOptions('h' => \my $display_help, 'in=s' => \my $input_fasta, 'out=s' => \my $output_fasta);

if($display_help)
{
	print "Command:\n\tperl shuffle_genome.pl [-h] -in INPUT_FASTA -out OUTPUT_FASTA\n";
	print "Function: \n\tGenerate a shuffled genome...\n\n";
	
	print "Usage:\n";
	print "\t-h\tPrint help info.\n";
	print "\t-in\tInput fasta.\n";
	print "\t-out\tOutput fasta. \n";
	
	exit 1;
}
if((!$input_fasta)||(!$output_fasta))
{
	print "Missing arguments.\n";
	print "Use -h for looking at usage...\n";
	exit 1;
}
open(OUT,">".$output_fasta);
my $genome_seqio=Bio::SeqIO->new(-file=>$input_fasta, -format=>"fasta");
while(my $chromosome_seq=$genome_seqio->next_seq)
{
	print OUT ">".$chromosome_seq->display_id."\n";
	my @tmp = split("",$chromosome_seq->seq);
	srand(149123);
	my $outseq=join("",shuffle(@tmp));
	while(length($outseq)>70)
	{
		print OUT substr($outseq,0,70),"\n";
		substr($outseq,0,70)="";
	}
	if($outseq)
	{
		print OUT $outseq,"\n";
	}
}