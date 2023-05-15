#!/usr/bin/env perl
#### Author: Shengdi Li ####
#### 2022 Nov. 2 ####
#### This is a script to prepare for input fasta sequence for ir program. The file contains N-nucleotide sequence centered at the PAM location of a target site. ####
## Author: Shengdi Li
use strict;
use Bio::SeqIO;
use Getopt::Long;
GetOptions('h' => \my $display_help, 'chr=s' => \my $chr, 'coord=s' => \my $coord, 'out=s' =>\my $out, 'id=s' =>\my $id, 'genome=s' =>\my $genome, 'window=s' =>\my $window);
if($display_help)
{
	print "Command: \n\tperl extract_fasta_target_window.pl [-h] -chr CHROMOSOME_ID -coord COORDINATE -out OUTPUT_FILE -window WINDOW_SIZE -genome REFERENCE_FASTA -id SAMPLE_ID\n\n";
	print "Function: \n\tCalculate Ir score that represent local sequence repetativeness.\n\n";
	print "Usage: \n";
	print "\t-h\tPrint help info.\n";
	print "\t-chr\tFrom which chromosome id the sequence should be extracted.\n";
	print "\t-coord\tFrom which coordinate on the chromosome id the sequence should be extracted.\n";
	print "\t-window\tLength of the sequence centered at the coordinate should be extracted.\n";
	print "\t-out\tOutput fasta file name.\n";
	print "\t-id\tName of the sequence.\n";
	print "\t-g\tGenome sequence file.\n";
	exit 1;
}

if((!$chr)||(!$out)||(!$window)||(!$genome)||(!$id)||(!$coord))
{
	print "Missing arguments.\n";
	print "Use -h for looking at usage...\n";
	exit 1;
}
sub max()
{
	my $a=$_[0];
	my $b=$_[1];
	if($a > $b)
	{
		return $a;
	}
	else
	{
		return $b;
	}
}

sub min()
{
	my $a=$_[0];
	my $b=$_[1];
	if($a > $b)
	{
		return $b;
	}
	else
	{
		return $a;
	}
}

my $seqio=Bio::SeqIO->new(-file=>$genome,-format=>"fasta");
my $chr_seq;
while(my $seq=$seqio->next_seq)
{
	if($seq->display_id eq $chr)
	{
		$chr_seq=$seq->seq;
		last;
	}
}
open(OUT,">".$out);
print OUT ">".$id."\n";
print OUT substr($chr_seq,&max(0,$coord-$window/2-1),&min(length($chr_seq),$coord+$window/2-1)-&max(0,$coord-$window/2-1)),"\n";
