#!/usr/bin/env perl
use strict;
use Bio::SeqIO;
use Getopt::Long;
GetOptions('h' => \my $display_help, 'in=s' => \my $input_meta_file, 'out=s' =>\my $output_file, 'var=s' =>\my $variant_file, 'w=s' =>\my $seq_length, 'g=s' =>\my $genome);
if($display_help)
{
	print "Command: \n\tperl edit_distance_at_variant.pl [-h] -in INPUT_META_FILE -var VARIANT_FILE -w SEQUENCE_LENGTH -out OUTPUT_FILE -g GENOME_FASTA\n\n";
	print "Function: \n\tFind best match from alignment and compute edit distance...\n\n";
	print "Usage: \n";
	print "\t-h\tPrint help info.\n";
	print "\t-in\tInput list file containing the meta data of samples.\n";
	print "\t-out\tOutput file name.\n";
	print "\t-var\tList of variants identified in each sample.\n";
	print "\t-w\tWindow size (length of region centered at PAM to be checked).\n";
	print "\t-g\tGenome sequence fasta file.\n";
	exit 1;
}

if((!$input_meta_file)||(!$variant_file)||(!$output_file)||(!$seq_length)||(!$genome))
{
	print "Missing arguments.\n";
	print "Use -h for looking at usage...\n";
	exit 1;
}
sub offtarget
{		
	my $seq1=$_[0];
	my $seq2=$_[1];
	print "seq1: ",$seq1,"\n";
	print "seq2: ",$seq2,"\n";
	$seq2=$seq2."NRG";
	my %scorematrix;
	my @nucl=qw{A T C G N};
	foreach my $x(@nucl)
	{
		foreach my $y(@nucl)
		{
			if(($x eq "N")||($y eq "N"))
			{
				$scorematrix{$x}{$y}{"distant_to_PAM"}=5;
				$scorematrix{$x}{$y}{"close_to_PAM"}=5;
				$scorematrix{$x}{$y}{"PAM"}=5;
			}
			elsif($x eq $y)
			{
#score for a match
				$scorematrix{$x}{$y}{"distant_to_PAM"}=5;
				$scorematrix{$x}{$y}{"close_to_PAM"}=5;
				$scorematrix{$x}{$y}{"PAM"}=5;
			}
			else
			{
		#score for a mismatch		
				$scorematrix{$x}{$y}{"distant_to_PAM"}=-5;
				$scorematrix{$x}{$y}{"close_to_PAM"}=-5;
				$scorematrix{$x}{$y}{"PAM"}=-10000;
			}
			
		}
			
		if(($x eq "A")||($x eq "G"))
		{
			$scorematrix{$x}{"R"}{"distant_to_PAM"}=5;
			$scorematrix{$x}{"R"}{"close_to_PAM"}=5;
			$scorematrix{$x}{"R"}{"PAM"}=5;
			$scorematrix{"R"}{$x}{"distant_to_PAM"}=5;
			$scorematrix{"R"}{$x}{"close_to_PAM"}=5;
			$scorematrix{"R"}{$x}{"PAM"}=5;
		}
		else
		{
			$scorematrix{$x}{"R"}{"distant_to_PAM"}=-5;
			$scorematrix{$x}{"R"}{"close_to_PAM"}=-5;
			$scorematrix{$x}{"R"}{"PAM"}=-5;
			$scorematrix{"R"}{$x}{"distant_to_PAM"}=-5;
			$scorematrix{"R"}{$x}{"close_to_PAM"}=-5;
			$scorematrix{"R"}{$x}{"PAM"}=-10000;
		}
	}
#score for a gap
	my %gap;
	$gap{"distant_to_PAM"}=-15;
	$gap{"close_to_PAM"}=-15;
	$gap{"PAM"}=-10000;

#initiation
	my @matrix; 
	$matrix[0][0]{"score"}   = 0; 
	$matrix[0][0]{"pointer"} = "none"; 

	for (my $i = 1; $i <= length($seq1); $i++) 
	{
		 $matrix[$i][0]{"score"}   = 0;
		 $matrix[$i][0]{"pointer"} = "Left"; 
	}  
	for(my $j = 1; $j <= length($seq2); $j++) 
	{
		my $distance_to_PAM;
		if($j>(length($seq2)-11))
		{
			$distance_to_PAM="close_to_PAM";
			if($j>(length($seq2)-3))
			{
				$distance_to_PAM="PAM";
			}
			
		}
		else
		{
			$distance_to_PAM="distant_to_PAM";
		}
		$matrix[0][$j]{"score"}   = $matrix[0][$j-1]{score}+$gap{$distance_to_PAM};
		$matrix[0][$j]{"pointer"} = "Right";
	} 

#fill the matrix
	for(my $i = 1; $i <= length($seq1); $i++) 
	{
		 for(my $j = 1; $j <= length($seq2); $j++) 
		 {
				 	
			my $distance_to_PAM;
			if($j>(length($seq2)-11))
			{
				$distance_to_PAM="close_to_PAM";
				if($j>(length($seq2)-3))
				{
					$distance_to_PAM="PAM";
				}
				
			}
			else
			{
				$distance_to_PAM="distant_to_PAM";
			}
			$matrix[$i][$j]{"score"}=$matrix[$i-1][$j-1]{"score"}+$scorematrix{substr($seq1,$i-1,1)}{substr($seq2,$j-1,1)}{$distance_to_PAM};

			$matrix[$i][$j]{"pointer"}="Match";
			if($j==length($seq2))
			{
				if($matrix[$i][$j]{"score"}<$matrix[$i-1][$j]{"score"})
				{
					$matrix[$i][$j]{"score"}=$matrix[$i-1][$j]{"score"};
					$matrix[$i][$j]{"pointer"}="Left";
				}
			}
			else
			{
				if($matrix[$i][$j]{"score"}<($matrix[$i-1][$j]{"score"}+$gap{$distance_to_PAM}))
				{
					$matrix[$i][$j]{"score"}=$matrix[$i-1][$j]{"score"}+$gap{$distance_to_PAM};
					$matrix[$i][$j]{"pointer"}="Left";
				}
			}
			
			if($matrix[$i][$j]{"score"}<($matrix[$i][$j-1]{"score"}+$gap{$distance_to_PAM}))
			{
				$matrix[$i][$j]{"score"}=$matrix[$i][$j-1]{"score"}+$gap{$distance_to_PAM};
				$matrix[$i][$j]{"pointer"}="Right";
			}
			 
		}
	}
#find best local alignment
	my $maxscore=$matrix[length($seq1)][length($seq2)]{score};

#trace back
	my $align1;
	my $align2;
	my $i=length($seq1);
	my $j=length($seq2);
		
#count mismatches,deletions
	while(1) 
	{
		last if $matrix[$i][$j]{"pointer"} eq "none";
		if($matrix[$i][$j]{"pointer"} eq "Match") 
		{
			$align1 .= substr($seq1, $i-1, 1);
			$align2 .= substr($seq2, $j-1, 1);
			$i--; $j--;				
		}
		elsif($matrix[$i][$j]{"pointer"} eq "Left") 
		{
			$align1 .= substr($seq1, $i-1, 1);
			$align2 .= "-";
			$i--;
		}
		elsif($matrix[$i][$j]{"pointer"} eq "Right")
		{
			$align1 .= "-";
			$align2 .= substr($seq2, $j-1, 1);
			$j--;
		}
	}  
	$align1=reverse $align1;
	$align2=reverse $align2;
	my $tmpaln1=$align1;
	my $tmpaln2=$align2;
	while(substr($tmpaln2,0,1) eq "-")
	{
		substr($tmpaln2,0,1)="";
		substr($tmpaln1,0,1)="";
	}
	while(substr($tmpaln2,length($tmpaln2)-1,1) eq "-")
	{
		substr($tmpaln1,length($tmpaln1)-1,1)="";
		substr($tmpaln2,length($tmpaln2)-1,1)="";
	}
	my $mm=0;
	for(my $i = 0; $i<length($tmpaln1);$i++)
	{
		if((substr($tmpaln1,$i,1) eq "-")||(substr($tmpaln2,$i,1) eq "-")||($scorematrix{substr($tmpaln1,$i,1)}{substr($tmpaln2,$i,1)}{"PAM"} < 0))
		{
			$mm++;
		}
	}
	return ($maxscore,$align1, $align2,$mm);
}

sub reversecomplement
{
	my $input=$_[0];
	my %rc=("A"=>"T","T"=>"A","C"=>"G","G"=>"C","N"=>"N");

	for(my $p=0;$p<length($input);$p++)
	{
		substr($input,$p,1)=$rc{substr($input,$p,1)};
			
	}
	$input=reverse $input;
	return $input;
}

my $seqio=Bio::SeqIO->new(-file=>$genome);
my %chr_reference;
while(my $seq=$seqio->next_seq)
{
	$chr_reference{$seq->display_id}=$seq->seq;
}

open(IN,$input_meta_file);
my %guides;
while(my $read_line=<IN>)
{
	chomp $read_line;
	my ($sample_name,$guide)=split("\t",$read_line);
	$guides{$sample_name}=$guide;
}
open(VAR,$variant_file);
open(OUT,">".$output_file);
while(my $read_line=<VAR>)
{
	chomp $read_line;
	my ($sample_name,$chr,$start,$end,$variant)=split("\t",$read_line);
	my $seq=substr($chr_reference{$chr},$start-$seq_length/2,$end-$start+$seq_length);
	my ($maxscore,$align1, $align2, $mm)=&offtarget($seq,$guides{$sample_name});
	my ($maxscore1,$align3, $align4, $mm1)=&offtarget(&reversecomplement($seq),$guides{$sample_name});
	if($maxscore>$maxscore1)
	{
		print OUT $read_line,"\t",$mm,"\t",$align1,"\t",$align2,"\n";
	}
	else
	{
		print OUT $read_line,"\t",$mm1,"\t",$align3,"\t",$align4,"\n";
	}
}

