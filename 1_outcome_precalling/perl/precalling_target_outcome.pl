#!/usr/bin/env perl
#### Author: Shengdi Li ####
#### 2022 Mar. 31 ####
#### Pre-call on-target edits/non-edits from VCF files ####
#### Output on-target status into: Edited, partially edited, unintended edits and unknown ####
#### Classification of intended vs. unintended variants are based on edit distances ####

use Getopt::Long;
use Bio::SeqIO;
use strict;

#start reading parameters
my $extension_donor_region=200;
my $non_edit_coverage_threshold=1;

GetOptions('h' => \my $display_help, 'in=s' => \my $input_meta_file, 'workdir=s' => \my $data_dir, 'out=s' =>\my $output_file_summary, 'donor_ext:i'=>\$extension_donor_region, 'NE_threshold:i'=>\$non_edit_coverage_threshold, 'genome=s'=>\my $reference_path);

if($display_help)
{
	print "Command:\n\tperl precalling_target_outcome.pl [-h] -in INPUT_LIST -workdir PATH_TO_DATA_FOLDER -out OUTPUT_SUMMARY [--donor_ext 200] [--NE_threshold 2]\n\n";
	print "Function: \n\tSummarizing on-target mutations from vcf files, classifying intended and unintended variants...\n\n";
	
	print "Usage:\n";
	print "\t-h\tPrint help info.\n";
	print "\t-in\tInput list file containing the meta data of samples.\n";
	print "\t-workdir\tThe path where WGS result sub-folders can be found (VCF and gVCF are required).\n";
	print "\t-out\tOutput file containing the summary of intended/unintended variants and edit status.\n";
	print "\t-donor_ext {200}\tNumber of basepairs up- and downstream the variant position to be considered for variant calling. The 2x{donor_ext} length sequence will be used to generate pseudo sequence for variant classifications.\n";
	print "\t-NE_threshold {1}\tThreshold of non-edited reference read depth for defining a non-edited sample; default requires 2 reads having the reference allele.\n";
	print "\t-genome\tReference genome file.\n";
	
	exit 1;
}

if((!$input_meta_file)||(!$data_dir)||(!$output_file_summary)||(!$reference_path))
{
	print "Missing arguments.\n";
	print "Use -h for looking at usage...\n";
	exit 1;
}

#define sub fuctions
sub levenshtein
{

    my ($s1, $s2) = @_;
    my ($len1, $len2) = (length $s1, length $s2);

    return $len2 if ($len1 == 0);
    return $len1 if ($len2 == 0);

    my %mat;

    for (my $i = 0; $i <= $len1; ++$i)
    {
        for (my $j = 0; $j <= $len2; ++$j)
        {
            $mat{$i}{$j} = 0;
            $mat{0}{$j} = $j;
        }

        $mat{$i}{0} = $i;
    }

    my @ar1 = split(//, $s1);
    my @ar2 = split(//, $s2);

    for (my $i = 1; $i <= $len1; ++$i)
    {
        for (my $j = 1; $j <= $len2; ++$j)
        {
            my $cost = ($ar1[$i-1] eq $ar2[$j-1]) ? 0 : 1;
            $mat{$i}{$j} = min([$mat{$i-1}{$j} + 1,
                                $mat{$i}{$j-1} + 1,
                                $mat{$i-1}{$j-1} + $cost]);
        }
    }
    return $mat{$len1}{$len2};
}

sub max
{
    my @list = @{$_[0]};
    my $max = $list[0];

    foreach my $i (@list)
    {
        $max = $i if ($i > $max);
    }

    return $max;
}
sub min
{
    my @list = @{$_[0]};
    my $min = $list[0];

    foreach my $i (@list)
    {
        $min = $i if ($i < $min);
    }

    return $min;
}

#loading yeast reference genome into MEM

my $seqio=Bio::SeqIO->new(-file=>$reference_path);
my %chr_reference;
while(my $seq=$seqio->next_seq)
{
	$chr_reference{$seq->display_id}=$seq->seq;
}

print "Finish loading the genome...\n";
print "Start processing samples...\n";
#start processing samples
open(INDEX,$input_meta_file);
my $header = <INDEX>;
open(OUT,">".$output_file_summary);
open(OUT2,">".$output_file_summary.".vars.txt");
print OUT "#sample_name\tediting_outcome\tdesigned_variants\tsynthetic_errors\tspontaneous_variants\n";
print OUT2 "#sample_name\tvariant\tpos_to_PAM\tvariant_type\toutcome\n";
while(my $read_line=<INDEX>)
{
	chomp $read_line;
	my ($sample_name,$path1,$path2,$chr,$strand,$pam_pos,$v_coord,$d_start,$d_end,$ref,$alt,$donor_wgs,$designed_v,$syn_v)=split("\t",$read_line);
	print "runing ",$sample_name,"...\t";
	#get the sequence of donor region from reference chromosome, donor_ext is the parameter to extend the region for scanning
	
	my $offset_long=&max([$d_start-$extension_donor_region-1,0]);
	my $offset_short=&max([$d_start-1,0]);
	my $region_end_long=&min([$d_end+$extension_donor_region-1,length($chr_reference{$chr})-1]);
	my $region_end_short=&min([$d_end-1,length($chr_reference{$chr})-1]);
	my $donor_wt_long=substr($chr_reference{$chr},$offset_long,$region_end_long-$offset_long+1);
	my $donor_wt_short=substr($chr_reference{$chr},$offset_short,$region_end_short-$offset_short+1);
	my $donor_mut_long=$donor_wt_long;
	my $donor_mut_short=$donor_wt_short;
	substr($donor_mut_long,$v_coord-$offset_long-1,length($ref))=$alt;
	substr($donor_mut_short,$v_coord-$offset_short-1,length($ref))=$alt;
	my %variants_sum;
	my @donor_wgs=split("",$donor_wt_short);
	my @designed_vs=split(",",$designed_v);
	my @syn_vs=split(",",$syn_v);
	foreach my $v(@designed_vs)
	{
		my ($tmpchr,$tmppos,$tmpref,$gap,$tmpalt)=split("_",$v);
		$variants_sum{$v}{"type"} = "designed";
		$variants_sum{$v}{"pos1"}=$tmppos;
		if(length($tmpref) eq length($tmpalt))
		{
			if($strand eq "+")
			{
				$variants_sum{$v}{"pos2"}=$pam_pos - $tmppos;
			}
			else
			{
				$variants_sum{$v}{"pos2"}=$tmppos - $pam_pos -2;
			}
		}
		elsif(length($tmpref) > length($tmpalt))
		{
			
			if($strand eq "+")
			{
				$variants_sum{$v}{"pos2"}=$pam_pos + length($tmpref) - $tmppos - 1;
			}
			else
			{
				$variants_sum{$v}{"pos2"}=$tmppos + 1 - $pam_pos -2;
			}
		}
		elsif(length($tmpref) < length($tmpalt))
		{
			if($strand eq "+")
			{
				$variants_sum{$v}{"pos2"}=$pam_pos - $tmppos - 1;
			}
			else
			{
				$variants_sum{$v}{"pos2"}=$tmppos + 1 - $pam_pos -2;
			}
		}
	}
	foreach my $v(@syn_vs)
	{
		my ($tmpchr,$tmppos,$tmpref,$gap,$tmpalt)=split("_",$v);
		$variants_sum{$v}{"type"} = "synthetic_error";
		$variants_sum{$v}{"pos1"}=$tmppos;
		if(length($tmpref) eq length($tmpalt))
		{
			
			if($strand eq "+")
			{
				$variants_sum{$v}{"pos2"}=$pam_pos - $tmppos;
			}
			else
			{
				$variants_sum{$v}{"pos2"}=$tmppos - $pam_pos -2;
			}
		}
		elsif(length($tmpref) > length($tmpalt))
		{
			if($strand eq "+")
			{
				$variants_sum{$v}{"pos2"}=$pam_pos + length($tmpref) - $tmppos - 1;
			}
			else
			{
				$variants_sum{$v}{"pos2"}=$tmppos + 1 - $pam_pos -2;
			}
		}
		elsif(length($tmpref) < length($tmpalt))
		{
			if($strand eq "+")
			{
				$variants_sum{$v}{"pos2"}=$pam_pos - $tmppos -1;
			}
			else
			{
				$variants_sum{$v}{"pos2"}=$tmppos + 1 - $pam_pos -2;
			}
		}
	}

	#scanning VCF file;
	open(VCF,$data_dir."/vcf/".$sample_name.".target.vcf");
	#flags
	my $flag_ifedit=0;
	my $flag_ifcovered=0;
	my (@pos_donor_called,%designed_variants,%synthetic_errors,%spontaneous_variants);
	
	while(my $read_line1=<VCF>)
	{
		if($read_line1!~/^#/)
		{
			chomp $read_line1;
			my ($v_chr,$v_pos,$gap,$v_ref,$v_alt)=split("\t",$read_line1);
			my $v_id=$v_chr."_".$v_pos."_".$v_ref."_to_".$v_alt;
			
			
			if((($v_pos-$offset_short-1)>=0)&&(($v_pos-$offset_short-1+length($v_ref))<length($donor_wt_short)))
			{
				my $donor_tmp_short=$donor_wt_short;
				substr($donor_tmp_short,$v_pos-$offset_short-1,length($v_ref))=$v_alt;
				#record positions that already called in vcf; so not double calling variants in gvcf
				my $i=length($v_ref);
				my $called_pos=$v_pos-$offset_long-1;
				while($i>0)
				{
					$pos_donor_called[$called_pos]="called";
					$called_pos++;
					$i--;
				}
				# calculate edit distances
				my $e1 = &levenshtein($donor_tmp_short,$donor_mut_short);
				my $e2 = &levenshtein($donor_wt_short,$donor_mut_short);
				my $e3 = &levenshtein($donor_tmp_short,$donor_wgs);
				my $e4 = &levenshtein($donor_wt_short,$donor_wgs);
				
				#if the variant has smaller leven distance than the wt sequence, compared with the design mutant, it's an intended variant
				if($e1 < $e2)
				{
					$designed_variants{$v_id}=$e2."->".$e1;
			#		$designed_variants_pos_sorting{$v_id}=$v_pos;
				}
				elsif(($e1 >= $e2)&&($e3 < $e4))
				{
					$synthetic_errors{$v_id}=$e2."->".$e1;
			#		$synthetic_errors_pos_sorting{$v_id}=$v_pos;
				}
				elsif(($e1 >= $e2)&&($e3 >= $e4))
				{
					$spontaneous_variants{$v_id}=$e2."->".$e1;
			#		$spontaneous_variants_pos_sorting{$v_id}=$v_pos;
				}
			}
			elsif((($v_pos-$offset_long-1)>=0)&&(($v_pos-$offset_long-1+length($v_ref))<length($donor_wt_long)))
			{
				my $donor_tmp_long=$donor_wt_long;
				substr($donor_tmp_long,$v_pos-$offset_long-1,length($v_ref))=$v_alt;
				#record positions that already called in vcf; so not double calling variants in gvcf
				my $i=length($v_ref);
				my $called_pos=$v_pos-$offset_long-1;
				while($i>0)
				{
					$pos_donor_called[$called_pos]="called";
					$called_pos++;
					$i--;
				}
				# calculate edit distances
				my $e1 = &levenshtein($donor_tmp_long,$donor_mut_long);
				my $e2 = &levenshtein($donor_wt_long,$donor_mut_long);
				if($e1 < $e2)
				{
					$designed_variants{$v_id}=$e2."->".$e1;
			#		$designed_variants_pos_sorting{$v_id}=$v_pos;
				}
				elsif($e1 >= $e2)
				{
					$spontaneous_variants{$v_id}=$e2."->".$e1;
			#		$spontaneous_variants_pos_sorting{$v_id}=$v_pos;
				}
			}
		}
	}
	
	#processing the gVCF file for extracting less confident intended variants (mark as LC)
	open(gVCF,$data_dir."/gvcf/".$sample_name.".target.g.vcf");
	while(my $read_line1=<gVCF>)
	{
		if($read_line1!~/^#/)
		{
			chomp $read_line1;
			my ($v_chr,$v_pos,$c3,$v_ref,$v_alts,$c6,$c7,$c8,$c9,$infos)=split("\t",$read_line1);
			my @v_alt=split(",",$v_alts);
			my @info=split(":",$infos);
			my @allele_depth=split(",",$info[1]);
			my $wt_read_depth=$allele_depth[0];
			
			#if the variant location has reference call over the defined threshold (default >=2 reads minimum), return a flag indicating the site is covered
			#version 2 take into account long reference allele
			if(($wt_read_depth >= $non_edit_coverage_threshold)&&(($v_pos >= $v_coord)&&($v_pos<=($v_coord + length($ref)-1))))
			{
				$flag_ifcovered++;
			}
			if(($wt_read_depth < $non_edit_coverage_threshold)&&(($v_pos >= $d_start)&&($v_pos <= $d_end)))
			{
				$donor_wgs[$v_pos - $d_start] = "N";
			}
			
			if((($v_pos-$offset_short-1)>=0)&&(($v_pos-$offset_short-1+length($v_ref))<length($donor_wt_short)))
			#if there is an alt call at the site , go through downstream check
			{
				if($v_alt[0] ne "<NON_REF>")
				{
					my $v_id=$v_chr."_".$v_pos."_".$v_ref."_to_".$v_alt[0];
					my $donor_tmp_short=$donor_wt_short;
					substr($donor_tmp_short,$v_pos-$offset_short-1,length($v_ref))=$v_alt[0];
					my $e1 = &levenshtein($donor_tmp_short,$donor_mut_short);
					my $e2 = &levenshtein($donor_wt_short,$donor_mut_short);
					my $e3 = &levenshtein($donor_tmp_short,$donor_wgs);
					my $e4 = &levenshtein($donor_wt_short,$donor_wgs);
					if($e1 < $e2)
					{
						#check if it's already called in vcf
						my $i=length($v_ref);
						my $called_pos=$v_pos-$offset_long-1;
						my $if_called=0;
						while(($i>0)&&(!$if_called))
						{
							if($pos_donor_called[$called_pos] eq "called")
							{
								$if_called=1;
							}
							$called_pos++;
							$i--;
						}
						#if it's a new mixture variant, report it as less confident call
						if(!$if_called)
						{
							$designed_variants{"(LC)".$v_id}=$e2."->".$e1;
					#		$designed_variants_pos_sorting{"(LC)".$v_id}=$v_pos;
						}
					
					}
					elsif(($e1 >= $e2)&&($e3 < $e4))
					{
						#if this is an synthetic error, check if it's already called in vcf
						my $i=length($v_ref);
						my $called_pos=$v_pos-$offset_long-1;
						my $if_called=0;
						while(($i>0)&&(!$if_called))
						{
							if($pos_donor_called[$called_pos] eq "called")
							{
								$if_called=1;
							}
							$called_pos++;
							$i--;
						}
						#if it's a new synthetic error, report it as less confident call
						if(!$if_called)
						{
							$synthetic_errors{"(LC)".$v_id}=$e2."->".$e1;
					#		$synthetic_errors_pos_sorting{"(LC)".$v_id}=$v_pos;
						}
					}
					elsif(($e1 >= $e2)&&($e3 >= $e4))
					{
						#if this is an spontaneous variant, check if it's already called in vcf
						my $i=length($v_ref);
						my $called_pos=$v_pos-$offset_long-1;
						my $if_called=0;
						while(($i>0)&&(!$if_called))
						{
							if($pos_donor_called[$called_pos] eq "called")
							{
								$if_called=1;
							}
							$called_pos++;
							$i--;
						}
						#if it's a new spontaneous variant, report it as less confident call
      						#Updated on June 19, 2025: stopped reporting LC spontaneous variants as they are easily confused with sequencing errors. LC variants can be useful for calling synthetic errors and designed variants, as the chance for sequencing error to match them is low. 
						#if(!$if_called)
						#{
						#	$spontaneous_variants{"(LC)".$v_id}=$e2."->".$e1;
						#	$spontaneous_variants_pos_sorting{"(LC)".$v_id}=$v_pos;
						#}
					}
					
				}
			}
			elsif((($v_pos-$offset_long-1)>=0)&&(($v_pos-$offset_long-1+length($v_ref))<length($donor_wt_long)))
			{
				if($v_alt[0] ne "<NON_REF>")
				{
					my $v_id=$v_chr."_".$v_pos."_".$v_ref."_to_".$v_alt[0];
					my $donor_tmp_long = $donor_wt_long;
					substr($donor_tmp_long,$v_pos-$offset_long-1,length($v_ref))=$v_alt[0];
					my $e1 = &levenshtein($donor_tmp_long,$donor_mut_long);
					my $e2 = &levenshtein($donor_wt_long,$donor_mut_long);
					if($e1 < $e2)
					{
						#check if it's already called in vcf
						my $i=length($v_ref);
						my $called_pos=$v_pos-$offset_long-1;
						my $if_called=0;
						while(($i>0)&&(!$if_called))
						{
							if($pos_donor_called[$called_pos] eq "called")
							{
								$if_called=1;
							}
							$called_pos++;
							$i--;
						}
						#if it's a new mixture variant, report it as less confident call
						if(!$if_called)
						{
						$designed_variants{"(LC)".$v_id}=$e2."->".$e1;
						#$designed_variants_pos_sorting{"(LC)".$v_id}=$v_pos;
						}
				
					}
					elsif($e1 >= $e2)
					{
						#if this is an spontaneous variant, check if it's already called in vcf
						my $i=length($v_ref);
						my $called_pos=$v_pos-$offset_long-1;
						my $if_called=0;
						while(($i>0)&&(!$if_called))
						{
							if($pos_donor_called[$called_pos] eq "called")
							{
								$if_called=1;
							}
							$called_pos++;
							$i--;
						}
						#if it's a new spontaneous variant, report it as less confident call
						#Updated on June 19, 2025: stopped reporting LC spontaneous variants as they are easily confused with sequencing errors. LC variants can be useful for calling synthetic errors and designed variants, as the chance for sequencing error to match them is low. 
						#if(!$if_called)
						#{
						#	$spontaneous_variants{"(LC)".$v_id}=$e2."->".$e1;
						#	print $v_id,"\n";
						#	#$spontaneous_variants_pos_sorting{"(LC)".$v_id}=$v_pos;
						#}
					}
					
				}
			}
		}
	}
	#combining all intended edits, distinguishing complete edits from partial edits
	my @donor_region_combine_array=split("",$donor_wt_long);
	foreach(keys %designed_variants)
	{
		my @spliter=split("_",$_);
		my $start_pos=$spliter[1]-$offset_long-1;
		my $i=length($spliter[2]);
		my $p1=$start_pos;
		while($i>0)
		{
			$donor_region_combine_array[$p1]="";
			$p1++;
			$i--;
		}
		$donor_region_combine_array[$start_pos]=$spliter[4];	
	}
	my $donor_region_combine=join("",@donor_region_combine_array);

	if(&levenshtein($donor_region_combine,$donor_mut_long) == 0)
	{
		$flag_ifedit="designed_variant_installed";
	}
	elsif(&levenshtein($donor_region_combine,$donor_mut_long) < &levenshtein($donor_wt_long,$donor_mut_long))
	{
		$flag_ifedit="designed_variant_partially_installed";
	}
	elsif($flag_ifcovered >= length($ref))
	{
		$flag_ifedit="designed_variant_unedited";
	}
	else
	{
		$flag_ifedit="unknown";
	}
	## process complex variant calling ##
	foreach(keys %designed_variants)
	{
		my @spliter=split("_",$_);
		my $start_pos=$spliter[1]-$offset_short-1;
		my $i=length($spliter[2]);
		my $p1=$start_pos;
		while($i>0)
		{
			$donor_wgs[$p1]="";
			$p1++;
			$i--;
		}
		$donor_wgs[$start_pos]=$spliter[4];	
	}
	foreach(keys %synthetic_errors)
	{
		my @spliter=split("_",$_);
		my $start_pos=$spliter[1]-$offset_short-1;
		my $i=length($spliter[2]);
		my $p1=$start_pos;
		while($i>0)
		{
			$donor_wgs[$p1]="";
			$p1++;
			$i--;
		}
		$donor_wgs[$start_pos]=$spliter[4];	
	}
	foreach(keys %spontaneous_variants)
	{
		my @spliter=split("_",$_);
		my $start_pos=$spliter[1]-$offset_short-1;
		if(($spliter[1] >= $d_start)&&(($spliter[1] + length($spliter[2]) -1) <= $d_end))
		{
			
			my $i=length($spliter[2]);
			my $p1=$start_pos;
			while($i>0)
			{
				$donor_wgs[$p1]="";
				$p1++;
				$i--;
			}
			$donor_wgs[$start_pos]=$spliter[4];	
		}
	}
	my $donor_wgs_seq = join("",@donor_wgs);
	
	### start processing each single variants on donor ###
	foreach my $xv(keys %variants_sum)
	{
		my $donor_tmp2 = $donor_wt_short;
		my ($tmpchr,$tmppos,$tmpref,$gap,$tmpalt)=split("_",$xv);
		substr($donor_tmp2,$tmppos-$d_start,length($tmpref))=$tmpalt;
		## if not covered
		if($donor_wgs[$tmppos-$d_start] eq "N")
		{
			$variants_sum{$xv}{"outcome"}="unknown";
		}
		else
		{
			## calculate e5,e6
			my $e5 = &levenshtein($donor_wgs_seq,$donor_wt_short);
			my $e6 = &levenshtein($donor_wgs_seq,$donor_tmp2);
			if($e5 > $e6)
			{
				$variants_sum{$xv}{"outcome"}="installed";
			}
			elsif($e5 < $e6)
			{
				$variants_sum{$xv}{"outcome"}="not_installed";
			}
			else
			{
				$variants_sum{$xv}{"outcome"}="unknown";
			}
		}
		print OUT2 $sample_name,"\t",$xv,"\t",$variants_sum{$xv}{"pos2"},"\t",$variants_sum{$xv}{"type"},"\t",$variants_sum{$xv}{"outcome"},"\n";		
	}
	#end of processing a sample; output line
	
	print OUT $sample_name,"\t";
	my $n_syn=(keys %synthetic_errors);
	my $n_sp=(keys %spontaneous_variants);
	print OUT $flag_ifedit;
	print $flag_ifedit;
	if($n_syn)
	{
		print OUT ";synthetic_error_installed";
		print ";synthetic_error_installed";
	}
	if($n_sp)
	{
		print OUT ";spontaneous_mutations";
		print ";spontaneous_mutations";
	}
	print OUT "\t";
	print "\n";
	
	my (@tmp_designed,@tmp_synthetic,@tmp_spontaneous);
	foreach(sort {$a cmp $b} (keys %designed_variants))
	{
		push @tmp_designed,$_;
	}
	print OUT join(",",@tmp_designed);
	print OUT "\t";
	foreach(sort {$a cmp $b} (keys %synthetic_errors))
	{
		push @tmp_synthetic,$_;
	}
	print OUT join(",",@tmp_synthetic);
	print OUT "\t";
	foreach(sort {$a cmp $b} (keys %spontaneous_variants))
	{
		push @tmp_spontaneous,$_;
	}
	print OUT join(",",@tmp_spontaneous);
	print OUT "\n";
}

