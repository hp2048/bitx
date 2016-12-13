#!/usr/bin/perl -w
use strict;
my ($infile, $outfile) = @ARGV;
open (IN, "<$infile") or die $!;
open (OUT, ">$outfile") or die $!;
my $seq = "";
my $header = "";
while (<IN>){
  chomp $_;
  if ($_ =~ />(.*)/){
	 if (length($seq)>0){
		for (my $i=-3;$i<=3;$i++){
		  if ($i != 0){
			 my $translation = get_translation($seq, $i);
			 print OUT ">$header".":$i\n";
			 print OUT "$translation\n";
		  }
		}
	 }
	 $header = $1;
	 $seq = "";
  }else{
	 $seq .= $_;
  }
}
if (length($seq)>0){
  for (my $i=-3;$i<=3;$i++){
	 if ($i != 0){
		my $translation = get_translation($seq, $i);
		print OUT ">$header".":$i\n";
		print OUT "$seq\n";
	 }
  }
}
sub get_translation {
  my %genetic_code  = (
  'TCA' => 'S',
  'TCC' => 'S',
  'TCG' => 'S',
  'TCT' => 'S',
  'TTC' => 'F',
  'TTT' => 'F',
  'TTA' => 'L',
  'TTG' => 'L',
  'TAC' => 'Y',
  'TAT' => 'Y',
  'TAA' => '*',
  'TAG' => '*',
  'TGC' => 'C',
  'TGT' => 'C',
  'TGA' => '*',
  'TGG' => 'W',
  'CTA' => 'L',
  'CTC' => 'L',
  'CTG' => 'L',
  'CTT' => 'L',
  'CCA' => 'P',
  'CAT' => 'H',
  'CAA' => 'Q',
  'CAG' => 'Q',
  'CGA' => 'R',
  'CGC' => 'R',
  'CGG' => 'R',
  'CGT' => 'R',
  'ATA' => 'I',
  'ATC' => 'I',
  'ATT' => 'I',
  'ATG' => 'M',
  'ACA' => 'T',
  'ACC' => 'T',
  'ACG' => 'T',
  'ACT' => 'T',
  'AAC' => 'N',
  'AAT' => 'N',
  'AAA' => 'K',
  'AAG' => 'K',
  'AGC' => 'S',
  'AGT' => 'S',
  'AGA' => 'R',
  'AGG' => 'R',
  'CCC' => 'P',
  'CCG' => 'P',
  'CCT' => 'P',
  'CAC' => 'H',
  'GTA' => 'V',
  'GTC' => 'V',
  'GTG' => 'V',
  'GTT' => 'V',
  'GCA' => 'A',
  'GCC' => 'A',
  'GCG' => 'A',
  'GCT' => 'A',
  'GAC' => 'D',
  'GAT' => 'D',
  'GAA' => 'E',
  'GAG' => 'E',
  'GGA' => 'G',
  'GGC' => 'G',
  'GGG' => 'G',
  'GGT' => 'G',
  );
  
  my $seq = shift;
  my $frame = shift;
  my $translation = "";
  if ($frame < 0){
	 $seq = reverse($seq);
	 $seq =~ tr/ACGT/TGCA/;
  }
  
  for (my $i=(abs($frame)-1);$i<length($seq)-3;$i+=3){
	 $translation .= $genetic_code{substr($seq, $i, 3)};
  }
  return $translation;
}
