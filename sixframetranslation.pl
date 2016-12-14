#!/usr/bin/perl -w
use strict;
my ($infile, $outbase, $minLen, $maxLen) = @ARGV;
$minLen = (($minLen - 2) > 2) ? $minLen - 2 : 2;
$maxLen = (($maxLen - 2) > 2) ? $maxLen - 2 : 2;
open (ALL, ">$outbase.all.fasta") or die $!;
open (ORF, ">$outbase.orf.fasta") or die $!;
my $seqs = readFasta($infile);

foreach my $seqcounter (sort {$a<=>$b} keys %$seqs){
  foreach my $frame (1,2,3,-1,-2,-3){ 
	my $translation = get_translation($$seqs{$seqcounter}{'sequence'}, $frame);
	
	print ALL ">$$seqs{$seqcounter}{'id'}".":$frame\n";
	print ALL "$translation\n";
	
	while ($translation=~/(M\w{$minLen,$maxLen}\*)/g){
     my $orf   = $1;
     my $start = (pos($translation) - length($orf)) * 3 + 1;
     my $end   = pos($translation) * 3;
     $start = $start + abs($frame) - 1;
     $end   = $end   + abs($frame) - 1;
     
     print ORF ">$$seqs{$seqcounter}{'id'}".":$frame:$start:$end\n$orf\n";
	 pos($translation)=pos($translation)-length($orf)+1; 
	 ###this sets the pos so the next while iteration starts searching from after the detected methionine in this loop
	 ###eg. MAMAC* will be reported as MAMAC* and MAC* as well, if pos is not set then only MAMAC* gets reported
    } 
  }
}


sub readFasta {
  my ($f) = @_;
  my $counter = 0;
  my %seqs = ();
  open (FASTA, "<$f") or die $!;
  my $id = "";
  while (<FASTA>){
    chomp $_;
    if ($_ =~ />\s*(\S+)\s*(.*)/){
      $id = $1;
      my $description = $2;
 	  $counter++;
      $seqs{$counter}{'id'} = $id;
      $seqs{$counter}{'description'} = $description;
    }
    else{
      $seqs{$counter}{'sequence'} .= $_;
    }
  }
  close FASTA;
  return \%seqs;
}

sub get_translation {
  my %genetic_code  = (
  'GCA' => 'A',  'GCC' => 'A',  'GCG' => 'A',  'GCT' => 'A',
  'TGC' => 'C',  'TGT' => 'C',
  'GAC' => 'D',  'GAT' => 'D',
  'GAA' => 'E',  'GAG' => 'E',
  'TTC' => 'F',  'TTT' => 'F',
  'GGA' => 'G',  'GGC' => 'G',  'GGG' => 'G',  'GGT' => 'G',
  'CAT' => 'H',  'CAC' => 'H',
  'ATA' => 'I',  'ATC' => 'I',  'ATT' => 'I',
  'AAA' => 'K',  'AAG' => 'K',
  'TTA' => 'L',  'TTG' => 'L',  'CTA' => 'L',  'CTC' => 'L',  'CTG' => 'L',  'CTT' => 'L',
  'ATG' => 'M',
  'AAC' => 'N',  'AAT' => 'N',
  'CCA' => 'P',  'CCC' => 'P',  'CCG' => 'P',  'CCT' => 'P',
  'CAA' => 'Q',  'CAG' => 'Q',
  'CGA' => 'R',  'CGC' => 'R',  'CGG' => 'R',  'CGT' => 'R',  'AGA' => 'R',  'AGG' => 'R',
  'TCA' => 'S',  'TCC' => 'S',  'TCG' => 'S',  'TCT' => 'S',  'AGC' => 'S',  'AGT' => 'S',
  'ACA' => 'T',  'ACC' => 'T',  'ACG' => 'T',  'ACT' => 'T',
  'GTA' => 'V',  'GTC' => 'V',  'GTG' => 'V',  'GTT' => 'V',
  'TGG' => 'W',
  'TAC' => 'Y',	 'TAT' => 'Y',
  'TAA' => '*',  'TAG' => '*',  'TGA' => '*',
  );
  
  my $seq = shift;
  my $frame = shift;
  my $translation = "";
  $seq = uc($seq);
  if ($frame < 0){
	 $seq = reverse($seq);
	 $seq =~ tr/ACGT/TGCA/;
  }
  
  for (my $i=(abs($frame)-1);$i<=length($seq)-3;$i+=3){
	 $translation .= $genetic_code{substr($seq, $i, 3)};
  }
  return $translation;
}
