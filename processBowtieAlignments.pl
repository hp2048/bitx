#!/usr/bin/perl -w
use strict;

die "USAGE: $0 <infile sam> <outfile sam>" unless (scalar(@ARGV) == 2);
my ($infile, $outfile) = @ARGV;

open (EBAM, ">$outfile") or die "Can't open $outfile for writing.\n";

my %alignments = ();
my $total_alignments = 0;
my $valid_alignments = 0;
my $reads_aligned = 0;

open (RBAM, "<$infile") or die "Can't open $infile for reading.\n";
while (my $samline = <RBAM>){
  if ($samline =~ /^@/){
    print EBAM $samline;
    next;
  }
  chomp $samline;
  my @a = split ("\t", $samline);
  next if ($a[1] & 4); ###skip unaligned reads
  $total_alignments++;
  $samline =~ /CM:i:(\d+)/;
  $alignments{$a[0]}{$1}{$samline}="";
  if (keys %alignments == 2){
    $reads_aligned++;
    foreach my $read (keys %alignments){
      if ($read ne $a[0]){
        foreach my $colormismatches (sort {$a<=>$b} keys %{$alignments{$read}}){
          $valid_alignments += keys (%{$alignments{$read}{$colormismatches}});
          foreach my $aln (keys %{$alignments{$read}{$colormismatches}}){
            print EBAM "$aln\n";
          }
          last;
        }
      }
    }
    %alignments = ();
    $alignments{$a[0]}{$1}{$samline}="";
  }
}
###tackle the last read
$reads_aligned++;

foreach my $read (keys %alignments){
  foreach my $colormismatches (sort {$a<=>$b} keys %{$alignments{$read}}){
    $valid_alignments += keys (%{$alignments{$read}{$colormismatches}});
    foreach my $aln (keys %{$alignments{$read}{$colormismatches}}){
      print EBAM "$aln\n";
    }
    last;
  }
}

print "LOG: TotalAlignments=$total_alignments, ValidAlignments=$valid_alignments, ReadsAligned=$reads_aligned\n";

exit;
