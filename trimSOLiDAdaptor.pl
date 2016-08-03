#!/usr/bin/perl -w
use strict;

die "USAGE: $0 csfastain csqualityin csfastaout csqualityout minimumlength minimumquality\n" unless (scalar(@ARGV) == 6);
my ($csin, $cqin, $trimmedcsfasta, $trimmedqvfasta, $minlength, $minqv) = @ARGV;

print "LOG: Started processing for adaptor removal and quality filter at ". (`date`);
######create file handles for the output
open (OUTCS, ">$trimmedcsfasta") or die $!;
open (OUTQV, ">$trimmedqvfasta") or die $!;
######Read in colors from the file
if ($csin =~ /\.gz$/){
  open (COLORS, "gunzip -c $csin |") or die $!;
}else{
  open (COLORS, "<$csin") or die $!;
}
if ($cqin =~ /\.gz$/){
  open (QUAL, "gunzip -c $cqin |") or die $!;
}else{
  open (QUAL, "<$cqin") or die $!;
}
####chew up lines before the header line in quality files
my ($quality_header, $quality_seq, $color_header, $color_seq);
while ($quality_header = <QUAL>){
  if ((substr $quality_header, 0, 1) ne "#"){
    $quality_seq = <QUAL>;
    last;
  }else{
    next;
  }
}
my $tags_counts = 0;
my $retained_tags_adapter = 0;
my $retained_tags_quality = 0;

while (my $colorline = <COLORS>){
  ###skip comment lines
  next if ((substr $colorline, 0, 1) eq "#");
  if ((substr $colorline, 0, 1) eq ">"){
    $color_header = $colorline;
    $color_seq    = <COLORS>;
    ####chew all the quality sequences that do not match the current color sequence
    while ($quality_header ne $color_header){
      $quality_header = <QUAL>;
      $quality_seq = <QUAL>;
    }
    $tags_counts++;
    chomp $color_seq;
    ######look for adaptor sequence
    my ($adaptor_begin, $adaptor_score, $adaptor_mismatches) = findAdaptor($color_seq);
    my $retained_seq = 0;
    if ($adaptor_begin > $minlength){
      ####should contain a distinguishable adaptor sequence, divide by four to keep error rate below 25%
      if ($adaptor_mismatches <= ((length($color_seq) - $adaptor_begin) / 4)){
        $retained_seq = substr($color_seq,0,$adaptor_begin) if (substr($color_seq,0,$adaptor_begin) !~ /\./);
        #$retained_tags{$line} = substr($color_seq,0,$adaptor_begin) if (substr($color_seq,0,$adaptor_begin) !~ /\./);
      }
      else{
        $retained_seq = $color_seq if ($color_seq !~ /\./);
        #$retained_tags{$line} = $cspace if ($color_seq !~ /\./);
      }
    }
    else{
      ###retain the whole sequence a distinguishable adaptor sequence is not found
      if ($adaptor_mismatches > ((length($color_seq) - $adaptor_begin) / 4)){
        $retained_seq = $color_seq if ($color_seq !~ /\./);
        #$retained_tags{$line} = $cspace if ($cspace !~ /\./);
      }
    }
    
    ####check for the quality values now
    if ($retained_seq ne "0"){
      $retained_tags_adapter++;
      chomp $quality_seq;
      my @fullQ = split(" ", $quality_seq);
      if (Mean(\@fullQ) >= $minqv){
        ###need to change to the quality values of the retained part only in the future release
        $retained_tags_quality++;
        print OUTCS "$color_header$retained_seq\n";
        print OUTQV "$color_header".join(" ", @fullQ[0..(length($retained_seq) - 2)])."\n";
      }
    }
  }
}
close COLORS;

print "LOG: $csin, TotalTags=$tags_counts, RetainedTagCountsAdapter=$retained_tags_adapter, RetainedTagCountsQuality=$retained_tags_quality\n";
print "LOG: Finished processing for adaptor removal and quality filter at ". (`date`);


#######subroutine to find start of the adaptor sequence
sub findAdaptor {
  my $tag_seq = shift;
  my $adaptor_seq = "330201030313112312"; ####first C is trimmed off
  ####Change this adaptor sequence trimming off the first color since the first color could be any of four possible colors
  ####strategy
  #### T1333002212102132203302010303131123  color-space read
  ####                   *================  adaptor sequence
  my %score4positions = ();
  my $tlen = length($tag_seq);
  my $alen = length($adaptor_seq);
  if ((index($tag_seq,$adaptor_seq) - 1) >= 0){ ###see if the entire adaptor is found as it is in the tag
    return (index($tag_seq,$adaptor_seq) - 1), ($alen + 1), 0;
  }else{  ###identify the high scoring position
    my $maxscore = -10000000;
    my $mismatches = 0;
    my ($adaptor_begin, $adaptor_mismatches);
    for (my $start = 2; $start < $tlen; $start++){
      my $t = $start;
      my $score = 0;
      $mismatches = 0;
      for (my $i = 0; $i < $alen; $i++){
        if (substr($tag_seq, $t, 1) ne substr($adaptor_seq, $i, 1)){
          $score += -1;
          $mismatches++;
          last if ($mismatches > $alen/4); ####this will save time if no. of mismatches exceeds this by not looping over the entire adaptor sequence
        }else{
          $score++;
        }
        $t++;
        last if ($t == $tlen);
      }
      if ($score > $maxscore){
        $maxscore = $score;
        $adaptor_begin = $start - 1;
        $adaptor_mismatches = $mismatches;
      }
    }
    return $adaptor_begin, $maxscore, $adaptor_mismatches;
  }
}

sub Mean {
  my $arrayref = shift;
  return undef unless defined $arrayref && @$arrayref > 0;
  my $result;
  foreach ( @$arrayref ) { $result += $_ }
  return $result / @$arrayref;
}


