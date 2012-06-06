#!/usr/bin/env perl

#
# This is a reference implementation of the *NAIVE* row-based partitioning
# of ADCIRC grids, which to minimize the number of communicating 
# interfaces between domains.  The balancing of grid nodes is accomplished
# using a k-step shortest path formulation.
#

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;

my $F14 = "./fort.14";
my $MINCOUNT = 5000;
my $DIVISOR = undef;

GetOptions("f=s" => \$F14, 
           "m=i" => \$MINCOUNT,
           "d=i" => \$DIVISOR,); # overrides -m, or $MINCOUNT

my $NE = 0;
my $NN = 0;
# number of significant figures to use when counting nodes as being in the same longitude
my $precision = 5; 
my %bins = ();

# creates the arithmetic mean of a given population
sub mean {
  my $bins_ref  = shift; # the main bin hash
  my $group_ref = shift; # the group members, a subset of the bin hash
  my $sample_size = 0;
  my $sum = 0;
  foreach my $group_member (@{$group_ref}) {
    $sum += $group_member * $bins_ref->{$group_member}; # !! Assumes that all group member scores are >0
    $sample_size += $bins_ref->{$group_member};
  }
  my $mean = $sum / $sample_size;
  return $mean;
}

sub mean_a {
  my $group_ref = shift; # the group members, a subset of the bin hash
  my $sample_size = @{$group_ref};
  my $sum = 0;
  foreach my $group_member (@{$group_ref}) {
    $sum += $group_member; 
  }
  my $mean = $sum / $sample_size;
  return $mean;
}

# computes the variance of a set of bins
sub variance {
  my $bins_ref  = shift; # the main bin hash
  my $group_ref = shift; # the group members, a subset of the bin hash
  my $sample_size = 0;
  my $sum = 0;
  my $mean = mean($bins_ref,$group_ref);
  my $diff = 0;
  foreach my $group_member (@{$group_ref}) {
    $diff = $group_member - $mean;
    $sum += $diff * $diff * $bins_ref->{$group_member}; # !! Assumes that all group member scores are >0
    $sample_size += $bins_ref->{$group_member};
  }
  my $variance = $sum / $sample_size;
  return $variance;
}

sub variance_a {
  my $group_ref = shift; # the group members, a subset of the bin hash
  my $sample_size = @{$group_ref};
  my $sum = 0;
  my $mean = mean_a($group_ref);
  my $diff = 0;
  foreach my $group_member (@{$group_ref}) {
    $diff = $group_member - $mean;
    $sum += $diff * $diff - $group_member;
  }
  my $variance = $sum / $sample_size;
  return $variance;
}

# open file - put in initial set of bins
open(F14,"<$F14");
while (<F14>) {
  if ( 2 == $. ) {
    # extact out number of elements and nodes
    $_ =~ m/(\d+)\s+(\d+)/;
    $NE = $1;
    $NN = $2;
    if(defined($DIVISOR)) {
      $MINCOUNT = $NN/($DIVISOR+1);
    }
  }
  # avoid line 1 /and/ line 2
  elsif ( 1 != $. ) {
    # extract out nodes and count in bin
    $_ =~ m/^\s+(\d+)\s+([+\-]*\d+\.\d+)\s+([+\-]*\d+\.\d+)/;
    my $node = $1;
    my $x    = $2;
    my $y    = $3;
    $bins{sprintf("%.2f",$y)}++;
  } 
  # jump out of loop when all nodes have been looked at
  last if ($NN+2 == $.);
}
close(F14);

my $sum = 0;
my $num_bins = 0;
foreach (keys(%bins)) {
 $sum += $bins{$_};
 ++$num_bins;
 #print "$_ .... $bins{$_} \n";
}
($sum == $NN) ? print STDERR "OK\n" : die "Not all nodes considered!";

# collapse %bins

my %collapsedbins = (); # maps which original bins in %bins are grouped together
my %newbins = ();
my $MAX = 0;
my $MIN = $NN+1;
my @I = ();
my $s = 0; my $n = 0; my $c = 0;
foreach $c (sort {$a <=> $b} keys(%bins)) {
  $s += $bins{$c};
  $collapsedbins{$c}=$n; # get last bin
  if ($s >= $MINCOUNT) {
    $newbins{$n} = $s;
    $MAX = $s if ($s > $MAX);
    $MIN = $s if ($s < $MIN);
    push(@I,$s);
    $s = 0; # reset $s
    ++$n;   # increment $n
  }
}
$collapsedbins{$c}=$n; # get last bin
$newbins{$n} = $s;
# stats
$MAX = $s if ($s > $MAX);
$MIN = $s if ($s < $MIN);
push(@I,$s);
my $mean = mean_a(\@I);
my $variance = variance_a(\@I);
# report
print STDERR "$n partitions created (min:$MIN, max:$MAX, mean:$mean, variance:$variance)\n";
  
$sum = 0;
$num_bins = 0;
foreach (sort {$a <=> $b} keys(%newbins)) {
 $sum += $newbins{$_};
 ++$num_bins;
 #print "$_ .... $newbins{$_} \n";
}
($sum == $NN) ? print STDERR "OK\n" : die "Not all nodes considered!";

#
# Where the partition happens after collapse...this is a naive approach that 
# simply dumps the number of beens created wrt min-number per bin - need to put in 
# kundu's algorith for real power and flexibility...this is just testing the
# creation of a partmesh.txt file
#

open(F14,"<$F14");
while (<F14>) {
  if ( 2 == $. ) {
    # extact out number of elements and nodes
    $_ =~ m/(\d+)\s+(\d+)/;
    $NE = $1;
    $NN = $2;
  }
  # avoid line 1 /and/ line 2
  elsif ( 1 != $. ) {
    # extract out nodes and count in bin
    $_ =~ m/^\s+(\d+)\s+([+\-]*\d+\.\d+)\s+([+\-]*\d+\.\d+)/;
    my $node = $1;
    my $x    = $2;
    my $y    = $3;
    my $b = sprintf("%.2f",$y);
    # find out what new bin $b is in, then output: <nodenum> <newbin>
    printf("\t%s\n",$collapsedbins{$b});  
  } 
  # jump out of loop when all nodes have been looked at
  last if ($NN+2 == $.);
}
close(F14);

1;
