#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;

#
# B. Estrade <estrabd@lsu.edu> 
# 17 June 2008
#
# This is a reference implementation of the kstep row-based partitioning
# of ADCIRC grids, which to minimize the number of communicating 
# interfaces between domains.  The balancing of grid nodes is accomplished
# using a k-step shortest path formulation based on the k-step shortest path
# algorithm presented in [Kundu98] below.
#
#@Article{Kundu98,
#  author =      "S. Kundu",
#  title =       "A Solution to Histogram-Equalization and Other Related Problems by Shortest-Path Methods",
#  journal =     "Pattern Recognition",
#  volume =      "31",
#  year =        "1998",
#  number =      "3",
#  month =       mar,
#  pages =       "231--234",
#  URL =         "http://www.sciencedirect.com/science/article/B6V14-3SYS091-R/2/8afabb8cfeb2a538c56f7b071e01f977",
#  bibsource =   "http://www.visionbib.com/bibliography/image-proc216.html#TT15790",
#}

my $F14 = "./fort.14";
my $SUBDOMAINS = undef;
my $P = 5;
my $V = undef;
my $MINCOUNT = 10; # adjusted with -m option
my $ABSOLUTEMIN = 1;

GetOptions("f=s" => \$F14, 
           "d=s" => \$SUBDOMAINS,
           "m=s" => \$MINCOUNT,
           "p=s" => \$P,);
           #"h"   => &show_help());

sub show_help {
  print<<END;

Reference implementation of row based grid partitioning for ADCIRC grid files. 

Input

  ADCIRC grid file.

Output

  ADCIRC partmesh.txt, used by adcprep >= 47.x to define subdomains.

Usage

  adcpart-smart.pl -d <numdomains> [-f ./fort.14] [-p <precision>] [-m <min_nodes_per_bin>] > partmesh.txt 

 ... then run adcprep, skipping the option (1) that creates partmesh.txt.

Options        Description
   -d             Number of subdomains
   -f             Specifies file (default is ./fort.14)
   -m             Defines minimum number of nodes in the new bins when reduced 
   -p             Defines precision in which to associate grid nodes to original bin
END
  exit;
}

my $NE = 0;
my $NN = 0;
# number of significant figures to use when counting nodes as being in the same longitude
my $precision = 5; 
my %bins = ();

# open file - put in initial set of bins
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
    $_ =~ m/^\s*(\d+)\s+([+\-]*\d+\.\d+)\s+([+\-]*\d+\.\d+)/;
    my $node = $1;
    my $x    = $2;
    my $y    = $3;
    $bins{sprintf("%.".$P."f",$y)}++;
  } 
  # jump out of loop when all nodes have been looked at
  last if ($NN+2 == $.);
}
close(F14);

my $sum = 0;
my $num_bins = 0;
my @orderedbins = sort{$a <=> $b} keys(%bins);
foreach (@orderedbins) {
 $sum += $bins{$_};
 ++$num_bins;
 #print "$_ .... $bins{$_} \n";
}
($sum == $NN) ? print STDERR "OK\n" : die "Not all nodes considered!";
my $bc = @orderedbins;
print STDERR "$bc bins found\n";

# collapse %bins
my %newbins = ();
my $MAX = 0;
my $MIN = $NN+1;
my @I = ();
my $s = 0; my $n = 0; my $c = 0;
foreach $c (@orderedbins) {
  $s += $bins{$c};
  if ($s >= $MINCOUNT) {
    $newbins{$orderedbins[$n]} = $s;
    $MAX = $s if ($s > $MAX);
    $MIN = $s if ($s < $MIN);
    push(@I,$s);
    $s = 0; # reset $s
  }
  ++$n;   # increment $n
}
$newbins{$orderedbins[$n-1]} = $s;
%bins = %newbins;
# stats
$MAX = $s if ($s > $MAX);
$MIN = $s if ($s < $MIN);
push(@I,$s);
# report
@orderedbins = sort{$a <=> $b} keys(%bins);
$bc = @orderedbins;
print STDERR "reduced to $bc partitions (min:$MIN, max:$MAX)\n";
%newbins = ();

# creates the arithmetic mean of a given population
sub mean {
  my @group_ref = @_; # the group members, a subset of the bin hash
  my $sample_size = 0;
  my $sum = 0;
  foreach my $group_member (@group_ref) {
    $sum += $group_member * $group_member; # !! Assumes that all group member scores are >0
    $sample_size += $group_member;
  }
  my $mean = $sum / $sample_size;
  return $mean;
}

# computes the variance of a set of bins
sub variance {
  my @group_ref = @_; # the group members, a subset of the bin hash
  my $sample_size = 0;
  my $sum = 0;
  my $mean = mean(@group_ref);
  my $diff = 0;
  foreach my $group_member (@group_ref) {
    $diff = $group_member - $mean;
    $sum += $diff * $diff; # !! Assumes that all group member scores are >0
    $sample_size += $group_member;
  }
  my $variance = $sum / $sample_size;
  return $variance;
}

sub kstep {
  my $C = shift;    # hash ref of frequency keyed by bin; bin => freq
  my $M = shift;    # the number of intervals in which to partition the bins
  my $DIFF = shift; # used to determine weight between grouping and ideal
  my @orderedbins = sort{$a <=> $b} keys(%{$C});
  my $N = @orderedbins;
     --$N;
  my %s = ();
  my %d = ();
  my %minLength = ();
  my %node = ();
  for my $i (0 .. $N) {
    $s{$i}->{$i} = 0;
    for my $j ($i+1 .. $N) {
      $s{$i}->{$j} = $s{$i}->{$j-1} + $C->{$orderedbins[$j]};
      $d{$i}->{$j} = abs($s{$i}->{$j} - $DIFF) ;
    }
  }
  for my $j (1 .. $N) {
    $minLength{1}->{$j} = $d{0}->{$j};
    $node{1}->{$j} = 0;
  }
  for my $k (2 .. $M) {
    for my $j ($k .. $N) {
      $minLength{$k}->{$j} = $minLength{$k-1}->{$k-1} + $d{$k-1}->{$j};
      $node{$k}->{$j} = $k-1;
      for my $i ($k .. $j-1) {
        if ($minLength{$k-1}->{$i} + $d{$i}->{$j} < $minLength{$k}->{$j}) {
          $minLength{$k}->{$j} = $minLength{$k-1}->{$i} + $d{$i}->{$j};
          $node{$k}->{$j} = $i;
        }
      }
    }
  }
  my @ret = ();
  my $lastnode = $N;
  for (my $k=1;$k <= $M; $k++) {
    push(@ret,$node{$k}->{$lastnode});
  }
  return @ret;
}

my $IDEAL = $NN / $SUBDOMAINS;
print STDERR "$NN / $SUBDOMAINS = $IDEAL\n";
my @partition = kstep(\%bins,$SUBDOMAINS,$IDEAL);

sub getpartition {
  my $partition = shift;   # array ref containing partition
  my $orderedbins = shift; # ordered array ref of actual bins
  my $b = shift;           # node bin 
  my $num_partitions = @{$partition};
  my $count = 1;
  for my $i (0 .. $num_partitions-2) {
    #printf( "%s >= %s && %s < %s\n",$b,$orderedbins[$partition->[$i]],$b,$orderedbins[$partition->[$i+1]]);
    if ($b >= $orderedbins[$partition->[$i]] && $b < $orderedbins[$partition->[$i+1]]) {
      # return if internal partition found
      return $count;
    }
    $count++;
  }
  # if gotten this far, return count, since it means it's in the last partition
  return $count;
}

my %partcounts = ();
# output partmesh 
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
    $_ =~ m/^\s*(\d+)\s+([+\-]*\d+\.\d+)\s+([+\-]*\d+\.\d+)/;
    my $node = $1;
    my $x    = $2;
    my $y    = $3;
    my $b = sprintf("%.".$P."f",$y);
    # find out what new bin $b is in, then output: <nodenum> <newbin>
    my $p = getpartition(\@partition,\@orderedbins,$b);
    $partcounts{$orderedbins[$partition[$p-1]]}++;
    printf("\t%s\n",$p);  
  } 
  # jump out of loop when all nodes have been looked at
  last if ($NN+2 == $.);
}
close(F14);

$sum = 0;
$c = 0;
$MAX = 0;
$MIN = $NN+1;
my @counts = ();
foreach (sort{$a<=>$b} keys(%partcounts)) {
  print STDERR "(".++$c.") $_ - count = $partcounts{$_}\n";
  $sum += $partcounts{$_};
  push(@counts,$partcounts{$_});
  $MAX = $partcounts{$_} if ($partcounts{$_} > $MAX);
  $MIN = $partcounts{$_} if ($partcounts{$_} < $MIN);
}
my $variance = variance(@counts);
my $mean = mean(@counts);
my $num = @orderedbins;
print STDERR "stat ($c partitions from $num bins); max: $MAX, min: $MIN, mean: $mean, variance: $variance\n";
#print STDERR "$c & $num & $MAX & $MIN & $mean & $variance\n";
($sum == $NN) ? print STDERR "OK\n" : die "Not all nodes considered!";

1;
