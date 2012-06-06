#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;

#
# B. Estrade <estrabd@lsu.edu> 
# 17 June 2008
#

my $F14 = "./fort.14";
my $SUBDOMAINS = undef;
my $P = 5;
my $V = undef;
my $ABSOLUTEMIN = 1;

GetOptions("f=s" => \$F14, 
           "n=s" => \$SUBDOMAINS,
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
   -n             Number of subdomains
   -f             Specifies file (default is ./fort.14)
   -p             Defines precision in which to associate grid nodes to original bin
END
  exit;
}

my $NE = 0;
my $NN = 0;
# number of significant figures to use when counting nodes as being in the same longitude
my $precision = 5; 
my %bins = ();

# this "higher order" function reads over each node in the fort.14 and applies 
# $sub_ref (user provided subroutine reference) to the longitudinal position, $y
sub fold_f14 {
  my $sub_ref = shift;
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
      $sub_ref->($y);
      #$bins{sprintf("%.".$P."f",$y)}++;
    } 
    # jump out of loop when all nodes have been looked at
    last if ($NN+2 == $.);
  }
  close(F14);
}

fold_f14(sub {my $y = shift; $bins{sprintf("%.".$P."f",$y)}++});

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

my $MAX = 0;
my $MIN = $NN+1;

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

sub greedy_step {
  my $C = shift;    # hash ref of frequency keyed by bin; bin => freq
  my $M = shift;    # the number of intervals in which to partition the bins
  my $DIFF = shift; # used to determine weight between grouping and ideal
  my @orderedbins = sort{$a <=> $b} keys(%{$C});
  my $N = @orderedbins;
     --$N;
  my @node = ();
  my $sum = 0;
  for my $i (0 .. $N) {
    $sum += $C->{$orderedbins[$i]};
    if ($sum >= $DIFF) {
      push(@node,$orderedbins[$i]);
      $sum = 0;
    }
  }
  push(@node,$orderedbins[$N]);
  return @node;
}

my $IDEAL = ($NN / $SUBDOMAINS);
my @partition = greedy_step(\%bins,$SUBDOMAINS,$IDEAL);

sub getpartition {
  my $partition = shift;   # array ref containing @partition
  my $b = shift;           # node bin 
  my $num_partitions = @{$partition};
  my $count = 1;
  for my $i (0 .. $num_partitions-2) {
    if ($b >= $partition->[$i] && $b < $partition->[$i+1]) {
      # return if internal partition found
      return $count;
    }
    $count++;
  }
  # if gotten this far, return count, since it means it's in the last partition
  return $count;
}

my %partcounts = ();

# define sub ref to pass into fold_f14
my $sub_ref = sub {
  my $y = shift;
  my $b = sprintf("%.".$P."f",$y);
  # find out what new bin $b is in, then output: <nodenum> <newbin>
  my $p = getpartition(\@partition,$b);
  $partcounts{$partition[$p-1]}++;
  printf("\t%s\n",$p);  
};

fold_f14($sub_ref);

$sum = 0;
$MAX = 0;
$MIN = $NN+1;
my $c = 0;
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
