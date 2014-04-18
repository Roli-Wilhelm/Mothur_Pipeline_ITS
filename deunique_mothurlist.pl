#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $in1;
my $in2;
my $out;

my $result=GetOptions(
  "names=s"=>\$in1,
  "mothur=s"=>\$in2,
  "output=s"=>\$out);

if (!defined $in1 or !defined $in2 or !defined $out) {
  die("Usage: $0 -n namesfile -m mothurfile -o outputfile

  This program will take a name file in the following format:
    uniquename1     name1,name2,name3
    uniquename2     name4
    uniquename3     name5,name6
  and a mothur file:
    0.99        1234    uniquename1,uniquename2 uniquename3;
  and create the output file:
    0.99        1234    name1,name2,name3,name4 name5,name6;
");
}

open(IN1,"<$in1") or die("$in1: $!");
my %r;
while (<IN1>) {
  chomp;
  my ($aa,$bb)=split(/\t/,$_); # m/^(\S+)\t(\S.*)/;
  $r{$aa}=$bb;
}
close(IN1);

open(IN2,"<$in2") or die("$in2: $!");
open(OUT,">$out") or die("$out: $!");
while(<IN2>) {
  chomp;
  my ($l1,$l2,$l3)=(m/^(\S+)\t(\S+)\t(.*);?$/);
  my @s=split(/\t/,$l3);
  print OUT "$l1\t$l2";
  foreach my $s (@s) {
    print OUT "\t";
    my @o;
    foreach my $i (split(/,/,$s)) {
      if (!defined $r{$i}) {push @o,$i;next;}
      push @o,$r{$i};
    }
    print OUT join(",",@o);
  }
  print OUT "\n";
}
close(OUT);
close(IN2);
