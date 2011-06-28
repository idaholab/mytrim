#!/usr/bin/perl

$min=0;
$max=2500;
#$step=0.2;
$step=1.0;

$dp = 0;
while(<STDIN>)
{
  $i = ($_-$min)/$step;
  if($i>=0) { $hist[$i]++; }
  $dp++;
}

for( $i=0; $i< ($max-$min)/$step; $i++ )
{
  $c = $hist[$i] / ($step*$dp);
  $x1 = $i*$step+$min;
  $x2 = $x1 + $step;
  print "$x1\t$c\n";
  print "$x2\t$c\n";
}

