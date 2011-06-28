#!/usr/bin/perl

$min=0;
$max=2500;
#$step=0.2;
$step=1.0;

# density of TiB2
$rho = 4.52;

# weight per formula unit
$utib2 = 2*11.0 + 48.0;

# radius in nm
$r = 1.0;

# volume in cm^3 * 1e21
$v = 4.0/3.0 * 3.141 * $r**3;

# formula units per bubble
$a = $rho*$v/$utib2 * 6.022e2; 
# 1e2 = 1e23*1e-21

# clusters per sample
$N = 75;

$dp = 0;
while(<STDIN>)
{
  $i = ($_-$min)/$step;
  if($i>=0) { $hist[$i]++; }
  $dp++;
}

for( $i=0; $i< ($max-$min)/$step; $i++ )
{
  $c = 3.0*$a*$N * $hist[$i] / ($step*$dp);
  $x1 = $i*$step+$min;
  $x2 = $x1 + $step;
  print "$x1\t$c\n";
  print "$x2\t$c\n";
}



print 3*$a*$N;
