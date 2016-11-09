# Scale fator test

run with
```
../../build/runmytrim < xe_on_uo2.json
../../build/runmytrim < xe_on_uo2_scale10.json
```

And plot using `gnuplot` with
```
pl 'xe_on_uo2_vac.dat' w l, 'xe_on_uo2_scale10_vac.dat' u ($1*10):($2/10) w l
```
