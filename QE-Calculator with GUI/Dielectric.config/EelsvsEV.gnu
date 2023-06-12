set term png
set output "Eels vs Energy (eV)"
set xlabel 'Energy (eV)'
set ylabel 'EELS Intensity (arb. units)'
set xtic nomirror
set ytic nomirror
set autoscale
unset key
set xr [*:*]
set yr [*:*]
set border 3
plot "eels_Si.dat" using 1:2 w lines
