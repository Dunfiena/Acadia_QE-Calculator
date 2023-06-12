set term png
set output "Total Energy vs Lattice Parameter"
set ylabel 'Total Energy (Ry)'
set xlabel 'Lattice Parameter (A)'
set xtic nomirror
set ytic nomirror
set autoscale
unset key
set xr [*:*]
set yr [*:*]
set border 3
plot "VolumeVsEnergy.dat" using 1:3 w lines
