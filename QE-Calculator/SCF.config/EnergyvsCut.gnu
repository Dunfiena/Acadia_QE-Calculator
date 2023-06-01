set term png
set output "Total Ry vs Cutoff Ry"
set xlabel 'Energy Cutoff (Ry)'
set ylabel 'Total Energy (Ry)'
set xtic nomirror
set ytic nomirror
set autoscale
unset key
set xr [*:*]
set yr [*:*]
set border 3
plot "energy.dat" using 1:2 w lines
