set term png
set output "Total Energy vs KPoints"
set xlabel 'Kpoints'
set ylabel 'Total energy (Ry)'
set xtic nomirror
set ytic nomirror
set autoscale
unset key
set xr [*:*]
set yr [*:*]
set border 3
plot "Kpoint.dat" using 1:2 w lines
