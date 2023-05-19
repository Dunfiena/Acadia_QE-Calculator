set term png
set output "Density of state of Si Crystal"
set xlabel 'Energy (meV)'
set ylabel 'Density of States'
set xtic nomirror
set ytic nomirror
set autoscale
unset key
set xr [*:*]
set yr [*:*]
set border 3
plot "Si.phdos" using 1:2 w lines
