set term png
set output "Runtime vs Cutoff Ry"
set xlabel 'Energy Cutoff (Ry)'
set ylabel 'Runtime (s)'
set xtic nomirror
set ytic nomirror
set autoscale
unset key
set xr [*:*]
set yr [*:*]
set border 3
plot "energy.dat" using 1:3 w lines
