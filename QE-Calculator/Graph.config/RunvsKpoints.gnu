set term png
set output "Runtime vs KPoints"
set xlabel 'Kpoint'
set ylabel 'Runtime (s)'
set xtic nomirror
set ytic nomirror
set autoscale
unset key
set xr [*:*]
set yr [*:*]
set border 3
plot "Kpoint.dat" using 1:3 w lines
