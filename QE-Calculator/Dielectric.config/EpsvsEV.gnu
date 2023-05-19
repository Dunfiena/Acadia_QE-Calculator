set term png
set output "Epsilon vs Energy (eV)"
set xlabel 'Energy (eV)'
set ylabel 'ε'
set xtic nomirror
set ytic nomirror
set autoscale
unset key
set xr [*:*]
set yr [*:*]
set border 3
plot "epsi_Si.dat" using 1:2 w lines, \
 "epsr_Si.dat" using 1:2 w lines
