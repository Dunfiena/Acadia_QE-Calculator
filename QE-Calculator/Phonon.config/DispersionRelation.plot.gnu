set term png
unset key
set xlabel "Kpoint path [2π/a]"
set xtics ("L" 0.0, "G" 0.866, "X" 1.866)
set ylabel "Energy (meV)"
set output "Dispertion relation graph"
plot [0:1.866] for [i=2:7] "Si.freqNU.gp" using 1:(column(i)/8.0655) w l lc 3 lw 2
