### 2.1 Graph Configuration
set term x11 persist size 900,600
set xtic nomirror
set ytic nomirror
set border 3
unset key
set size 1.0,1.0
set origin 0,0

### 2.2 Multiplot configuration ###
set multiplot layout 2,3 columnsfirst

# --- Graph 1,1 --- #
set size 0.4,0.35
set origin 0.0,0.6
set xlabel 'Energy Total (Ry)'
set ylabel 'Energy Cutoff (Ry)'
set xr [*:*]
set yr [*:*]
plot "./Graphs/energy.dat" using 1:2 w lines

# --- Graph 1,2 --- #
set size 0.4,0.35
set origin 0.5,0.6
set xlabel 'Energy Total (Ry)'
set ylabel 'Runtime (s)'
set xr [*:*]
set yr [*:*]
plot "./Graphs/energy.dat" using 1:3 w lines

# --- Graph 2,1 --- #
set size 0.4,0.35
set origin 0,0.3
set xlabel 'KPoints'
set ylabel 'Energy Cutoff (Ry)'
set xr [*:*]
set yr [*:*]
plot "./Graphs/Kpoint.dat" using 1:2 w lines

# --- Graph 2,2 --- #
set size 0.4,0.35
set origin 0.5,0.3
set xlabel 'Kpoints'
set ylabel 'Runtime (s)'
set xr [*:*]
set yr [*:*]
plot "./Graphs/Kpoint.dat" using 1:3 w lines

# --- Graph 3,1 --- #
set size 0.4,0.35
set origin 0.0,0.0
set ylabel 'Total Energy (Ry)'
set xlabel 'Lattice Parameter (A)'
set xr [*:*]
set yr [*:*]
plot "./Graphs/Datafiles/celldmtmp.dat" using 1:2 w lines

# --- Graph 3,2 --- #
set size 0.4,0.35
set origin 0.5,0.0
set ylabel 'Runtime (s)'
set xlabel 'Lattice Parameter (A)'
set xr [*:*]
set yr [*:*]
plot "./Graphs/Datafiles/celldmtmp.dat" using 1:3 w lines

unset multiplot

### 2.3 Live update configuration
pause 1
reread
