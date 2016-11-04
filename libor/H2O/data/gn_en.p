set terminal postscript color eps enhanced "Palatino" 20
set output "energy.eps"

set border 31 lw 2
set key spacing 2.0
set pointsize 1.5
set key top right

set xlabel "r/r_0" font "Palatino-bold, 20"
set ylabel "Energy (a.u.)" font "Palatino-bold, 20"

plot "data_en" u 1:2 with points pt 7 t "", "data_en" u 1:2 smooth csplines lt 1 lc 1 lw 3 t ""

