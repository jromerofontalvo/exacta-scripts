set terminal postscript color eps enhanced "Palatino" 20
set output "overlap.eps"

set border 31 lw 2
set key spacing 2.0
set pointsize 1.5
set key top right

set yrange [0.2:1.1]
set xlabel "r/r_0" font "Palatino-bold, 20"
set ylabel "{/Symbol=26 \174}{/Symbol=26 \341}{/Symbol y}_{guess}{/Symbol=26 \174}{/Symbol y}_{exact}{/Symbol=26 \361}{/Symbol=26 \174}^2" font "Palatino-bold, 20"

plot "data_overlap" u 1:2 with points pt 7 t "", "data_overlap" u 1:2 smooth csplines lt 1 lc 1 lw 3 t "HF", "data_overlap_mps" u 1:2 with points pt 7 t "", "data_overlap_mps" u 1:2 smooth csplines lt 1 lc 3 lw 3 t "MPS_{M=2}"

