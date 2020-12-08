#!/usr/bin/gnuplot

set terminal png enhanced

set output "orbit.png"

set xrange [-3: 3]
set yrange [-3: 3]

set title "PO3" textcolor rgb "white"

set size square

set nokey

set xlabel "x"
set ylabel "y"

plot "orbit.dat" using 2:3 with lines lc rgb "cyan" lw 3, \
     "orbit.dat" using 4:5 with lines lc rgb "blue" lw 3, \
     "orbit.dat" using (-$2-$4):(-$3-$5) with lines lc rgb "green" lw 3, \
