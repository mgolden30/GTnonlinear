#!/usr/bin/gnuplot

set terminal png enhanced background rgb "black"

set output "orbit.png"

set xrange [-3: 3]
set yrange [-3: 3]

set title "PO3" textcolor rgb "white"

set size square

set nokey

set xlabel "x"
set ylabel "y"

plot "orbit.dat" using 1:2 with lines lc rgb "cyan" lw 3, \
     "orbit.dat" using 3:4 with lines lc rgb "blue" lw 3, \
     "orbit.dat" using (-$1-$3):(-$2-$4) with lines lc rgb "green" lw 3, \
