#!/usr/bin/gnuplot

set terminal png enhanced

set output "orbit.png"

#set xrange [-2.5: 2.5]
#set yrange [-2.5: 2.5]

set size square

set nokey

set xlabel "x"
set ylabel "y"

plot "orbit.dat" using 1:2 with lines, \
     "orbit.dat" using 3:4 with lines, \
     "orbit.dat" using (-$1-$3):(-$2-$4) with lines, \
      
