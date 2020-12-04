#!/usr/bin/gnuplot

set terminal png enhanced background rgb "black"

num_frames = 240
do for [ii=1:num_frames] {
  name = sprintf("frames/%03d_frame.png", ii)
  set output name

  set xrange [-5: 5]
  set yrange [-5: 5]

  set title "PO1" textcolor rgb "white"

  set size square

  set nokey

  set xlabel "x"
  set ylabel "y"

  num_lines = 23767
  jj = num_lines / num_frames * ii

  plot \
       "orbit.dat" using 2:3 every ::1::jj with lines lc rgb "cyan" lw 2, \
       "orbit.dat" using 4:5 every ::1::jj with lines lc rgb "blue" lw 2, \
       "orbit.dat" using (-$2-$4):(-$3-$5) every ::1::jj with lines lc rgb "green" lw 2, \
       "orbit.dat" using 2:3 every ::jj::jj with points pt 7 lc rgb "white", \
       "orbit.dat" using 4:5 every ::jj::jj with points pt 7 lc rgb "white", \
       "orbit.dat" using (-$2-$4):(-$3-$5) every ::jj::jj with points pt 7 lc rgb "white"
}
