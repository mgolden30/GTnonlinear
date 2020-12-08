#!/usr/bin/gnuplot

set terminal png enhanced

do for [ii=1:120]{

  file = sprintf("shapesphere/%03d.png", ii)
  set output file

  set title "PO on shape-sphere"

  set size square

  set nokey

  set xrange [-1:1]
  set yrange [-1:1]
  set zrange [-1:1]

  theta = 6.28318530718/120*ii

  #I will define the Jacobi coordinates
  
  rho_x(x1,y1,x2,y2,x3,y3) = x1 - x2
  rho_y(x1,y1,x2,y2,x3,y3) = y1 - y2
  lambda_x(x1,y1,x2,y2,x3,y3) = x1 + x2 - 2*x3
  lambda_y(x1,y1,x2,y2,x3,y3) = y1 + y2 - 2*y3

  nx(x1,y1,x2,y2,x3,y3) = 2*( (x1 - x2)*(x1 + x2 - 2*x3) + (y1 - y2)*(y1 + y2 - 2*y3) )/\
     (rho_x(x1,y1,x2,y2,x3,y3)**2 + \
      rho_y(x1,y1,x2,y2,x3,y3)**2 + \
      lambda_x(x1,y1,x2,y2,x3,y3)**2 + \
      lambda_y(x1,y1,x2,y2,x3,y3)**2)




  ny(x1,y1,x2,y2,x3,y3) = 2*( (x1 - x2)*(y1 + y2 - 2*y3) - (y1 - y2)*(x1 + x2 - 2*x3) )/\
     (rho_x(x1,y1,x2,y2,x3,y3)**2 + \
      rho_y(x1,y1,x2,y2,x3,y3)**2 + \
      lambda_x(x1,y1,x2,y2,x3,y3)**2 + \
      lambda_y(x1,y1,x2,y2,x3,y3)**2)




  nz(x1,y1,x2,y2,x3,y3) = ( - rho_x(x1,y1,x2,y2,x3,y3)**2 - \
      rho_y(x1,y1,x2,y2,x3,y3)**2 + \
      lambda_x(x1,y1,x2,y2,x3,y3)**2 + \
      lambda_y(x1,y1,x2,y2,x3,y3)**2) / \
     (rho_x(x1,y1,x2,y2,x3,y3)**2 + \
      rho_y(x1,y1,x2,y2,x3,y3)**2 + \
      lambda_x(x1,y1,x2,y2,x3,y3)**2 + \
      lambda_y(x1,y1,x2,y2,x3,y3)**2)

  set parametric
  set urange [0:2*pi]
  set vrange [-pi/2:pi/2]
  # Parametric functions for the sphere
  fx(v,u) = cos(v)*cos(u)
  fy(v,u) = cos(v)*sin(u)
  fz(v)   = sin(v)

  splot "orbit.dat" using (cos(theta)*nx($2,$3,$4,$5,(-$2-$4),(-$3-$5)) + sin(theta)*ny($2,$3,$4,$5,(-$2-$4),(-$3-$5)) ):( cos(theta)*ny($2,$3,$4,$5,(-$2-$4),(-$3-$5)) - sin(theta)*nx($2,$3,$4,$5,(-$2-$4),(-$3-$5)) ):(nz($2,$3,$4,$5,(-$2-$4),(-$3-$5))) with lines lc rgb "blue" lw 3, \
      "collisions.dat" using ($1*cos(theta) + sin(theta)*$2):($2*cos(theta) - $1*sin(theta)):3 w p lc rgb 'red' pt 9,\
      fx(v,u),fy(v,u),fz(v)
}
