#This is a script to play with parareal and analyze convergence properties
#

def parareal
  n  = 10
  dt = 1.0/n
  points = [1.0]
  n.times{
    state = points[-1]
    points.push( state - dt*state*state )
  }

  puts points
  

end

parareal
