s = 0.5
m = s*sqrt(pi/2)
m2 = m*m+(4-pi)/2*s*s

model = {
  x = {
    init_val = 1,
    density = function(x)
      return drayleigh(x, s)
    end
  }
}

functional = function()
  return {x-m, x*x-m2}
end

para = {
  niter = 100e3, nburn = 10e3,
}
