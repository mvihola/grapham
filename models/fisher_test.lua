d1=1; d2=7
mean = d2/(d2-2)
mom2 = 2*d2^2*(d1+d2-2)/(d1*(d2-2)^2*(d2-4)) + mean^2

model = {
  x = {
    init_val = 1,
    density = function(x)
      return dfisher(x, d1, d2)
    end
  }
}

functional = function()
  return {x-mean, x^2-mom2}
end

para = {
  niter = 200e3, nburn = 10e3 
}
