theta = 2; k = 3;
m = theta*k
m2 = (theta*k)^2+k*theta^2

model = {
  x = {
    init_val = 1,
    density = function(x)
      return dgamma(x, k, theta)
    end
  }
}

functional = function()
  return {x-m,x^2-m2}
end

para = {
  niter = 500e3, nburn = 10e3 
}
