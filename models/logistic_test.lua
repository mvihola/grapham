mu = 2; s = 2;
m2 = mu*mu+pi^2/3*s^2

model = {
  x = {
    init_val = 1,
    density = function(x)
      return dlogistic(x, mu, s)
    end
  }
}

functional = function()
  return {x-mu, x*x-m2}
end

para = {
   niter = 100e3, nburn = 10e3 , 
   algorithm="am"
}
