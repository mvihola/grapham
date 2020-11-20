mu = 3;
beta = 4; 
-- Euler-Mascheroni constant
gamma = 0.5772156649015328606
m = mu+beta*gamma
m2 = m^2+pi^2/6*beta^2

model = {
  x = {
    init_val = 1,
    density = function(x)
      return dgumbel(x, mu, beta)
    end
  }
}

functional = function()
  return {x-m,x^2-m2}
end

para = {
  niter = 100e3, nburn = 10e3 
}
