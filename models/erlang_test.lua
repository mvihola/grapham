k = 2;
lambda = 2; 
m = k/lambda
m2 = k^2/lambda^2+k/lambda^2

model = {
  x = {
    init_val = 1,
    density = function(x)
      return derlang(x, k, lambda)
    end
  }
}

functional = function()
  return {x-m,x^2-m2}
end

para = {
  niter = 100e3, nburn = 10e3 
}
