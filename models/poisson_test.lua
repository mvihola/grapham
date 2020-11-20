lambda = pi
m2 = lambda^2+lambda

model = {
   x = {
      kind = "integer",
      init_val = 1,
      density = function(k)
                   return dpoisson(k, lambda)
                end
   }
}

function functional()
   return {x-lambda, x*x-m2}
end

para = {
   niter = 100e3,
   nburn = 10e3,
}
