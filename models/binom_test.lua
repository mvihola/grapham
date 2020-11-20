p = 1/3
n = 100000
m = n*p

model = {
   x = {
      kind = "integer",
      init_val = 0,
      density = function(k)
                   return dbinom(k, n, p)
                end
   }
}

function functional()
   return {x-m}
end

para = {
   niter = 100e3,
   nburn = 10e3,
}
