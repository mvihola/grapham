p = 1/3
r = 0.2

m = r*(1-p)/p
m2 = m^2 + r*(1-p)/p^2

model = {
   x = {
      kind = "integer",
      density = function(k)
                   return dnbinom(k, r, p)
                end
   }
}

function functional()
   return {x-m, x*x-m2}
end

para = {
   niter = 100e3,
   nburn = 10e3,
}
