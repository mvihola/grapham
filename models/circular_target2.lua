model = {
   x = {
      dim = 2,
      density = function(x)
        nsqx = (15-sqrt(x[1]^2+x[2]^2))^2
        return -nsqx + log((1+cos(5*x[1]))*(1+cos(5*x[2])))
      end,
   }
}

para = {
   niter = 50e3,
   nburn = 10e3,
   algorithm = "asm",
   blocking = "full"
--   outfile = "circular.csv",
}

functional = function() return x end
