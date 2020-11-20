-- A hidden Markov model; the Markov chain is a
-- one-dimensional "almost constant velocity model".
-- The "process noise" variable rho has an
-- improper flat prior in (0,\infty)
-- The initial position x0 has also
-- an improper flat prior.

m0 = {0,0}
c0 = 100; 
cyk = 0.05;

dt = 1
l11 = pow(dt, 3/2)/sqrt(3)
l12 = sqrt(3)/2*sqrt(dt)
l22 = sqrt(dt-3/4*dt)

xpred = {0,0}

ydata = {-3, -2, 0, 1, 1, 0, -2, -5};

model = {
   rho_var = {
      init_val = 1,
      density = function(rho)
                   if rho<=0 then 
                      return NINF
                   else
                      return 1
                   end
                end
   },
   x0 = {
      dim = 2,
      density = function(x)
                   return 1
                end
   },
   x = {
      dim = 2,
      parents = {"x-","rho_var"},
      density = function(x,xp,rho)
                   xpred[1] = xp[1]+xp[2]*dt
                   xpred[2] = xp[2]
                   return dmvnorm_chol(x, xpred, {{rho*l11},
                                       {rho*l12,rho*l22}})
                end
   },
   y = {
      parents = {"x"},
      density = function(y,x)
                   return dnorm(y, x[1], cyk)
                end
   }
}

repeat_block({"y","x"}, ydata)

function functional() 
   return {rho_var,x0[1],x1[1],x2[1],x3[1],x4[1],x5[1],x6[1],x7[1],x8[1]}
end

para = {
   nburn = 10e3, 
   niter = 30e3, 
   nthin = 10,
   algorithm = "aswam",
   blocking = "node"
}
