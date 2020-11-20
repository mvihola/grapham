const = {
   mu = 0,
   v = 1
}

mm = exp(const.mu+const.v/2)
mm2 = (exp(const.v)-1)*exp(2*const.mu+const.v) + mm^2
lmm = const.mu
lmm2 = const.v + const.mu^2

model = {
   x = {
      init_val = 1,
      parents = {"mu","v"},
      density = "dlognorm"
   }
}

function functional() 
   local lnx = log(x)
   return {x-mm,x*x-mm2,lnx-lmm,lnx*lnx-lmm2}
end

para = {
   nburn = 10e3,
   niter = 100e3,
}
