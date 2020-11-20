y = {0.395, 0.375, 0.355, 0.334, 0.313, 0.313, 0.291,
0.269, 0.247, 0.247, 0.224, 0.224, 0.224, 0.224, 0.224, 0.200,
0.175, 0.148}
v = 0.00434
a1=-1; a2=-1
b1=2;  b2=2
mu0=0; s0=1

ndata = #y

t0 = {};
for k = 1,ndata+2 do 
   t0[k] = 1 
end

-- as a single functional
model = {
   t = {
      dim = ndata+2,
      init_val = t0,
      density = function(t)
                   local mu = t[ndata+1];
                   local a = t[ndata+2];
                   if (a<0) then return -huge end
                   
                   local result = -(mu-mu0)*(mu-mu0) / 2.0 / s0 / s0 -
                   b1/a - (a1+1) * log(a) - b2/v - (a2+1) * log(v);
                   for kk=1,ndata do
                      result = result 
                      - 0.5 * log(a) - (t[kk]-mu)*(t[kk]-mu) / 2.0 / a 
                      - 0.5 * log(v) - (y[kk]-t[kk])*(y[kk]-t[kk]) / 2.0 / v;
                   end
                   return result
                end
   }
}

function functional()
   return {t[1],t[2],t[ndata+1],t[ndata+2]}
   --return t
end

para = {   
   niter = 30e3,
   nburn = 10e3,
   algorithm = "aswam",
   blocking = "sc",
   scaling_adapt = function(sc, alpha, dim, k)
                      local delta
                      if alpha>0.234 then
                         delta = min(0.01, 1/sqrt(k+1))
                      else
                         delta = -min(0.01, 1/sqrt(k+1))
                      end
                      return sc*exp(delta)
                   end
}
