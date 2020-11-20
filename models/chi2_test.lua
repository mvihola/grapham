k = 3;
m = k;
m2 = m^2+2*k;

model = {
   x = {
      init_val = 1,
      density = function(x)
                   return dchi2(x, k)
                end
   }
}

para = {
   niter = 50e3,
   nburn = 10e3,
}

function functional() 
   return {x-m, x^2-m2}
end
