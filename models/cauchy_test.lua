x0 = 0; gamma = 1

model = {
   x = {
      density = function(x)
                   return dlaplace(x, x0, gamma)
                end
   }
}

para = {
   niter = 500e3,
   nburn = 10e3,
   -- nthin = 5,
   -- outfile = "cauchy.bin"
}

_a = 0; _b = 1;
function functional()
   _a = min(_a,x);
   _b = max(_b,x);
   return {abs(x)}
end

para.close_hook = function()
                     print("min: ".. _a .. ", max: " .._b)
                  end                    
                  