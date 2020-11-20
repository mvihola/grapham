const = {
   m = 2,
   v = 3
}

-- These are half-normals and a symmetric truncation
model = {
   x = {
      parents = {"m","v"},
      density = "dnorm",
      init_val = const.m,
      limits = {const.m,nil}
   },
   y = {
      parents = {"m","v"},
      density = "dnorm",
      init_val = const.m,
      limits = {nil, const.m}
   },
   z = {
      parents = {"m","v"},
      density = "dnorm",
      init_val = const.m,
      limits = {const.m-sqrt(const.v), const.m+sqrt(const.v)}
   }
}

para = {
   niter = 100e3,
   nburn = 10e3,
   algorithm = "asm",
}

functional = function()
  x0 = x-const.m-sqrt(const.v*2/pi)
  y0 = y-const.m+sqrt(const.v*2/pi)  
  return 
      {x0, x0^2 - const.v*(1-2/pi), 
       y0, y0^2 - const.v*(1-2/pi), z-const.m}
end
