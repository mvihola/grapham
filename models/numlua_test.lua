--require("numlua.matrix")

N = 10; 
x0 = {}; y0 = {};
for k = 1,N do
   x0[k] = 0; y0[k] = 1;
end
if matrix ~= nil then
   print "Using Numlua matrix class"
   x0 = matrix(N):fill(0);
   y0 = matrix(N):fill(1);
else
   print "Not using Numlua matrix class"
end

model = {
   x = {
      init_val = x0,
      dim = N,
      type = "vector",
      density = function(x)
                   local sum = 0
                   for k = 1,N do
                      sum = sum - .5*x[k]^2;
                   end
                   return sum
                end
   },
   y = {
      init_val = y0,
      dim = N,
      type = "custom",
      parents = {"x"},
      density = function(y,x)
                   local sum = 0
                   for k = 1,N do
                      sum = sum - .5*(y[k]-x[k])^2
                   end
                   return sum
                end
   }
}

if matrix ~= nil then
   model.x.type = "matrix"
   model.x.density = function(x)
                   return -.5*matrix.norm(x)^2
                     end
   model.y.type = "matrix"
   model.y.density = function(y,x)
                   return -.5*matrix.norm(y-x)^2
                end
end

function functional()
   return {x[1],x[1]^2-1,y[1],y[1]^2-2}
end

para = {
   niter = 100e3,
   nburn = 10e3,
   init = "freeze",
   algorithm = "am",
   blocking = "full"
}
