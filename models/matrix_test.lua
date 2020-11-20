requires_numlua()

N = 3;
x0 = matrix(N,N):fill(1);

model = {
   x = {
      dim = {N,N},
      type = "matrix",
      init_val = x0,
      density = function(x)
                   return -.5*matrix.norm(x)^2
                end
   },
   y = {
      init_val = x0,
      dim = {N,N},
      type = "matrix",
      parents = {"x"},
      density = function(y,x)
                   return -.5*matrix.norm(y-x)^2
                end
   }
}

function functional()
   return {x[1][1],x[1][1]^2-1,y[1][1],y[1][1]^2-2}
end

para = {
   niter = 20e3,
   nburn = 10e3,
   init = "freeze",
   algorithm = "am",
   blocking = "full"
}
