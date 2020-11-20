requires_numlua()

N = 3;
v = matrix.eye(N);
v[2][1]=0.1
n = 10;
x0 = matrix.eye(N);
--x0 = {{1,0,0},{0,1,0},{0,0,1}}
m = n*v
m11 = m[1][1]
m21 = m[2][1]
m22 = m[2][2]
m31 = m[3][1]
m32 = m[3][2]
m33 = m[3][3]

model = {
   x = {
      dim = {N,N},
      type = "matrix",
      kind = "symmetric",
      init_val = x0,
      density = function(x)
                   return dwishart(x, v, n)
                end
   },
}

function functional()
   return {x[1][1]-m11,x[2][1]-m21,x[2][2]-m22,
      x[3][1]-m31,x[3][2]-m32,x[3][3]-m33}
end

para = {
   niter = 100e3,
   nburn = 10e3,
   algorithm = "am",
   --outvars = {"x[3][3]"},
   --outfile = "wishart.bin"
}
