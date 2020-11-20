mu = {1,2};
v = {{.2,0},{.1,.8}}
p = 2.1
m21 = mu[1]*mu[2]+p/(p-2)*v[2][1]
m11 = mu[1]*mu[1]+p/(p-2)*v[1][1]
m22 = mu[2]*mu[2]+p/(p-2)*v[2][2]

model = {
   x = {
      dim = 2,
      density = function(x)
                   return dmvstudent(x, mu, v, p)
                end
   }
}

function functional()
   return {x[1]-mu[1],x[2]-mu[2], x[1]*x[1]-m11, x[2]*x[2]-m22, x[2]*x[1]-m21}
end

para = {
   niter = 100e3,
   nburn = 10e3
}

                