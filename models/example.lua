model = {
   x1 = {
      density = function(x) 
                   return dnorm(x,0,1)
                end
   },
   x2 = {
      parents = {"x1"},
      density = function(x,x1)
                   return dexp(x,1+x1^2)
                end
   },
   x3 = {
      parents = {"x1"},
      density = function(x,x1)
                   return dnorm(x,x1,0.1)
                end
   }
}

function functional()                
   return {x1,x2,x3}
end

para = {
   niter = 1e5
}
