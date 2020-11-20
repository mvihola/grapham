sigma = 1e9;
-- the intensity parameter of the exponential distribution
lambda = 1;

model = {
  x = {
    parents = {"a","b","c"},
    dim = 2,
    density = function(x,a,b,c)
      return(dmvnorm(x, {a,c}, {{b,b/2},{b/2,b}}))
    end,
  },
  b = {
    parents = {},
    dim = 1,
    init_val = 1,
    density = function(b)
      return(dexp(b,lambda));
    end
  },
  a = {
    parents = {},
    dim = 1,
    density = function(a)
      return(dnorm(a,1,sigma))
    end
  },
  c = {
     parents = {},
     dim = 1,
     density = function(c)
       return(dnorm(c,1,sigma))
     end
  },
  z = {
     parents = {"a","b"},
     dim = 2,
     density = function(z,a,b)
       b2 = exp(-b)
       return(dmvnorm(z,{a,a},{{1,b2},{b2,1}}))
       --return(d.mvnorm_full(z,{a,a},{{1,b2},{b2,1}}))
     end
  },
}

data = {
  x = {1,20}
}

para = {
niter = 50e3, algorithm="am", blocks = {{"a","b","c","z"}}
}

functional = function() 
   return {a,b,c,x[1],x[2],z[1],z[2]}
end
