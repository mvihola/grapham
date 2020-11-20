requires_numlua()

-- the parameters
N = 8;
a = 1; b = 1;
-- some shortcuts and precalculated stuff
aa = a*a; ab = a*b;
c = matrix(N):fill(1)
c[1] = 100; sqrtc1 = sqrt(c[1])
zer = matrix(N):fill(0)
tmp = zer:copy()
pp = {3.4895,5.0706,7.3441,10.2189,13.3616 }
true_pp = {.1, .25, .5, .75, .9};
Npp = #pp
res = matrix(Npp);
model = {
  x = {
    dim = N,
    density = function(x)
        tmp[1] = x[1]/a;
        tmp[2] = x[2]*a + ab*(x[1]*x[1] + aa);
        for k=3,N do tmp[k] = x[k] end
      return(dmvnorm_diag(tmp, zer, c))
   end,
  }
}

para = {
   niter = 50e3, 
   nburn = 10e3,
}

function functional()
   -- "banana analysis", i.e. 
  tmp[1] = x[1]/(a*sqrtc1);
  tmp[2] = x[2]*a + ab*(x[1]*x[1] + aa);
  for k=3,N do tmp[k] = x[k] end
  z = matrix.inner(tmp,tmp);
  for k=1,Npp do res[k] = 1-true_pp[k] end
  for k=1,Npp do 
     if z<=pp[k] then break else res[k] = res[k]-1 end
  end
  return res
end

