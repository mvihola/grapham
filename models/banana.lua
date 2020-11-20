-- the para
a = 1; b = 1;
aa = a*a; ab = a*b;

-- The covariance matrix
c = {{1, .9},{.9, 1}}
zer = {0,0}
tmp = {0,0}
pp = {0.21072103131565,
   0.57536414490356,
   1.38629436111989,
   2.77258872223979,
   4.60517018598809}
Npp = #pp
true_pp = {.1, .25, .5, .75, .9};

res = {}

model = {
  x = {
     dim = 2,
     density = function(x)
        tmp[1] = x[1]/a;
        tmp[2] = x[2]*a + ab*(x[1]*x[1] + aa);
        return(dmvnorm(tmp, zer, c))
   end,
  }
}

para = {
   niter = 50e3, 
   nburn = 10e3,
}

function functional()
  local tmp1, tmp2, tmp3
  tmp1 = x[1]/a;
  tmp2 = x[2]*a + ab*(x[1]*x[1] + aa);
  tmp2 = 2.294157338705618*tmp2 -2.064741604835056*tmp1;
  z = tmp1^2 + tmp2^2;
  for k=1,Npp do res[k] = 1-true_pp[k] end
  for k=1,Npp do 
     if z<=pp[k] then break else res[k] = res[k]-1 end
  end
  return res
end
