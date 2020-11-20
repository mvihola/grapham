requires_numlua()

if n == nil then
   n = 5
end
--require "numlua.rng" and use a special function
r = rng(systime_usec())

-- generate random positive definite matrix
c = matrix(n, n);
x0 = matrix(n);
x0:foreach(function(e) x0[e]=r:norm(0, 1) end)
c:foreach(function(i,j) c[i][j]=r:norm(0,1) end)
c = matrix.chol(c * c:transpose())
x0 = c:transpose() * x0;

-- compute the bounds for HPDs
th_true = {0.9,0.75,0.5,0.25,0.1}
th = {}
for i=1,#th_true do
   th[i] = spfun.qchisq(1-th_true[i], n)
end

const = {
   m = matrix(n):fill(0),
   thres = th,
   thres_true = th_true,
}

model = {
   cchol = {
      type = "matrix",
      dim = {n,n},
      init_val = c,
      instantiated = 1,
      density = "duniform",
   },
   x = {
      dim = n,
      init_val = x0,
      type = "matrix",
      density = "dmvnorm_chol",
      parents = {"m","cchol"}
   }
}

para = {
   niter = 400e3,
   nburn = 100e3,
   clib = "clib/thresholds.so"
}

functional = {
   name = "thresholds",
   dim = #th_true,
   args = {"x","m","cchol","thres","thres_true"}
}
