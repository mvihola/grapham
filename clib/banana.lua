const = {
   a = 1, 
   b = 1, 
   m = {0,0},
   -- The Cholesky factor of
   --   [ 1.0 0.9 ]
   --   [ 0.9 1.0 ]
   c = {{1, 0.9},
        {0, 0.435889894354067}},
}

model = {
  x = {
     dim = 2,
     density = "dbanana",
     parents = {"m","c","a","b"}
  }
}

para = {
   niter = 50e3, 
   nburn = 10e3,
   clib = "clib/dbanana.so",
}

functional = {
   name = "banana_fun",
   dim = 5,
   args = {"x", "c", "a", "b"}
}
