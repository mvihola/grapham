const = {
   c = 15
}

model = {
   x = {
      dim = 2,
      density = "dcircular",
      parents = {"c"}
   }
}

para = {
   niter = 50e3,
   nburn = 10e3,
   algorithm = "ascm",
   blocking = "full",
   clib = "clib/dcircular.so"
--   outfile = "circular.csv",
}

functional = {
   name = "mean",
   dim = 2,
   args = {"x"}
}
