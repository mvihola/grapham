model = {
  x = {
     dim = 2,
     density = "dcross"
  }
}

para = {
  niter = 100e3,
  nburn = 10e3,
  algorithm = "am",
  clib = "clib/dcross.so"
}

functional = {
   name = "cross_fun",
   dim = 2,
   args = {"x"}
}
