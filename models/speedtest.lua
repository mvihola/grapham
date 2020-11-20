requires_numlua()

n = 300;

const = {
   c = matrix(n):fill(1),
   m = matrix(n):fill(0)
}

model = {
   x = {
      dim = n,
      parents = {"m","c"},
      type = "matrix",
      density = "dmvnorm_diag"
   }
}

para = {
   niter = 10e3,
   algorithm = "asm"
}
