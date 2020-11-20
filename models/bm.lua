const = {
   x0 = 0, v = 1
}
model = {
   x = {
      parents = {"x-", "v"},
      density = "dnorm"
   }
}

repeat_block({"x"},10)

functional = {
   name = "mom2",
   args = {"x10"},
   dim = 1
}