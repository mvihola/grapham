const = {
   v = 0.00434
}
model = {
   mu = {
      density = "duniform"
   },
   t = {
      parents = {"mu","a"}, density = "dnorm"
   },
   y = {
      parents = {"t", "v"}, density = "dnorm"
   },
   a = {
      init_val = 1,
      density = function(a_)
        return dexp(1/a_, 1/2)
      end
   },
}
_, y = read_csv("models/baseball.data")
repeat_block({"y","t"}, y[1])
function functional()
  return {t1, mu, a}
end
para = {
   niter = 30000, nburn = 10000, algorithm = "asm",
}
