nu = 3
m2 = nu/(nu-2)

model = {
  x = {
    density = function(x)
      return dstudent(x, nu)
    end
  }
}

functional = function()
  return {x, x*x-m2}
end

para = {
  niter = 200e3, nburn = 10e3 
}
