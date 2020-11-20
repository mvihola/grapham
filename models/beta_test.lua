--a = .5; b = .5;
a = 2; b = 5;

mean = a/(a+b);
m2 = mean^2-a*b/((a+b)^2*(a+b+1))

model = {
  x = {
    init_val = .5,
    density = function(x)
      return dbeta(x, a, b)
    end
  }
}

functional = function() 
    return {x-mean,x^2-m2}
end

para = {
  niter = 500e3, nburn = 10e3, acc_opt = 0.44
}
