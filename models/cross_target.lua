requires_numlua()

m1 = matrix({4,4}); m2 = matrix({0,0})
D = matrix({{0.01, 0},
     {0,2}})
Q1 = matrix({{cos(pi/4), cos(pi*3/4)},
      {sin(pi/4), sin(pi*3/4)}})
Q2 = matrix({{cos(pi*3/4), cos(pi/4)},
      {sin(pi*3/4), sin(pi/4)}})

-- C1 = {{1.0050 , -0.9950},
--      {-0.9950 , 1.0050}}
-- C2 = {{1.0050 , 0.9950},
--      {0.9950 , 1.0050}}
C1 = Q1*D*(#Q1)
C2 = Q2*D*(#Q2)
CC1 = #(matrix.chol(C1));
CC2 = #(matrix.chol(C2));

model = {
  x = {
     dim = 2,
     type = "matrix",
     density = function(x)
                  local p1,p2
      p1 = dmvnorm_chol(x, m1, CC1)
      p2 = dmvnorm_chol(x, m2, CC2)
      return log(exp(p1)+exp(p2))
    end
  }
}

para = {
  niter = 100e3,
  nburn = 10e3,
  algorithm = "am"
}

function functional()
   return {x[1]-2, x[2]-2}
end
