

tri_fun = function(s){
  0.5 - 0.5*exp(-abs(s))
}

tri_fun(1)


s = seq(0.1, 3, 0.01)

plot(s, tri_fun(9*s))

#plot(s, 1/2*tanh(s^3))