using SchumakerSpline
using LinearAlgebra

x_row = [0.1,0.3]
x_col = [0.2,0.3,0.5]
ygrid = [1.0 2.0 3.0;
         18.0 -2.0 3.3]
sp = Schumaker2d(x_row,x_col,ygrid; extrapolation = (Linear, Linear))
sp(0.4, 0.3)
