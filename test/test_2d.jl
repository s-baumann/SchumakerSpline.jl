using SchumakerSpline
using Test

gridx = collect(1:10)
gridy = collect(1:10)
grid = gridx * gridy'

schum = Schumaker2d(gridx, gridy, grid; extrapolation = (Linear, Linear))

abs(schum(5,5) - 25) < 100 * eps()
abs(schum(5.5,5.5) - 30.25) < 100 * eps()
