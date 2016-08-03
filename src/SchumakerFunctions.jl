"""
Imputes gradients based on a vector of x and y coordinates.
### Takes
 * x - A Float64 vector of x coordinates
 * y - A Float64 vector of y coordinates
### Returns
 * A Float64 vector of gradients for each input point
"""
function imputeGradients(x,y)
     n = length(x)
     # Judd (1998), page 233, second last equation
     L = sqrt( (x[2:n]-x[1:(n-1)]).^2 + (y[2:n]-y[1:(n-1)]).^2)
     # Judd (1998), page 233, last equation
     d = (y[2:n]-y[1:(n-1)])./(x[2:n]-x[1:(n-1)])
     # Judd (1998), page 234, Eqn 6.11.6
     Conditionsi = d[1:(n-2)].*d[2:(n-1)] .> 0
     MiddleSiwithoutApplyingCondition = (L[1:(n-2)].*d[1:(n-2)]+L[2:(n-1)].* d[2:(n-1)]) ./ (L[1:(n-2)]+L[2:(n-1)])
     sb = Conditionsi .* MiddleSiwithoutApplyingCondition
     # Judd (1998), page 234, Second Equation line plus 6.11.6 gives this array of slopes.
     ff = [((-sb[1]+3*d[1])/2);  sb ;  ((3*d[n-1]-sb[n-2])/2)]
     return ff
 end

 """
 Creates a spline defined by interval starts IntStarts and quadratic coefficients SpCoefs which evaluates an input point.
 ### Takes
  * IntStarts - A Float64 vector that gives the starting points of intervals (in the x plane)
  * SpCoefs - A 3 column matrix with the same number of rows as the length of the IntStarts vector. The first column is the coefficient of the quadratic term, the second column for the linear term. The third column is the constant.
 ### Returns
  * A spline function that takes a single Float64 input and returns the spline value at that point.
 """
function ppmak(IntStarts,SpCoefs)
  function sp(PointToExamine)
    IntervalNum = searchsortedlast(IntStarts, PointToExamine)
    IntervalNum = max(IntervalNum, 1)
    xmt = PointToExamine - IntStarts[IntervalNum]
    Coefs = SpCoefs[ IntervalNum , :]
    return reshape(Coefs * [xmt^2 xmt 1]', 1)[1]
  end
  function Vsp(PointToExamine)
    return map(x -> sp(x), PointToExamine)
  end
  return Vsp
end


 """
 Creates the derivative function of the spline defined by interval starts IntStarts and quadratic coefficients SpCoefs which evaluates an input point.
 ### Takes
  * IntStarts - A Float64 vector that gives the starting points of intervals (in the x plane)
  * SpCoefs - A 3 column matrix with the same number of rows as the length of the IntStarts vector. The first column is the coefficient of the quadratic term, the second column for the linear term. The third column is the constant.
 ### Returns
 * The derivative function that takes a single Float64 input and returns the derivative at that point.
 """
function ppmakDeriv(IntStarts,SpCoefs)
  function sp(PointToExamine)
    IntervalNum = searchsortedlast(IntStarts, PointToExamine)
    IntervalNum = max(IntervalNum, 1)
    xmt = PointToExamine - IntStarts[IntervalNum]
    Coefs = SpCoefs[ IntervalNum , :]
    return reshape(Coefs * [2*xmt 1 0]', 1)[1]
  end
  function Vsp(PointToExamine)
    return map(x -> sp(x), PointToExamine)
  end
  return Vsp
end

"""
Creates the second derivative function of the spline defined by interval starts IntStarts and quadratic coefficients SpCoefs which evaluates an input point.
### Takes
 * IntStarts - A Float64 vector that gives the starting points of intervals (in the x plane)
 * SpCoefs - A 3 column matrix with the same number of rows as the length of the IntStarts vector. The first column is the coefficient of the quadratic term, the second column for the linear term. The third column is the constant.
### Returns
 * The second derivative function that takes a single Float64 input and returns the second derivative at that point.
"""
function ppmak2Deriv(IntStarts,SpCoefs)
  function sp(PointToExamine)
    IntervalNum = searchsortedlast(IntStarts, PointToExamine)
    IntervalNum = max(IntervalNum, 1)
    xmt = PointToExamine - IntStarts[IntervalNum]
    Coefs = SpCoefs[ IntervalNum , :]
    return reshape(Coefs * [2 0 0]', 1)[1]
  end
  function Vsp(PointToExamine)
    return map(x -> sp(x), PointToExamine)
  end
  return Vsp
end

"""
Splits an interval into 2 subintervals and creates the quadratic coefficients
### Takes
 * s - A 2 entry Float64 vector with gradients at either end of the interval
 * z - A 2 entry Float64 vector with y values at either end of the interval
 * Smallt - A 2 entry Float64 vector with x values at either end of the interval
### Returns
 * A 2 x 5 matrix. The first column is the x values of start of the two subintervals. The second column is the ends. The last 3 columns are quadratic coefficients in two subintervals.
"""
function schumakerIndInterval(s,z,Smallt)
   # The SchumakerIndInterval function takes in each interval individually
   # and returns the location of the knot as well as the quadratic coefficients in each subinterval.

   # Judd (1998), page 232, Lemma 6.11.1 provides this if condition:
   if (sum(s)*(Smallt[2]-Smallt[1]) == 2*(z[2]-z[1]))
     tsi = Smallt[2]
   else
     # Judd (1998), page 233, Algorithm 6.3 along with equations 6.11.4 and 6.11.5 provide this whole section
     delta = (z[2] -z[1])/(Smallt[2]-Smallt[1])
     Condition = ((s[1]-delta)*(s[2]-delta) >= 0)
     Condition2 = abs(s[2]-delta) < abs(s[1]-delta)
     if (Condition)
       tsi = sum(Smallt)/2
     elseif (Condition2)
       tsi = (Smallt[1] + (Smallt[2]-Smallt[1])*(s[2]-delta)/(s[2]-s[1]))
     else
       tsi = (Smallt[2] + (Smallt[2]-Smallt[1])*(s[1]-delta)/(s[2]-s[1]))
     end
   end

   # Judd (1998), page 232, 3rd last equation of page.
   alpha = tsi-Smallt[1]
   beta = Smallt[2]-tsi
   # Judd (1998), page 232, 4th last equation of page.
   sbar = (2*(z[2]-z[1])-(alpha*s[1]+beta*s[2]))/(Smallt[2]-Smallt[1])
   # Judd (1998), page 232, 3rd equation of page. (C1, B1, A1)
   Coeffs1 = [ (sbar-s[1])/(2*alpha)  s[1]  z[1] ]
   if (beta == 0)
     Coeffs2 = Coeffs1
   else
     # Judd (1998), page 232, 4th equation of page. (C2, B2, A2)
     Coeffs2 = [ (s[2]-sbar)/(2*beta)  sbar  Coeffs1 * [alpha^2, alpha, 1] ]
   end
   Machine4Epsilon = 4*eps()
     if (tsi  <  Smallt[1] + Machine4Epsilon )
         return [Smallt[1] Smallt[2] Coeffs2]
     elseif (tsi + Machine4Epsilon > Smallt[2] )
         return [Smallt[1] Smallt[2] Coeffs1]
     else
         return [Smallt[1] tsi Coeffs1 ; tsi Smallt[2] Coeffs2]
     end
 end

 """
 Calls SchumakerIndInterval many times to get full set of spline intervals and coefficients. Then calls extrapolation for out of sample behaviour
 ### Takes
 * gradients - A Float64 vector of gradients at each point
 * x - A Float64 vector of x coordinates
 * y - A Float64 vector of y coordinates
 * extrapolation - A string in ("Curve", "Linear", "Constant") that gives behaviour outside of interpolation range.
 ### Returns
  * A vector of interval starts
  * A vector of interval ends
  * A matrix of all coefficients
  """
 function getCoefficientMatrix(gradients,y,x, extrapolation)
   n = length(x)
   fullMatrix = schumakerIndInterval([gradients[1] gradients[2]], [y[1] y[2]], [x[1] x[2]] )
    for intrval = 2:(n-1)
      Smallt = [ x[intrval] , x[intrval + 1] ]
      s = [ y[intrval], y[intrval + 1] ]
      z = [ gradients[intrval], gradients[intrval + 1] ]
      intMatrix = schumakerIndInterval(z,s,Smallt)
      fullMatrix = vcat(fullMatrix,intMatrix)
    end
    fullMatrix = extrapolate(fullMatrix, extrapolation, x, y)
   return fullMatrix[:,1], fullMatrix[:,2], fullMatrix[:,3:5]
 end

 """
 Adds a row on top and bottom of coefficient matrix to give out of sample prediction.
 ### Takes
 * fullMatrix - output from GetCoefficientMatrix first few lines
 * extrapolation - A string in ("Curve", "Linear", "Constant") that gives behaviour outside of interpolation range.
 * x - A Float64 vector of x coordinates
 * y - A Float64 vector of y coordinates
 ### Returns
  * A new version of fullMatrix with out of sample prediction built into it.
  """
function extrapolate(fullMatrix, extrapolation, x, y)

  if (extrapolation == "Curve")
    return fullMatrix
  end

  dim = size(fullMatrix)[1]

  Botx   = fullMatrix[1,1]
  Boty   = y[1]

  if (extrapolation == "Linear")
    BotB = fullMatrix[1 , 4]
    BotC   = Boty - BotB
  else
    BotB = 0.0
    BotC = Boty
  end

  BotRow = [ Botx-1, Botx, 0.0, BotB, BotC]

  Topx = fullMatrix[dim,2]
  Topy = y[length(y)]

  if (extrapolation == "Linear")
    TopB = fullMatrix[dim ,4]
    TopC = Topy
  else
    TopB = 0.0
    TopC = Topy
  end

  TopRow = [ Topx, Topx + 1, 0.0 ,TopB ,TopC]

  fullMatrix = vcat(BotRow' , fullMatrix,  TopRow')

  return fullMatrix
end

"""
Creates splines for a given set of x and y values (and optionally gradients) and the first and second derivatives of this spline.
### Takes
* x - A Float64 vector of x coordinates.
* y - A Float64 vector of y coordinates.
* gradients (optional)- A Float64 vector of gradients at each point. If not supplied these are imputed from x and y.
* extrapolation (optional) - This should be a string in ("Curve", "Linear", "Constant") specifying how to interpolate outside of the sample domain. By default it is "curve" which extends out the first and last quadratic curves. The other options are "Linear" which extends the line (from first and last curve) out from the first and last point and "Constant" which extends out the y value at the first and last point.
### Returns
* A spline which takes and input value and returns the spline y value.
* The derivative of this spline.
* The second derivative of this spline.
 """
function schumaker(x,y,gradients = "Not-Supplied", extrapolation = ("Curve", "Linear", "Constant"))
  # This is the main function of the package that creates and returns a schumaker spline.
  # Inputs : An x array, the corresponding y array and optionally the gradients at each point.
  #           If gradients are not input then they are estimated.
  # Outputs: A schumaker spline function. Its derivative and its second derivative.

  if (extrapolation == ("Curve", "Linear", "Constant"))
    extrapolation = "Curve"
  end

  if (gradients == "Not-Supplied")
     gradients = imputeGradients(x,y)
  end

  IntStarts, IntEnds, SpCoefs = getCoefficientMatrix(gradients,y,x, extrapolation)

  Sp      = ppmak(IntStarts,SpCoefs)
  SpDeriv = ppmakDeriv(IntStarts,SpCoefs)
  SpDeriv2= ppmak2Deriv(IntStarts,SpCoefs)
  return Sp, SpDeriv, SpDeriv2
end
