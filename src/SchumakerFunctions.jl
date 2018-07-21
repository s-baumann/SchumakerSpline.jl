"""
Creates a Schumaker spline.
### Takes
* x - A Float64 vector of x coordinates.
* y - A Float64 vector of y coordinates.
* gradients (optional)- A Float64 vector of gradients at each point. If not supplied these are imputed from x and y.
* extrapolation (optional) - This should be a string in ("Curve", "Linear", "Constant") specifying how to interpolate outside of the sample domain. By default it is "curve" which extends out the first and last quadratic curves. The other options are "Linear" which extends the line (from first and last curve) out from the first and last point and "Constant" which extends out the y value at the first and last point.

### Returns
* A spline object
 """
struct Schumaker
    IntStarts_::Array{Float64}
    IntEnds_::Array{Float64}
    coefficient_matrix_::Array{Float64}
    function Schumaker(x::Array{Float64},y::Array{Float64},gradients::Array{Any} = [], extrapolation::String = "Curve")
        if length(gradients) == 0
           gradients = imputeGradients(x,y)
        end
        IntStarts, IntEnds, SpCoefs = getCoefficientMatrix(x,y,gradients, extrapolation)
        new(IntStarts, IntEnds, SpCoefs)
     end
    function Schumaker(x::Array{Int},y::Array{Float64},gradients::Array{Any} = [], extrapolation::String = ("Curve"))
        x_as_floats = convert.(Float64, x)
        return Schumaker(x_as_floats , y , gradients , extrapolation)
    end
    function Schumaker(x::Array{Date},y::Array{Float64},gradients::Array{Any} = [], extrapolation::String = ("Curve"))
        days_as_ints = Dates.days.(x)
        return Schumaker(days_as_ints , y , gradients , extrapolation)
    end
end

"""
Evaluates the spline at a point. The point can be specified as a Float64, Int or Date.
Derivatives can also be taken.
### Takes
 * spline - A Schumaker type spline
 * PointToExamine - The point at which to evaluate the integral
 * derivative - The derivative being sought. This should be 0 to just evaluate the spline, 1 for the first derivative or 2 for a second derivative.
 Higher derivatives are all zero (because it is a quadratic spline). Negative values do not give integrals. Use evaluate_integral instead.
### Returns
 * A Float64 value of the spline or appropriate derivative.
"""
function evaluate(spline::Schumaker,PointToExamine::Float64,  derivative::Int = 0)
    # Derivative of 0 means normal spline evaluation.
    # Derivative of 1, 2 are first and second derivatives respectively.
    IntervalNum = searchsortedlast(spline.IntStarts_, PointToExamine)
    IntervalNum = max(IntervalNum, 1)
    xmt = PointToExamine - spline.IntStarts_[IntervalNum]
    Coefs = spline.coefficient_matrix_[ IntervalNum , :]
    if derivative == 0
        return reshape(Coefs' * [xmt^2 xmt 1]', 1)[1]
    elseif derivative == 1
        return reshape(Coefs' * [2*xmt 1 0]', 1)[1]
    elseif derivative == 2
        return reshape(Coefs' * [2 0 0]', 1)[1]
    elseif derivative < 0
        error("This function cannot do integrals. Use evaluate_integral instead")
    else
        return 0.0
    end
end
function evaluate(spline::Schumaker,PointToExamine::Int,  derivative::Int = 0)
    point_as_float = convert.(Float64, PointToExamine)
    return evaluate(spline,point_as_float,  derivative)
end
function evaluate(spline::Schumaker,PointToExamine::Date,  derivative::Int = 0)
    days_as_int = Dates.days.(PointToExamine)
    return evaluate(spline,days_as_int,  derivative)
end

"""
Estimates the integral of the spline between lhs and rhs. These end points can be input
as Float64s, Ints or Dates.
### Takes
 * spline - A Schumaker type spline
 * lhs - The left hand limit of the integral
 * rhs - The right hand limit of the integral

### Returns
 * A Float64 value of the integral.
"""
function evaluate_integral(spline::Schumaker,lhs::Float64, rhs::Float64)
    first_interval = searchsortedlast(spline.IntStarts_, lhs)
    last_interval = searchsortedlast(spline.IntStarts_, rhs)
    number_of_intervals = last_interval - first_interval
    if number_of_intervals == 0
        return section_integral(spline , lhs, rhs)
    elseif number_of_intervals == 1
        first = section_integral(spline , lhs , spline.IntStarts_[first_interval + 1])
        last  = section_integral(spline , spline.IntStarts_[last_interval]  , rhs)
        return first + last
    else
        interior_areas = 0.0
        first = section_integral(spline , lhs , spline.IntStarts_[first_interval + 1])
        for i in 1:(number_of_intervals-1)
            interior_areas = interior_areas + section_integral(spline ,  spline.IntStarts_[first_interval + i] , spline.IntStarts_[first_interval + i+1] )
        end
        last  = section_integral(spline , spline.IntStarts_[last_interval]  , rhs)
        return first + interior_areas + last
    end
end
function evaluate_integral(spline::Schumaker,lhs::Int, rhs::Int)
    lhs_as_Float = convert(Float64, lhs)
    rhs_as_Float = convert(Float64, rhs)
    return evaluate_integral(spline,lhs_as_Float, rhs_as_Float)
end
function evaluate_integral(spline::Schumaker,lhs::Date, rhs::Date)
    lhs_as_int = Dates.days.(lhs)
    rhs_as_int = Dates.days.(rhs)
    return evaluate_integral(spline,lhs_as_int, rhs_as_int)
end
function section_integral(spline::Schumaker,lhs::Float64,  rhs::Float64)
    # Note that the lhs is used to infer the interval.
    IntervalNum = searchsortedlast(spline.IntStarts_, lhs)
    IntervalNum = max(IntervalNum, 1)
    Coefs = spline.coefficient_matrix_[ IntervalNum , :]
    r_xmt = rhs - spline.IntStarts_[IntervalNum]
    l_xmt = lhs - spline.IntStarts_[IntervalNum]
    return reshape(Coefs' * [(1/3)*r_xmt^3 0.5*r_xmt^2 r_xmt]', 1)[1] - reshape(Coefs' * [(1/3)*l_xmt^3 0.5*l_xmt^2 l_xmt]', 1)[1]
end

"""
Imputes gradients based on a vector of x and y coordinates.
### Takes
 * x - A Float64 vector of x coordinates
 * y - A Float64 vector of y coordinates

### Returns
 * A Float64 vector of gradients for each input point
"""
function imputeGradients(x::Array{Float64},y::Array{Float64})
     n = length(x)
     # Judd (1998), page 233, second last equation
     L = sqrt.( (x[2:n]-x[1:(n-1)]).^2 + (y[2:n]-y[1:(n-1)]).^2)
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
Splits an interval into 2 subintervals and creates the quadratic coefficients
### Takes
 * s - A 2 entry Float64 vector with gradients at either end of the interval
 * z - A 2 entry Float64 vector with y values at either end of the interval
 * Smallt - A 2 entry Float64 vector with x values at either end of the interval

### Returns
 * A 2 x 5 matrix. The first column is the x values of start of the two subintervals. The second column is the ends. The last 3 columns are quadratic coefficients in two subintervals.
"""
function schumakerIndInterval(s::Array{Float64},z::Array{Float64},Smallt::Array{Float64})
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
 function getCoefficientMatrix(x::Array{Float64},y::Array{Float64},gradients::Array{Float64}, extrapolation::String)
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
function extrapolate(fullMatrix::Array{Float64,2}, extrapolation::String, x::Array{Float64}, y::Array{Float64})
  if extrapolation == "Curve"
    return fullMatrix
  end

  dim = size(fullMatrix)[1]
  Botx   = fullMatrix[1,1]
  Boty   = y[1]
  if extrapolation == "Linear"
    BotB = fullMatrix[1 , 4]
    BotC   = Boty - BotB
  else
    BotB = 0.0
    BotC = Boty
  end
  BotRow = [ Botx-1, Botx, 0.0, BotB, BotC]
  Topx = fullMatrix[dim,2]
  Topy = y[length(y)]
  if extrapolation == "Linear"
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
