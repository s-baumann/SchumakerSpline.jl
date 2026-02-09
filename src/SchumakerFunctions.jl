"""
This creates an enum which details how extrapolation from the interpolation domain should be done.
The possible enum values are:

* `Curve` - Curve extrapolation extends out the quadratic form at the edges. This can lead to a nonmonotonic result (as the curve can eventually change direction)
* `Linear` - Linear extrapolation extends out the gradient from at the edges. This will always lead to a monotonic result.
* `Constant` - Linear extrapolation extends out the value from at the edges. This leads to flat values being extended out.
"""
@enum Schumaker_ExtrapolationSchemes begin
    Curve = 0
    Linear = 1
    Constant = 2
end


"""
    Schumaker(x::Array{T,1},y::Array{T,1} ; gradients::Union{Missing,Array{T,1}} = missing, extrapolation::Schumaker_ExtrapolationSchemes = Curve,
                  left_gradient::Union{Missing,T} = missing, right_gradient::Union{Missing,T} = missing)
    Schumaker(x::Array{Int,1},y::Array{T,1} ; gradients::Union{Missing,Array{T,1}} = missing, extrapolation::Schumaker_ExtrapolationSchemes = Curve,
                  left_gradient::Union{Missing,T} = missing, right_gradient::Union{Missing,T} = missing)
    Schumaker(x::Array{Date,1},y::Array{T,1} ; gradients::Union{Missing,Array{T,1}} = missing, extrapolation::Schumaker_ExtrapolationSchemes = Curve,
                  left_gradient::Union{Missing,T} = missing, right_gradient::Union{Missing,T} = missing)
Creates a Schumaker spline.
### Inputs
* `x` - A vector of x coordinates.
* `y` - A vector of y coordinates.
* `extrapolation` - This should be `Curve`, `Linear` or `Constant` specifying how to interpolate outside of the sample domain.
* `gradients` - A vector of gradients at each point. If not supplied these are imputed from x and y.
* `left_gradient` - The gradient at the lowest value of x in the domain. This will override the gradient imputed or submitted in the gradients optional argument (if it is submitted there)
* `right_gradient` - The gradient at the highest value of x in the domain. This will override the gradient imputed or submitted in the gradients optional argument (if it is submitted there)

### Returns
* A `Schumaker` object which contains the spline. This object can then be evaluated with evaluate or evaluate_integral.
 """
struct Schumaker{T<:AbstractFloat}
    IntStarts_::Array{T,1}
    coefficient_matrix_::Array{T,2}
    function Schumaker(x::Array{T,1},y::Array{<:Real,1} ; gradients::Union{Missing,Array{<:Real,1}} = missing, left_gradient::Union{Missing,<:Real} = missing, right_gradient::Union{Missing,<:Real} = missing,
                       extrapolation::Tuple{Schumaker_ExtrapolationSchemes,Schumaker_ExtrapolationSchemes} = (Curve,Curve)) where T<:Real
        if length(x) == 0
            error("Zero length x vector is insufficient to create Schumaker Spline.")
        elseif length(x) == 1
            G = T<:AbstractFloat ? T : Float64
            IntStarts = Array{G,1}(x)
            @inbounds SpCoefs = G[0 0 y[1]]
            # Note that this hardcodes in constant extrapolation. This is only
            # feasible one as we do not have derivative or curve information.
            return new{G}(IntStarts, SpCoefs)
        elseif length(x) == 2
            G = T<:AbstractFloat ? T : Float64
            @inbounds IntStarts = Array{G,1}([x[1]])
            @inbounds linear_coefficient = (y[2]- y[1]) / (x[2]-x[1])
            @inbounds SpCoefs = G[0 linear_coefficient y[1]]
            # In this case it defaults to curve extrapolation (which is same as linear here)
            # So we just alter in case constant is specified.
            if extrapolation[1] == Constant || extrapolation[2] == Constant
                matrix_without_extrapolation = hcat(IntStarts , SpCoefs)
                @inbounds matrix_with_extrapolation    = extrapolate(matrix_without_extrapolation, extrapolation, x[2], y)
                return @inbounds new{G}(G.(matrix_with_extrapolation[:,1]), G.(matrix_with_extrapolation[:,2:4]))
            else
                return new{G}(IntStarts, SpCoefs)
            end
        end
        if ismissing(gradients)
           gradients = imputeGradients(x,y)
        end
        if !ismissing(left_gradient)
            @inbounds gradients[1] = left_gradient
        end
        if !ismissing(right_gradient)
            @inbounds gradients[length(gradients)] = right_gradient
        end
        IntStarts, SpCoefs = getCoefficientMatrix(x,y,gradients, extrapolation)
        if T<:AbstractFloat
            G = T
        else
            G = Float64
        end
        return new{G}(G.(IntStarts), G.(SpCoefs))
     end
    function Schumaker(IntStarts_::Array{T,1}, coefficient_matrix_::Array{T,2}) where {T<:Real}
        return new{T}(IntStarts_, coefficient_matrix_)
    end
    function Schumaker{Q}(IntStarts_::Array{<:Real,1}, coefficient_matrix_::Array{<:Real,2}) where Q<:Real
        return new{Q}(Q.(IntStarts_), Q.(coefficient_matrix_))
    end
    function Schumaker(x::Array{Date,1},y::Array{<:Real,1} ; gradients::Union{Missing,Array{<:Real,1}} = missing, left_gradient::Union{Missing,<:Real} = missing, right_gradient::Union{Missing,<:Real} = missing,
                       extrapolation::Tuple{Schumaker_ExtrapolationSchemes,Schumaker_ExtrapolationSchemes} = (Curve,Curve))
        days_as_ints = Dates.days.(x)
        T = promote_type(eltype(y), eltype(days_as_ints))
        if T<:AbstractFloat
            G = T
        else
            G = Float64
        end
        return Schumaker{G}(days_as_ints , y; gradients = gradients , extrapolation = extrapolation, left_gradient = left_gradient, right_gradient = right_gradient)
    end
    function Schumaker{Q}(x::AbstractArray, y::AbstractArray ; gradients::Union{Missing,AbstractArray} = missing, left_gradient::Union{Missing,Real} = missing, right_gradient::Union{Missing,Real} = missing,
                          extrapolation::Tuple{Schumaker_ExtrapolationSchemes,Schumaker_ExtrapolationSchemes} = (Curve,Curve)) where Q<:Real
        got_both = (!).(ismissing.(x) .| ismissing.(y))
        @inbounds new_x = convert.(Q, x[got_both])
        @inbounds new_y = convert.(Q, y[got_both])
        @inbounds new_gradients = ismissing(gradients)        ? missing : convert.(Q, gradients[got_both])
        converted_left  = ismissing(left_gradient)  ? missing : convert(Q, left_gradient)
        converted_right = ismissing(right_gradient) ? missing : convert(Q, right_gradient)
        new_left  = got_both[1]                     ? converted_left  : missing
        new_right = got_both[length(got_both)]      ? converted_right : missing
        if length(new_x) == 0 error("After removing missing elements there are no points left to estimate schumaker spline") end
        return Schumaker(new_x , new_y; gradients = new_gradients , extrapolation = extrapolation, left_gradient = new_left, right_gradient = new_right)
    end
    function Schumaker(x::Union{AbstractArray{T,1},AbstractArray{Union{Missing,T},1}},y::Union{AbstractArray{R,1},AbstractArray{Union{Missing,R},1}} ; gradients::Union{Missing,AbstractArray{<:Real,1}} = missing, left_gradient::Union{Missing,Real} = missing, right_gradient::Union{Missing,Real} = missing,
                       extrapolation::Tuple{Schumaker_ExtrapolationSchemes,Schumaker_ExtrapolationSchemes} = (Curve,Curve)) where T<:Real where R<:Real
       promo_type = promote_type(T,R)
       return Schumaker{promo_type}(x, y; gradients = gradients, extrapolation = extrapolation, left_gradient = left_gradient, right_gradient = right_gradient)
   end
end
Base.broadcastable(e::Schumaker) = Ref(e)

"""
Evaluates the spline at a point. The point can be specified as a Real number (Int, Float, etc) or a Date.
Derivatives can also be taken.
### Inputs
 * `spline` - A `Schumaker` type spline
 * `PointToExamine` - The point at which to evaluate the integral
 * `derivative` - The derivative being sought. This should be 0 to just evaluate the spline, 1 for the first derivative or 2 for a second derivative.
 Higher derivatives are all zero (because it is a quadratic spline). Negative values do not give integrals. Use evaluate_integral instead.
### Returns
 * A value of the spline or appropriate derivative in the same format as specified in the spline.
"""
function evaluate(spline::Schumaker, PointToExamine, derivative::Integer = 0)
    IntervalNum = searchsortedlast(spline.IntStarts_, PointToExamine)
    IntervalNum = max(IntervalNum, 1)
    @inbounds xmt = PointToExamine - spline.IntStarts_[IntervalNum]
    @inbounds a = spline.coefficient_matrix_[IntervalNum, 1]
    @inbounds b = spline.coefficient_matrix_[IntervalNum, 2]
    @inbounds c = spline.coefficient_matrix_[IntervalNum, 3]
    if derivative == 0
        return a * xmt * xmt + b * xmt + c
    elseif derivative == 1
        return 2 * a * xmt + b
    elseif derivative == 2
        return 2 * a
    elseif derivative < 0
        error("This function cannot do integrals. Use evaluate_integral instead")
    else
        return zero(c)
    end
end
function evaluate(spline::Schumaker, PointToExamine::Date,  derivative::Int = 0)
    days_as_int = Dates.days.(PointToExamine)
    return evaluate(spline,days_as_int,  derivative)
end
function (s::Schumaker)(x, deriv::Integer = 0)
    return evaluate(s, x, deriv)
end


"""
Estimates the integral of the spline between lhs and rhs. These end points can be input
as Reals or Dates.
### Inputs
 * `spline` - A Schumaker type spline
 * `lhs` - The left hand limit of the integral
 * `rhs` - The right hand limit of the integral
### Returns
 * A `Float64` value of the integral.
"""
function evaluate_integral(spline::Schumaker, lhs::Real, rhs::Real)
    first_interval = searchsortedlast(spline.IntStarts_, lhs)
    first_interval = max(first_interval, 1)
    last_interval = searchsortedlast(spline.IntStarts_, rhs)
    last_interval = max(last_interval, 1)
    number_of_intervals = last_interval - first_interval
    if number_of_intervals == 0
        return _section_integral(spline, first_interval, lhs, rhs)
    elseif number_of_intervals == 1
        @inbounds first = _section_integral(spline, first_interval, lhs, spline.IntStarts_[first_interval + 1])
        @inbounds lst   = _section_integral(spline, last_interval, spline.IntStarts_[last_interval], rhs)
        return first + lst
    else
        interior_areas = 0.0
        @inbounds first = _section_integral(spline, first_interval, lhs, spline.IntStarts_[first_interval + 1])
        @simd for i in 1:(number_of_intervals-1)
            @inbounds sec_int = _section_integral(spline, first_interval + i, spline.IntStarts_[first_interval + i], spline.IntStarts_[first_interval + i + 1])
            interior_areas += sec_int
        end
        @inbounds lst = _section_integral(spline, last_interval, spline.IntStarts_[last_interval], rhs)
        return first + interior_areas + lst
    end
end
function evaluate_integral(spline::Schumaker, lhs::Date, rhs::Date)
    return evaluate_integral(spline, Dates.days.(lhs) , Dates.days.(rhs))
end
function _section_integral(spline::Schumaker, IntervalNum::Int, lhs::Real, rhs::Real)
    @inbounds a = spline.coefficient_matrix_[IntervalNum, 1]
    @inbounds b = spline.coefficient_matrix_[IntervalNum, 2]
    @inbounds c = spline.coefficient_matrix_[IntervalNum, 3]
    @inbounds start = spline.IntStarts_[IntervalNum]
    r = rhs - start
    l = lhs - start
    return a * (r^3 - l^3) / 3 + b * (r^2 - l^2) / 2 + c * (r - l)
end

function section_integral(spline::Schumaker, lhs::Real, rhs::Real)
    IntervalNum = searchsortedlast(spline.IntStarts_, lhs)
    IntervalNum = max(IntervalNum, 1)
    return _section_integral(spline, IntervalNum, lhs, rhs)
end

"""
find_derivative_spline(spline::Schumaker)
    Returns a SchumakerSpline that is the derivative of the input spline
"""
function find_derivative_spline(spline::Schumaker{V}) where V
    coefficient_matrix = Array{V,2}(undef, size(spline.coefficient_matrix_, 1), 3)
    coefficient_matrix[:,3] .= spline.coefficient_matrix_[:,2]
    coefficient_matrix[:,2] .= 2 .* spline.coefficient_matrix_[:,1]
    coefficient_matrix[:,1] .= 0.0
    return Schumaker(spline.IntStarts_, coefficient_matrix)
end



"""
    imputeGradients(x::Vector{T}, y::Vector{T})
Imputes gradients based on a vector of x and y coordinates.
"""
function imputeGradients(x::Vector{<:Real}, y::Vector{<:Real})
     n = length(x)
     @views begin
         # Judd (1998), page 233, second last equation
         L = sqrt.((x[2:n] .- x[1:(n-1)]).^2 .+ (y[2:n] .- y[1:(n-1)]).^2)
         # Judd (1998), page 233, last equation
         d = (y[2:n] .- y[1:(n-1)]) ./ (x[2:n] .- x[1:(n-1)])
         # Judd (1998), page 234, Eqn 6.11.6
         Conditionsi = d[1:(n-2)] .* d[2:(n-1)] .> 0
         MiddleSiwithoutApplyingCondition = (L[1:(n-2)] .* d[1:(n-2)] .+ L[2:(n-1)] .* d[2:(n-1)]) ./ (L[1:(n-2)] .+ L[2:(n-1)])
     end
     sb = Conditionsi .* MiddleSiwithoutApplyingCondition
     # Judd (1998), page 234, Second Equation line plus 6.11.6 gives this array of slopes.
     ff = Vector{Float64}(undef, n)
     @inbounds ff[1] = (-sb[1] + 3 * d[1]) / 2
     @inbounds ff[n] = (3 * d[n-1] - sb[n-2]) / 2
     @inbounds for i in 1:length(sb)
         ff[i+1] = sb[i]
     end
     return ff
 end

"""
Writes 1 or 2 rows into `buffer` starting at row `row`. Returns the number of rows written (1 or 2).
Buffer columns are [x_start, a, b, c].
"""
function schumakerIndInterval!(buffer::Matrix{<:Real}, row::Int, g1::Real, g2::Real, y1::Real, y2::Real, x1::Real, x2::Real)
   # Judd (1998), page 232, Lemma 6.11.1 provides this if condition:
   if ((g1 + g2) * (x2 - x1) == 2 * (y2 - y1))
     tsi = x2
   else
     # Judd (1998), page 233, Algorithm 6.3 along with equations 6.11.4 and 6.11.5
     delta = (y2 - y1) / (x2 - x1)
     Condition = ((g1 - delta) * (g2 - delta) >= 0)
     Condition2 = abs(g2 - delta) < abs(g1 - delta)
     if Condition
       tsi = (x1 + x2) / 2
     elseif Condition2
       tsi = x1 + (x2 - x1) * (g2 - delta) / (g2 - g1)
     else
       tsi = x2 + (x2 - x1) * (g1 - delta) / (g2 - g1)
     end
   end

   # Judd (1998), page 232, 3rd last equation of page.
   alpha = tsi - x1
   beta = x2 - tsi
   # Judd (1998), page 232, 4th last equation of page.
   sbar = (2 * (y2 - y1) - (alpha * g1 + beta * g2)) / (x2 - x1)
   # Judd (1998), page 232, 3rd equation of page. (C1, B1, A1)
   c1_a = (sbar - g1) / (2 * alpha)
   c1_b = g1
   c1_c = y1
   # Value at tsi: computed without intermediate division to avoid Inf*0 = NaN when alpha ≈ 0.
   # Mathematically: c1_a * alpha^2 + c1_b * alpha + c1_c = (sbar - g1)*alpha/2 + g1*alpha + y1
   val_at_tsi = (sbar - g1) * alpha / 2 + g1 * alpha + y1
   if beta == 0
     c2_a = c1_a
     c2_b = c1_b
     c2_c = c1_c
   else
     # Judd (1998), page 232, 4th equation of page. (C2, B2, A2)
     c2_a = (g2 - sbar) / (2 * beta)
     c2_b = sbar
     c2_c = val_at_tsi
   end
   Machine4Epsilon = 4 * eps()
   if tsi < x1 + Machine4Epsilon
       buffer[row, 1] = x1
       buffer[row, 2] = c2_a
       buffer[row, 3] = c2_b
       buffer[row, 4] = c2_c
       return 1
   elseif tsi + Machine4Epsilon > x2
       buffer[row, 1] = x1
       buffer[row, 2] = c1_a
       buffer[row, 3] = c1_b
       buffer[row, 4] = c1_c
       return 1
   else
       buffer[row, 1] = x1
       buffer[row, 2] = c1_a
       buffer[row, 3] = c1_b
       buffer[row, 4] = c1_c
       buffer[row + 1, 1] = tsi
       buffer[row + 1, 2] = c2_a
       buffer[row + 1, 3] = c2_b
       buffer[row + 1, 4] = c2_c
       return 2
   end
end

function schumakerIndInterval(gradients::Vector{<:Real}, y::Vector{<:Real}, x::Vector{<:Real})
   @inbounds g1 = gradients[1]; @inbounds g2 = gradients[2]
   @inbounds y1 = y[1]; @inbounds y2 = y[2]
   @inbounds x1 = x[1]; @inbounds x2 = x[2]
   buffer = Matrix{Float64}(undef, 2, 4)
   nrows = schumakerIndInterval!(buffer, 1, g1, g2, y1, y2, x1, x2)
   return buffer[1:nrows, :]
end

 """
 Calls `SchumakerIndInterval` many times to get full set of spline intervals and coefficients. Then calls extrapolation for out of sample behaviour
### Inputs
 * `gradients` - A vector of gradients at each point
 * `x` - A vector of `x` coordinates
 * `y` - A vector of `y` coordinates
 * `extrapolation` - A string in (Curve, Linear or Constant) that gives behaviour outside of interpolation range.
### Returns
 * A vector of interval starts
 * A vector of interval ends
 * A matrix of all coefficients
  """
 function getCoefficientMatrix(x::Array{<:Real,1}, y::Array{<:Real,1}, gradients::Array{<:Real,1},
                               extrapolation::Tuple{Schumaker_ExtrapolationSchemes,Schumaker_ExtrapolationSchemes})
   n = length(x)
   # Each of (n-1) intervals produces at most 2 rows. Extrapolation adds at most 2 more.
   max_rows = 2 * (n - 1) + 2
   buffer = Matrix{Float64}(undef, max_rows, 4)
   # Leave row 1 open for possible left extrapolation row
   data_start = extrapolation[1] != Curve ? 2 : 1
   current_row = data_start
   @inbounds nr = schumakerIndInterval!(buffer, current_row, gradients[1], gradients[2], y[1], y[2], x[1], x[2])
   current_row += nr
   for intrval = 2:(n-1)
       @inbounds nr = schumakerIndInterval!(buffer, current_row, gradients[intrval], gradients[intrval+1], y[intrval], y[intrval+1], x[intrval], x[intrval+1])
       current_row += nr
   end
   # data rows are in buffer[data_start : current_row-1, :]
   last_data_row = current_row - 1

   # Left extrapolation: write into row 1 (before data rows)
   if extrapolation[1] != Curve
       Botx = buffer[data_start, 1]
       Boty = y[1]
       if extrapolation[1] == Linear
           BotB = buffer[data_start, 3]
           BotC = Boty - BotB * 1e-10
       else # Constant
           BotB = 0.0
           BotC = Boty
       end
       buffer[1, 1] = Botx - 1e-10
       buffer[1, 2] = 0.0
       buffer[1, 3] = BotB
       buffer[1, 4] = BotC
   end

   # Right extrapolation: write into current_row (after data rows)
   # When both sides are Curve, no extrapolation rows are needed (quadratic extends naturally).
   # The right Curve row is only added when left is NOT Curve (mixed mode).
   if !(extrapolation[1] == Curve && extrapolation[2] == Curve)
       Topx = x[n]
       Topy = y[n]
       top_row = current_row
       if extrapolation[2] == Linear
           last_interval_width = Topx - buffer[last_data_row, 1]
           grad_at_right = 2 * buffer[last_data_row, 2] * last_interval_width + buffer[last_data_row, 3]
           buffer[top_row, 1] = Topx
           buffer[top_row, 2] = 0.0
           buffer[top_row, 3] = grad_at_right
           buffer[top_row, 4] = Topy
           current_row += 1
       elseif extrapolation[2] == Constant
           buffer[top_row, 1] = Topx
           buffer[top_row, 2] = 0.0
           buffer[top_row, 3] = 0.0
           buffer[top_row, 4] = Topy
           current_row += 1
       elseif extrapolation[2] == Curve
           Gap = Topx - buffer[last_data_row, 1]
           TopA = buffer[last_data_row, 2]
           TopB = buffer[last_data_row, 3] + 2 * TopA * Gap
           TopC = buffer[last_data_row, 4] + TopB * Gap - TopA * Gap * Gap
           buffer[top_row, 1] = Topx
           buffer[top_row, 2] = TopA
           buffer[top_row, 3] = TopB
           buffer[top_row, 4] = TopC
           current_row += 1
       end
   end

   # Determine final range
   first_row = extrapolation[1] != Curve ? 1 : data_start
   last_row = current_row - 1
   result = buffer[first_row:last_row, :]
   return result[:, 1], result[:, 2:4]
 end

"""
 Adds a row on top and bottom of coefficient matrix to give out of sample prediction.
### Inputs
 * `fullMatrix` - output from `GetCoefficientMatrix` first few lines
 * `extrapolation` - A tuple with two enums in (Curve, Linear, Constant) that gives behaviour outside of interpolation range.
 * `x` - A vector of x coordinates
 * `y` - A vector of y coordinates
### Returns
  * A new version of fullMatrix with out of sample prediction built into it.
"""
function extrapolate(fullMatrix::Array{<:Real,2}, extrapolation::Tuple{Schumaker_ExtrapolationSchemes,Schumaker_ExtrapolationSchemes}, Topx::Real, y::Array{<:Real,1})
    # In the fullMatrix the first column are the x values and then the next 3 columns are a, b, c in the expression a(x-start)^2 + b(x-start) + c
  if (extrapolation[1] == Curve) && (extrapolation[2] == Curve)
    return fullMatrix
  end
  # Preliminaries used throughout
  dim = size(fullMatrix, 1)
  Botx   = fullMatrix[1,1]
  Boty   = y[1]
  Topy = y[length(y)]
  # Initialising variable so their domain is not restricted to the if statement spaces.
  BotB = 0.0
  BotC = 0.0
  TopA = 0.0
  TopB = 0.0
  TopC = 0.0
  # Now creating the extrapolation to the left.
  if extrapolation[1] == Linear
    BotB = fullMatrix[1 , 3]
    BotC   = Boty - BotB * 1e-10
  elseif extrapolation[1] == Constant
    BotB = 0.0
    BotC = Boty
  end # Note for the curve case we will simply not append the new block.
  BotRow = [Botx-1e-10, 0.0, BotB, BotC]

  # Now doing the extrapolation to the right.
  if extrapolation[2] == Linear # This is a bit more complicated than before because the
                                # coefficients by themselves give the gradients at the left
                                # of the interval. Here we want the gradient at the right.
    last_interval_width = Topx - fullMatrix[dim , 1]
    @inbounds grad_at_right = 2 * fullMatrix[dim , 2] * last_interval_width + fullMatrix[dim , 3]
    TopB = grad_at_right
    TopC = Topy
  elseif extrapolation[2] == Constant
    TopB = 0.0
    TopC = Topy
  elseif extrapolation[2] == Curve # We are just going to add this one on regardless as otherwise end of data interval information is lost.
    Gap  =  Topx - fullMatrix[dim ,1]
    TopA = fullMatrix[dim ,2]
    @inbounds TopB = fullMatrix[dim ,3] + 2 * TopA * Gap
    @inbounds TopC = fullMatrix[dim ,4] + TopB * Gap - TopA * Gap*Gap
  end
  TopRow = [ Topx, TopA ,TopB ,TopC]
  # Appending blocks and returning.
  if extrapolation[1] != Curve fullMatrix = vcat(BotRow' , fullMatrix) end
  fullMatrix = vcat(fullMatrix,  TopRow')
  return fullMatrix
end
