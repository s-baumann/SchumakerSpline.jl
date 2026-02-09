import Base.+, Base.-, Base./, Base.*

function +(spl::Schumaker, num::Real)
    new_coefficients_ = copy(spl.coefficient_matrix_)
    @views new_coefficients_[:, 3] .+= num
    return Schumaker(spl.IntStarts_, new_coefficients_)
end
function -(spl::Schumaker, num::Real)
    new_coefficients_ = copy(spl.coefficient_matrix_)
    @views new_coefficients_[:, 3] .-= num
    return Schumaker(spl.IntStarts_, new_coefficients_)
end
function *(spl::Schumaker, num::Real)
    new_coefficients_ = spl.coefficient_matrix_ .* num
    return Schumaker(spl.IntStarts_, new_coefficients_)
end
function /(spl::Schumaker, num::Real)
    new_coefficients_ = spl.coefficient_matrix_ ./ num
    return Schumaker(spl.IntStarts_, new_coefficients_)
end
function +(num::Real, spl::Schumaker)
    return +(spl, num)
end
function -(num::Real, spl::Schumaker)
    new_coefficients_ = spl.coefficient_matrix_ .* -1
    @views new_coefficients_[:, 3] .+= num
    return Schumaker(spl.IntStarts_, new_coefficients_)
end
function *(num::Real, spl::Schumaker)
    return *(spl, num)
end
