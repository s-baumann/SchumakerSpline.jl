using SchumakerSpline
from = 1
to = 1000
x1 = collect(range(from, stop=to, length=1000))
y1 = log.(x1)
x2 = vcat(x1[1:500] .* 0.01, x1[500]  .+ (0.05 .* x1[500:1000]))
y2 = sqrt.(x2)
s1 = Schumaker(x1,y2)
s2 = Schumaker(x2,y2)

function get_crossover_in_interval(s1::Schumaker{T}, s2::Schumaker{R}, interval::Tuple{G,G}) where T<:Real where R<:Real where G<:Real

end

function loop_schumakers_to_detect_crossover(s1::Schumaker{T}, s2::Schumaker{R}) where T<:Real where R<:Real
    # What x locations to loop over
    all_starts = sort(unique(vcat(s1.IntStarts_, s2.IntStarts_)))
    start_of_overlap = maximum([minimum(s1.IntStarts_), minimum(s2.IntStarts_)])
    overlap_starts = all_starts[all_starts .> start_of_overlap]
    # Getting a container to return results
    promo_type = promote_type(T,R)
    locations_of_crossovers = Array{promo_type,1}()
    # For the first part what function is higher.
    last_one_greater = evaluate(s1, overlap_starts[1]) > evaluate(s2, overlap_starts[1])
    for start in overlap_starts[2:length(overlap_starts)]
        val_1 = evaluate(s1, start)
        val_2 = evaluate(s2, start)
        # Need to take into account the ordering and record when it flips
        one_greater = val_1 > val_2
        if one_greater != last_one_greater
            indices_containing_crossover = (i, j)
            push!(locations_of_crossovers, indices_containing_crossover)
        end
        last_one_greater = one_greater
    end
    return indices_of_flips
end


locations_of_crossovers(s1,s2)

typeof((1,2))
