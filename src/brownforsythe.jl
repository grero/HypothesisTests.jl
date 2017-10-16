struct BrownForsytheTest <: HypothesisTests.HypothesisTest
    F::Float64      #test statistics
    ng::Int64       #number of groups
    N::Int64        #number of observations
end

testname(BrownForsytheTest) = "Brown-Forsythe test for equality of group variances"
default_tail(x::BrownForsytheTest) = :right

function show_params(io::IO, x::BrownForsytheTest, indent)
    println(io, indent, "F-test statistic: ", x.F)
    println(io, indent, "Number of groups: ", x.ng)
    println(io, indent, "Number of observations: ", x.N)
end


function pvalue(x::BrownForsytheTest; tail=:right)
    pv = NaN
    if tail == :right
        pv = 1-cdf(FDist(x.ng, x.N-x.ng+1), x.F)
    end
    return pv
end

"""
Performs the Brown-Forsythe test that all groups in `x` have the same variance
"""
function BrownForsythe(x, group)
    N = length(x)
    groups = unique(group)
    sort!(groups)
    ngroups = length(groups)
    medians = zeros(ngroups)
    ng = zeros(Int64,ngroups)
    for i in 1:ngroups
        idx = group.==groups[i]
        medians[i] = median(x[idx])
        ng[i] += sum(idx)
    end
    z = zeros(length(x))
    zm = 0.0
    zgm = zeros(ngroups)
    for i in 1:length(x)
        gi = group[i]
        z[i] = abs(x[i] - medians[gi])
        zm += z[i]
        zgm[gi] += z[i]
    end
    zgm ./= ng
    zm /= N
    s = 0.0
    for i in 1:N
        gi = group[i]
        ss = z[i] - zgm[gi]
        s += ss*ss
    end
    a = 0.0
    for i in 1:ngroups
        aa = zgm[i] - zm
        a += ng[i]*aa*aa
    end
    F1 = (N-ngroups)*a
    F2 = (ngroups-1)*s
    F = F1/F2
    BrownForsytheTest(F, ngroups, N)
end
