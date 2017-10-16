using HypothesisTests, Base.Test
const HT = HypothesisTests

srand(1234)

#test1: very different variances
σ = [0.1, 2.0, 1.5]
group = rand(1:3, 1000)
x = σ[group].*randn(1000)
F = HT.BrownForsyth(x, group)
@test  pvalue(F) ≈ 0.0

#test2: similar variances
σ = [2.1, 2.0, 1.8]
group = rand(1:3, 1000)
x = σ[group].*randn(1000)
F = HT.BrownForsyth(x, group)
@test pvalue(F) ≈ 0.006959798508421544

#test3: very similar variances
σ = [2.1, 2.0, 2.05]
group = rand(1:3, 1000)
F = HT.BrownForsyth(x, group)
@test pvalue(F) ≈ 0.5715114374730406
