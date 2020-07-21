using DelimitedFiles
using LinearAlgebra
using Plots

kn_range = 0.1:0.1:1.0

V100 = readdlm("V_100.mat", ' ', Float64, '\n');
F100 = svd(V100);
F100SN = normalize(F100.S, Inf)

svArray = []

for R=10:10:100
  # low rank approximation with rank R
  # R = 50
  V100LR = F100.U[:,1:R] * diagm(F100.S[1:R]) * F100.V[:,1:R]'

  # sanity check: Frobenius norm should be S[R+1] is not the case but still low
  # ! FIXME
  println(norm(V100 - V100LR ,2))

  L = F100.U[:,1:R]
  OPCompl = I(100) - L * L'  # projector to complement of LL'

  # ansatz
  # p(t1, t2, deg) = t1^deg  # monomials
  p(t1, t2, deg) = (1/t1)^deg  # monomials
  maxdeg = 1
  P = vcat(
    [[p(kn1, kn2, d) for d=1:1]' for kn1 in kn_range, kn2 in kn_range]...
  )
  Q = Matrix(qr(P).Q)
  objective = OPCompl * Q  # --> MIN!!!

  push!(svArray, svd(objective).S...)
end

plot([svd(OPCompl).S, svArray], yaxis=:log, title="Singular Values, max ansatz degree $maxdeg, R=$R", xlabel="index", ylabel="value", label=["(I - LL')" "rank=10,20,...,100"], marker = (:dot))



if false # comparison
  V361 = readdlm("V_361.mat", ' ', Float64, '\n');
  F361 = svd(V361);
  F361SN = normalize(F361.S, Inf)

  V400 = readdlm("V_400.mat", ' ', Float64, '\n');
  F400 = svd(V400);
  F400SN = normalize(F400.S, Inf)

  max = 100
  plot([F100SN, F400SN[1:max]], yaxis=:log, title="Normlized Singular Values", xlabel="index", ylabel="value", label=["h=0.1" "h=0.05"], marker = (:dot))
end


