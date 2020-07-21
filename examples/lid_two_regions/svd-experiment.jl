using DelimitedFiles
using LinearAlgebra
using Plots
using SymPy

kn_range = 0.1:0.1:1.0

V100 = readdlm("V_100.mat", ' ', Float64, '\n');
F100 = svd(V100);
F100SN = normalize(F100.S, Inf)

# low rank approximation with rank R
R = 20
V100LR = F100.U[:,1:R] * diagm(F100.S[1:R]) * F100.V[:,1:R]'

# sanity check: Frobenius norm should be S[R+1] is not the case but still low
# ! FIXME
println(norm(V100 - V100LR ,2))

L = F100.U[:,1:R]
OPCompl = I(100) - L * L'  # projector to complement of LL'

# ansatz
# p(t1, t2, deg) = t1^deg  # monomials
# p(t1, t2, deg) = (1/t1)^deg  # 1/x
# p(t1, t2, deg) = exp(t1*deg)  # exponentials
# p(t1, t2, deg) = log(t1)^deg  # logrithms
# p(t1, t2, deg) = (t1*t2)^deg  # logrithms
# p(t1, t2, deg) = sin(deg * pi * ((t2-0.1)*100+t1*10) / 100)  # logrithms
# p(t1, t2, deg) = (t2-0.1)*100+t1*10  # logrithms
# p(t1, t2, deg) = 1/(t1+t2)
# p(t1, t2, deg) = t2 + 1/(t1)^deg
p(t1, t2, deg) = 1/(t2)^deg + 1/(t1)^deg
# function p(t1, t2, deg)
#   x = Sym("x")
#   f(x) = x^float(deg-1)
#   G(y) = subs(integrate(f(x), x), x, y)
#   return N(G(float(t1)))
# end

degrange = -3:3
P = vcat(
  [[p(kn1, kn2, d) for d in degrange]' for kn1 in kn_range, kn2 in kn_range]...,
)
Q = Matrix(qr(P).Q)
objective = OPCompl * Q  # --> MIN!!!

display(plot([svd(OPCompl).S, svd(objective).S], yaxis=:log, title="Singular Values, max ansatz degree, R=$R", xlabel="index", ylabel="value", label=["(I - LL')" "(I - LL') * Q"], marker = (:dot)))



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


