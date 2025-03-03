using LinearAlgebra
using Plots

#################################################################
# Looking at Points
#################################################################

# @testset "constructors" begin 

#     function get_A_from_abg(p::EPoint) 
#         MB = expV(p.beta/2)
#         MAG = [exp(p.alpha) 0.0; 0.0 exp(p.gamma)] 
#         return MB * MAG * MB 
#     end

#     function get_A_from_computations(p::EPoint) 
#         MB = [p.coshfrac p.sinhfrac; p.sinhfrac p.coshfrac]
#         MAG = [p.expalpha 0.0; 0.0 p.expgamma] 
#         return MB * MAG * MB 
#     end

#     # Matrix representation

#     AA = [
#         [1.0 0.0; 0.0 1.0], 
#         [2.0 0.0; 0.0 1.0], 
#         [1.0 0.0; 0.0 2.0], 
#         [2.0 1.0; 1.0 2.0], 
#         [2.0 0.5; 0.5 0.5], 
#         [2 1; 1 1], 
#     ]

#     for A in AA
#         EA = EPoint(A)

#         # two forms of getting A
#         @test A ≈ [EA.c11 EA.c12; EA.c12 EA.c22]  
#         @test A ≈ get_A_from_abg(EA)
#         @test A ≈ get_A_from_computations(EA)

#         # other constructor 
#         EB = EPoint([A[1,1], A[1,2], A[2,2]])
#         @test EA == EB
#     end

#     # Alpha, beta, gamma representation 

#     ABG = [
#         [0.0, 0.0, 0.0],
#         [1.0, 0.0, 0.0], 
#         [0.0, 1.0, 0.0], 
#         [0.0, 0.0, 1.0], 
#         [1.0, 1.0, 0.0], 
#         [1.0, 0.0, 1.0], 
#         [0.0, 1.0, 1.0], 
#         [1.0, 1.0, 1.0], 
#         [0.5, 3.0, -2.0],
#         5.0 * (rand(3) .- 0.5)
#     ]

#     for abg in ABG 
#         EA = EPoint(abg[1], abg[2], abg[3])

#         MB = expV(EA.beta/2)
#         MAG = [exp(EA.alpha) 0.0; 0.0 exp(EA.gamma)] 
#         BB = MB * MAG * MB 

#         MBtilda = [EA.coshfrac EA.sinhfrac; EA.sinhfrac EA.coshfrac]
#         MAGtilda = [EA.expalpha 0.0; 0.0 EA.expgamma] 
#         BBtilda = MBtilda * MAGtilda * MBtilda 

#         @test abg ≈ [EA.alpha, EA.beta, EA.gamma]
#         @test BB ≈ [EA.c11 EA.c12; EA.c12 EA.c22]
#         @test BBtilda ≈ BB
    
#         # Another representation 
#         ch = cosh(abg[2]/sqrt(2))
#         sh = sinh(abg[2]/sqrt(2))
#         th = tanh((abg[1]-abg[3])/2)
#         CC = exp((abg[1]+abg[3])/2) * cosh((abg[1]-abg[3])/2) * [ch + th sh; sh ch - th]
#         @test CC ≈ BB
    
#     end

# end

#################################################################
# Checking computations 
#################################################################

# @testset "computations" begin

#     function get_params(a::Float64, b::Float64, g::Float64)
#         ch = cosh(b/sqrt(2))
#         sh = sinh(b/sqrt(2))
#         tau = tanh((a-g)/2)
#         rho = exp((a+g)/2) * cosh((a-g)/2)
#         AA = rho * [ch + tau sh; sh ch - tau]
#         return ch, sh, tau, rho, AA
#     end

#     # a0, b0, g0 = 5.0 * (rand(3) .- 0.5)
#     # a1, b1, g1 = 5.0 * (rand(3) .- 0.5)
#     # ch0, sh0, tau0, rho0, AA0 = get_params(a0, b0, g0)
#     # ch1, sh1, tau1, rho1, AA1 = get_params(a1, b1, g1)

#     # @test det(AA0) ≈ rho0^2 * (1.0 - tau0^2)
#     # @test AA0[1,1] + AA0[2,2] ≈ rho0 * 2.0 * ch0 

#     # invsqrtAA0 = inv(sqrt(AA0))
#     # midpoint = invsqrtAA0 * AA1 * invsqrtAA0
#     # trace = 2.0 * rho1 / (rho0 * (1-tau0^2)) * (cosh((b0-b1)/sqrt(2)) - tau0 * tau1)
#     # midet = rho1^2 * (1.0 - tau1^2) / (rho0^2 * (1.0 - tau0^2))
#     # @test midpoint[1,1] + midpoint[2,2] ≈ trace 
#     # @test det(midpoint) ≈ midet

#     # eigs = eigen(midpoint)
#     # lp = 0.5*(trace + sqrt(trace^2 - 4.0*midet))
#     # lm = 0.5*(trace - sqrt(trace^2 - 4.0*midet))
#     # @test sort(eigs.values) ≈ [lm, lp]
#     # @test get_distance(Ell, EPoint(AA0), EPoint(AA1)) ≈ sqrt(log(lm)^2 + log(lp)^2)

#     function getlog(B::Matrix{Float64})
#         traceB = B[1,1] + B[2,2]
#         deterB = det(B) 

#         # the eigenvalues of B are (traceB +- sqrt(traceB^2 - 4 * deterB))/2
#         # computing the coefficients kI and kB of the decomposition log(B) = kI * I + kB * B 
#         if abs(traceB^2 - 4 * deterB) <= 1e-8 # if both eigenvalues are "the same"
#             # then the unique eigenvalue is traceB/2 
#             kI = log(traceB/2.0) - 1.0
#             kB = 2.0 / traceB
#         else 
#             lBp = (traceB + sqrt(traceB^2 - 4 * deterB)) / 2.0
#             lBm = (traceB - sqrt(traceB^2 - 4 * deterB)) / 2.0
#             kI = (log(lBm) * lBp - lBm * log(lBp)) / (lBp - lBm)
#             kB = (log(lBp) - log(lBm)) / (lBp - lBm)
#         end

#         return [kI 0; 0 kI] .+ kB .* B
#     end

#     function getexp(B::Matrix{Float64})
#         trace = (B[1,1] + B[2,2]) / 2
#         deter = trace^2 - det(B)

#         return exp(trace) .* (cosh(sqrt(deter)) .* I(2) .+ sinh(sqrt(deter)) / sqrt(deter) .* (B .- trace * I(2)))
#     end

#     for (a,b,g) in [rand(3) for _ in 1:10]
#         AA = EPoint(a,b,g)

#         log0A = log([AA.c11 AA.c12; AA.c12 AA.c22])
#         log1A = getlog([AA.c11 AA.c12; AA.c12 AA.c22])
#         @test log0A ≈ log1A

#         exp0A = exp([AA.c11 AA.c12; AA.c12 AA.c22])
#         exp1A = getexp([AA.c11 AA.c12; AA.c12 AA.c22])
#         @test exp0A ≈ exp1A
#     end

# end

#################################################################
# Looking at distances 
#################################################################

# function get_dp()
#     V = rand(2,2)
#     V[1,1] = abs(V[1,1]) + 0.1
#     # det = V11 * V22 - V12 * V21 
#     V[2,2] = (0.001 + 0.2*rand() + V[1,2] * V[2,1]) / V[1,1]
#     return V
# end

# function get_sdp()
#     V = get_dp()
#     return 0.5 * (V + V')
# end

# function theD(AA::Matrix{Float64}, BB::Matrix{Float64})
#     return get_distance(Ell, EPoint(AA), EPoint(BB))
# end

# Q = get_sdp() # def pos mais pas sym 
# V = get_sdp()

#################################################################
# Looking at geodesics
#################################################################

# @testset "geodesics" begin 

#     oldtime = 0.0
#     newtime = 0.0
#     Nrand = 100

#     for (a0,b0,g0,a1,b1,g1,h) in [rand(7) for _ in 1:Nrand]
#         for param in [a0,b0,g0,a1,b1,g1]
#             param = 5.0 * (param - 0.5)
#         end
#         AA0 = EPoint(a0,b0,g0)
#         AA1 = EPoint(a1,b1,g1)
#         V01 = Velocity(AA1, 1.0, 1.0)

#         oldtime += @elapsed expold = get_exponential_slow(Ell, V01, AA0, h)
#         newtime += @elapsed expnew = get_exponential(Ell, V01, AA0, h)

#         @test expold.alpha ≈ expnew.alpha
#         @test expold.beta  ≈ expnew.beta 
#         @test expold.gamma ≈ expnew.gamma
#     end

#     println("Moyenne old : ", oldtime / Nrand, ", new : ", newtime / Nrand)

#     # alpha0 = 0.0
#     # alpha1 = 3.0
#     # beta0 = 3.0
#     # beta1 = 1.0
#     # gamma0 = 0.0
#     # gamma1 = -1.0

#     # AA = EPoint(alpha0, beta0, gamma0) # starting point 
#     # BB = EPoint(alpha1, beta1, gamma1) # target 

#     # vv = Velocity(BB, 1.0, 1.0)
#     # hh = LinRange(0.0, 1.0, 10)

#     # oldex = [get_exponential_slow(Ell, vv, AA, h) for h in hh]
#     # expis = [get_exponential(Ell, vv, AA, h) for h in hh]

#     # # p = plot_EPoints(Ell, oldex)
#     # # q = plot_EPoints(Ell, expis)
#     # # display(plot(p,q, layout=[1,1]))

# end

# tt = LinRange(0.0, 10.0, 200)
# doplot(M::Matrix{Float64}, W::Matrix{Float64}) = plot(tt, [theD(M, exp(t*W') * M * exp(t*W)) for t in tt]);
# doplot!(p, M::Matrix{Float64}, W::Matrix{Float64}) = plot!(p, tt, [theD(M, exp(t*W') * M * exp(t*W)) for t in tt]);

# AA = EPoint([2.0, 0.0, 1.0]) # starting point 
# BB = EPoint([1.0, 1.0, 2.0]) # target 

# vv = Velocity(BB, 1.0, 1.0)
# expis = [get_exponential(Ell, vv, AA, h) for h in LinRange(0.0,1.0,30)]
# p = plot_EPoints(Ell, expis)
# eigis = [eigen([expi.c11 expi.c12; expi.c12 expi.c22]) for expi in expis]
# x1is = [ee.vectors[1,1] for ee in eigis]
# x2is = [ee.vectors[2,1] for ee in eigis]

# q = scatter(x1is, label="x1", markersize=2, legend=:outerbottom)
# scatter!(q, x2is, label="x2", markersize=2)

# display(plot(p,q,layout=[1,1]))

# alpha0 = 0.0
# alpha1 = 3.0
# beta0 = 3.0
# beta1 = 1.0
# gamma0 = 0.0
# gamma1 = -1.0

# AA = EPoint(alpha0, beta0, gamma0) # starting point 
# BB = EPoint(alpha1, beta1, gamma1) # target 

# hh = LinRange(0.0,1.0,30)
# vv = Velocity(BB, 1.0, 1.0)
# expis = [get_exponential(Ell, vv, AA, h) for h in hh]
# alphas = [expi.alpha for expi in expis]
# betas  = [expi.beta  for expi in expis]
# gammas = [expi.gamma for expi in expis]
# sumags = [aa + gg for (aa,gg) in zip(alphas, gammas)]
# eigies = [eigen([expi.c11 expi.c12; expi.c12 expi.c22]) for expi in expis]
# eigiep = [maximum(eigy.values) for eigy in eigies]
# eigiem = [minimum(eigy.values) for eigy in eigies]
# traces = [expi.c11 + expi.c22 for expi in expis]
# prods = [([1.0 -1.0] * [expi.c11 expi.c12; expi.c12 expi.c22] * [1.0; -1.0])[1] for expi in expis]

# poly = fit(hh, betas, 6)
# println("poly : ", poly)

# p = plot(hh, betas)
# plot!(p, hh, [beta0 + h * (beta1 - beta0) for h in hh], ls=:dash)

# p = plot(hh, traces)
# plot!(p, hh, [exp((1-h) * AA.alpha + h * BB.alpha) + exp((1-h) * AA.gamma + h * BB.gamma) for h in hh], ls=:dash)
# display(p)

# p = plot(hh, prods)
# plot!(hh, traces)
# display(p)

# p = plot(hh, alphas)
# plot!(hh, gammas)
# plot!(hh, sumags)
# plot!(hh, [(1-h) * (alpha0 + gamma0) + h * (alpha1 + gamma1) for h in hh])
# display(p)

# p = plot(hh, eigiem)
# plot!(hh, eigiep)
# display(p)

# q = plot_EPoints(Ell, expis)
# display(plot(p,q, layout=[1,1]))

#################################################################
# Searching for an intersection of a geodesic
#################################################################

# @testset "intersection" begin

#     function get_params(a::Float64, b::Float64, g::Float64)
#         ch = cosh(b/sqrt(2))
#         sh = sinh(b/sqrt(2))
#         Tau = tanh((a-g)/2)
#         Rho = exp((a+g)/2) * cosh((a-g)/2)
#         AA = Rho * [ch + Tau sh; sh ch - Tau]
#         return ch, sh, Tau, Rho, AA
#     end

#     function rho(at::Float64, gt::Float64)
#         return exp(0.5*(at+gt)) * cosh(0.5*(at-gt))
#     end

#     function rhoprime(at::Float64, gt::Float64, atprime::Float64, gtprime::Float64)
#         meant  = 0.5 * (at + gt)
#         difft  = 0.5 * (at - gt)
#         meantp = 0.5 * (atprime + gtprime)
#         difftp = 0.5 * (atprime - gtprime)
#         return meantp * exp(meant) * cosh(difft) + difftp * exp(meant) * sinh(difft)
#     end

#     function tau(at::Float64, gt::Float64)
#         return tanh(0.5*(at-gt))
#     end

#     function tauprime(at::Float64, gt::Float64, atprime::Float64, gtprime::Float64)
#         difft  = 0.5 * (at - gt)
#         difftp = 0.5 * (atprime - gtprime)
#         return difftp / cosh(difft)^2
#     end

#     function lambdapm(sign::Float64, beta0::Float64, olbeta::Float64, at::Float64, gt::Float64, rho0::Float64, tau0::Float64)
#         coshfrac = cosh((beta0-olbeta)/sqrt(2))
#         taut = tau(at,gt)
#         return rho(at,gt) / (rho0 * (1.0-tau0^2)) * (coshfrac - tau0 * taut + sign * sqrt((coshfrac - tau0 * taut)^2 - (1.0 - tau0^2)*(1.0 - taut^2)))
#     end

#     function lambdapmprime(sign::Float64, beta0::Float64, olbeta::Float64, at::Float64, gt::Float64, atprime::Float64, gtprime::Float64, rho0::Float64, tau0::Float64)
#         coshfrac = cosh((beta0-olbeta)/sqrt(2))
#         taut = tau(at,gt)
#         tautp = tauprime(at,gt,atprime,gtprime)
#         bigsqrt = sqrt((coshfrac - tau0 * taut)^2 - (1.0-tau0^2)*(1.0-taut^2))
#         return rhoprime(at,gt,atprime,gtprime) / (rho0 * (1.0-tau0^2)) * 
#                 (coshfrac - tau0 * taut + sign * bigsqrt) + rho(at,gt) / (rho0 * (1.0-tau0^2)) * 
#                 (-tau0 * tautp + sign * ((coshfrac - tau0 * taut) * (-tau0 * tautp) + (1.0 - tau0^2)*taut*tautp) / bigsqrt)
#     end

#     function Phi(t::Float64, dist01::Float64, olbeta::Float64, a0::Float64, beta0::Float64, g0::Float64, a1::Float64, g1::Float64)
#         return t^2 * dist01^2 - log(lambdapm( 1.0, beta0, olbeta, a0+t*(a1-a0), g0+t*(g1-g0), rho(a0,g0), tau(a0,g0)))^2 - 
#                                 log(lambdapm(-1.0, beta0, olbeta, a0+t*(a1-a0), g0+t*(g1-g0), rho(a0,g0), tau(a0,g0)))^2 
#     end

#     function Phiprime(t::Float64, dist01::Float64, olbeta::Float64, a0::Float64, beta0::Float64, g0::Float64, a1::Float64, g1::Float64)
#         lp = lambdapm( 1.0, beta0, olbeta, a0+t*(a1-a0), g0+t*(g1-g0), rho(a0,g0), tau(a0,g0))
#         lm = lambdapm(-1.0, beta0, olbeta, a0+t*(a1-a0), g0+t*(g1-g0), rho(a0,g0), tau(a0,g0))
#         lpprime = lambdapmprime( 1.0, beta0, olbeta, a0+t*(a1-a0), g0+t*(g1-g0), a1-a0, g1-g0, rho(a0,g0), tau(a0,g0))
#         lmprime = lambdapmprime(-1.0, beta0, olbeta, a0+t*(a1-a0), g0+t*(g1-g0), a1-a0, g1-g0, rho(a0,g0), tau(a0,g0))
#         return 2.0 * (t * dist01^2 - lpprime * log(lp) / lp - lmprime * log(lm) / lm)
#     end

#     a0, b0, g0 = 5.0 * (rand(3) .- 0.5)
#     a1, b1, g1 = 5.0 * (rand(3) .- 0.5)
#     if abs(b1 - b0) <= 1e-6
#         b1 = b0 + 1.0
#     end
#     ch0, sh0, tau0, rho0, AA0 = get_params(a0, b0, g0)
#     ch1, sh1, tau1, rho1, AA1 = get_params(a1, b1, g1)
#     dist01 = get_distance(Ell, EPoint(AA0), EPoint(AA1))

#     vv = Velocity(EPoint(AA1), 1.0, 1.0)
#     tt = LinRange(0.1, 0.9, 5) # the true times 
#     geod = [get_exponential(Ell, vv, EPoint(AA0), t) for t in tt] # true geodesic
#     betas = [AAt.beta for AAt in geod] # true betas 

#     # game : recover tt from beta 
#     for (t, AAt, olbeta) in zip(tt, geod, betas)

#         # test
#         lp = lambdapm( 1.0, beta0, olbeta, a0+t*(a1-a0), g0+t*(g1-g0), rho(a0,g0), tau(a0,g0))
#         lm = lambdapm(-1.0, beta0, olbeta, a0+t*(a1-a0), g0+t*(g1-g0), rho(a0,g0), tau(a0,g0))
#         @test a0+t*(a1-a0) ≈ AAt.alpha
#         @test g0+t*(g1-g0) ≈ AAt.gamma
#         # @test t * dist01 ≈ sqrt(log(lp)^2 + log(lm)^2) 
#         # @test t * dist01 ≈ get_distance(Ell, EPoint(AA0), AAt)

#         tn = (olbeta - b0) / (b1 - b0) # initial point of the Newton method 
#         print("Erreur init ", abs(tn - t))
#         ceinture = 0; bretelles = 100; goon = true
#         while goon && (ceinture <= bretelles)
#             # tnp1 = tn - Phi(tn, dist01, olbeta, a0, b0, g0, a1, g1) / Phiprime(tn, dist01, olbeta, a0, b0, g0, a1, g1)
#             tnp1 = tn - 2.0 * Phi(tn, dist01, olbeta, a0, b0, g0, a1, g1) * Phiprime(tn, dist01, olbeta, a0, b0, g0, a1, g1)
#             if abs(tnp1 - tn) <= 1e-8
#                 goon = false
#             end
#             tn = tnp1
#             ceinture += 1
#         end
#         println(", $ceinture it pour ", abs(tn - t))

#     end

# end

#################################################################
# Looking at distances again
#################################################################

# # MA = Matrix([cosh(1) sinh(1); sinh(1) cosh(1)])
# # AA = EPoint(MA)
# # AAs = [EPoint(Matrix(MA^t)) for t in LinRange(-2.0, 2.0, 11)]

# fA(a::Float64, b::Float64, c::Float64) = EPoint(exp([a b/2; b/2 c]))
# a = 1.0
# c = 2.0
# AAs = [fA(a, b, c) for b in LinRange(-1.0, 1.0, 10)]
# # vv = Velocity(EPoint(MA^2), 1.0, 1.0)
# # BB = EPoint(MA^(-2))
# # AAs = [get_exponential(Ell, vv, BB, t) for t in LinRange(0.0, 1.0, 11)]
# p = plot_EPoints(Ell, AAs)
# display(p)

# N = 30
# sqlen = 2.0
# # betas = [0.0, 0.25, 0.5, 0.75, 1.0, 3.0]
# # betas = LinRange(0.0,1.0,10)
# betas = [0.0,1.0]
# # cols = [:red, :orange, :green, :brown, :blue, :black]
# cols = distinguishable_colors(length(betas), [RGB(1,0,0), RGB(0,0,1)])

# aas = LinRange(-sqlen, sqlen, N) # x
# ggs = LinRange(-sqlen, sqlen, N) # y

# plane_x = [aa for aa in aas for gg in ggs]
# plane_y = [gg for aa in aas for gg in ggs]

# p = plot(legend=false, xlabel="a", ylabel="g", camera=(-45,0))

# lambdap(a,b1,b2,g) = cosh(0.5*(a-g))^2 * cosh((b1-b2)/sqrt(2)) - sinh(0.5*(a-g))^2 + sqrt((cosh((b1-b2)/sqrt(2)) * cosh(0.5*(a-g))^2 - sinh(0.5*(a-g))^2)^2 - 1.0)
# lambdam(a,b1,b2,g) = cosh(0.5*(a-g))^2 * cosh((b1-b2)/sqrt(2)) - sinh(0.5*(a-g))^2 - sqrt((cosh((b1-b2)/sqrt(2)) * cosh(0.5*(a-g))^2 - sinh(0.5*(a-g))^2)^2 - 1.0)

# thetry(a,b1,b2,g) = sqrt(log(lambdap(a,b1,b2,g))^2 + log(lambdam(a,b1,b2,g))^2)

# for (beta,col) in zip(betas,cols)
#     plane_z = [get_distance(Ell, EPoint(aa,0.0,gg), EPoint(aa,beta,gg)) for aa in aas for gg in ggs]
#     # tryin_z = [cosh(ff*(aa-gg))*get_distance(Ell, EPoint(-3.,0.0,3.), EPoint(-3.,beta,3.))/cosh(-6.0*ff) for aa in aas for gg in ggs]
#     scatter3d!(p, plane_x, plane_y, plane_z, markersize=1, msc=col, color=col)
#     if beta > 0
#         tryin_z = [thetry(aa,beta,0.0,gg) for aa in aas for gg in ggs]
#         scatter3d!(p, plane_x, plane_y, tryin_z, markersize=1, msc=:brown, color=:brown)
#     end
# end

# display(p)

#################################################################
# Testing the mesh
#################################################################

# function ibetaFromImesh(imesh::Int64, nag::Int64)
#     return 1 + floor(Int, (imesh - 1) / nag^2)
# end

# function ialphaFromImesh(imesh::Int64, nag::Int64)
#     return 1 + floor(Int, (imesh - floor(Int,(imesh - 1) / nag^2)*nag^2 - 1) / nag)
# end

# function igammaFromImesh(imesh::Int64, nag::Int64)
#     return imesh - floor(Int,(imesh - 1) / nag^2)*nag^2 - floor(Int, (imesh - floor(Int,(imesh - 1) / nag^2)*nag^2 - 1) / nag) * nag
# end

# @testset "mesh_3D" begin 

#     domain = Ellipses("test", 1.0, 1.0, false) # do not change (1.0, 1.0)
    
#     # step = 0.7 # do not change (it is 0.7) 
#     step = 0.5 # do not change (it is 0.7) 
#     mesh = get_mesh(domain, step, true)

#     # p = plot_EPoints(domain, mesh.points)
#     # display(p)

#     lp = 12.0
#     lm = 1.0
#     AA = [(lp+lm)/2 (lp-lm)/2; (lp-lm)/2 (lp+lm)/2]
#     EA = EPoint(AA)
#     println("EA : ", EA)
#     p = plot_EPoint(domain, EA)
#     display(p)

#     # function format(fl) 
#     #     if abs(fl) <= 1e-8 
#     #         return "0.0  " 
#     #     else 
#     #         return "$fl"[1:min(end,5)]
#     #     end
#     # end

#     # ip = 1
#     # for ib in 1:mesh.nb 
#     #     for ia in 1:mesh.nag
#     #         for ig in 1:mesh.nag
#     #             println("Point $ip : (b,a,g) = ", format(mesh.points[ip].beta), ", ", format(mesh.points[ip].alpha), ", ", format(mesh.points[ip].gamma))
#     #             ip += 1
#     #         end
#     #     end
#     # end

#     # p = plot_mesh(domain, mesh)
#     # display(p)

#     # randlooks = 1 .+ round.(Int, (mesh.npoints-1) * rand(5))
#     # for randlook in randlooks 
#     #     # target = mesh.points[randlook]
#     #     # indice = search_in_mesh(Ell, mesh, mesh.points[randlook])
#     #     # result = mesh.points[indice]
#     #     # println("Searching ", target.alpha, ", ", target.beta, ", ", target.gamma)
#     #     # println("Found     ", result.alpha, ", ", result.beta, ", ", result.gamma)
#     #     @test search_in_mesh(Ell, mesh, mesh.points[randlook]) == randlook
#     #     ibeta  = ibetaFromImesh(domain.ctedet, randlook, mesh.nag)
#     #     ialpha = ialphaFromImesh(domain.ctedet, randlook, mesh.nag)
#     #     igamma = igammaFromImesh(domain.ctedet, randlook, mesh.nag)
#     #     @test randlook == 1 + (ibeta-1)*mesh.nag^2 + (ialpha-1)*mesh.nag + (igamma-1)
#     #     hatalpha = -domain.alga + 2.0*domain.alga * (ialpha-1)/(mesh.nag-1)
#     #     hatbeta  = -domain.beta + 2.0*domain.beta * (ibeta -1)/(mesh.nb -1)
#     #     hatgamma = -domain.alga + 2.0*domain.alga * (igamma-1)/(mesh.nag-1)
#     #     # println("ibeta $ibeta for $(mesh.nb) betas, hat ", hatbeta, " for ", mesh.points[randlook].beta)
#     #     # println("ialpha $ialpha for $(mesh.nag) alphas, hat ", hatalpha, " for ", mesh.points[randlook].alpha)
#     #     # println("igamma $igamma for $(mesh.nag) gammas, hat ", hatgamma, " for ", mesh.points[randlook].gamma)
#     #     @test abs(mesh.points[randlook].alpha - hatalpha) <= 1e-8 
#     #     @test abs(mesh.points[randlook].beta  - hatbeta)  <= 1e-8 
#     #     @test abs(mesh.points[randlook].gamma - hatgamma) <= 1e-8 
#     # end

#     # seeds = [
#     #     [mesh.points[1]], # one-point

#     #     [mesh.points[1],mesh.points[2]], # two adjacent points, same strata
#     #     [mesh.points[1],mesh.points[9]], # two non-adjacent points, same strata
#     #     [mesh.points[1],mesh.points[8]], # two non-adjacent points, same strata, gathering a neighbor
#     #     [mesh.points[1],mesh.points[10]], # two vertically aligned points, adjacent stratae
#     #     [mesh.points[1],mesh.points[19]], # two vertically aligned points, non-adjacent stratae
#     #     [mesh.points[1],mesh.points[13]], # diagonal points, adjacent stratae
#     #     [mesh.points[1],mesh.points[18]], # diagonal points, adjacent stratae
#     #     [mesh.points[1],mesh.points[22]], # diagonal points, non-adjacent stratae
#     #     [mesh.points[1],mesh.points[31]], # diagonal points, non-adjacent stratae
#     #     [mesh.points[1],mesh.points[end]], # diagonal points, non-adjacent stratae

#     #     [mesh.points[1],mesh.points[2],mesh.points[4]], # three adjacent, same strata
#     #     [mesh.points[1],mesh.points[3],mesh.points[7]], # three non adjacent, same strata
#     #     [mesh.points[1],mesh.points[6],mesh.points[7]], # three non adjacent, same strata
#     #     [mesh.points[1],mesh.points[2],mesh.points[10]], # three adjacent, different strata
#     #     [mesh.points[1],mesh.points[2],mesh.points[19]], # three nonadjacent, different strata
#     #     [mesh.points[1],mesh.points[2],mesh.points[37]], # three nonadjacent, different strata

#     #     [mesh.points[1], mesh.points[7], mesh.points[3], mesh.points[37]], # nawak
#     #     [mesh.points[1], mesh.points[7], mesh.points[3], mesh.points[9], 
#     #     mesh.points[37], mesh.points[39], mesh.points[43], mesh.points[45]], # whole mesh
#     #     [mesh.points[2], mesh.points[3], mesh.points[4], mesh.points[7], mesh.points[9], 
#     #     mesh.points[10], mesh.points[37], mesh.points[39], mesh.points[43], mesh.points[45]], # whole mesh but 1
#     # ]

#     # goals = [
#     #     [1],
        
#     #     [1,2],
#     #     [1,5,9],
#     #     [1,4,5,8],
#     #     [1,10],
#     #     [1,10,19],
#     #     [1,13],
#     #     [1,18],
#     #     [1,10,22],
#     #     [1,10,22,31],
#     #     [1,10,23,32,45], 

#     #     [1,2,4],
#     #     [1,2,3,4,5,7],
#     #     [1,2,4,5,6,7,8],
#     #     [1,2,10],
#     #     [1,2,10,19],
#     #     [1,2,10,11,19,28,37],

#     #     [1,2,3,4,5,7,10,11,13,19,20,22,28,37],
#     #     collect(1:mesh.npoints),
#     #     collect(2:mesh.npoints)
#     # ]

#     # for (seed, goal) in zip(seeds, goals)
#     #     convex = sort(get_mesh_convex(domain, mesh, seed))
#     #     @test convex == goal

#     #     if seed == seeds[end]
#     #         println("convex : ", convex)
#     #         p = plot_convex_in_mesh(domain, mesh, convex)
#     #         display(p)

#     #         # vv = Velocity(mesh.points[40],1.0,1.0)
#     #         # geod = [get_exponential(domain, vv, mesh.points[3], h) for h in LinRange(0.0, 1.0, 40)]
#     #         # aas = [geo.alpha for geo in geod]
#     #         # ggs = [geo.gamma for geo in geod]
#     #         # bbs = [sign(geo.beta) * get_distance(domain, geo, EPoint(geo.alpha,0.0,geo.gamma)) for geo in geod]
#     #         # scatter3d!(aas,ggs,bbs,color=:brown, msc=:brown, ms=3)
#     #         # display(p)
#     #     end
#     # end

#     # step = 0.3
#     # mesh = get_mesh(domain, step, true)
#     # println("Mesh : ", mesh)

#     # convex = get_mesh_convex(domain, mesh, [mesh.points[1], mesh.points[end]])
#     # println("convex : ", convex)
#     # p = plot_convex_in_mesh(domain, mesh, convex)
#     # display(p)

# end

# @testset "mesh_2D" begin 

#     domain = Ellipses("test", 1.0, 1.0, true) # do not change (1.0, 1.0)
    
#     # step = 0.7 # do not change (it is 0.7) 
#     step = 1.5 # do not change (it is 0.7) 
#     mesh = get_mesh(domain, step, true)

#     # p = plot_mesh(domain, mesh)
#     # for ip in 1:mesh.npoints
#     #     ipa = FLagHada.ialphaFromImesh(domain.ctedet, ip, mesh.nag)
#     #     ipb = FLagHada.ibetaFromImesh(domain.ctedet, ip, mesh.nag)
#     #     for jp in 1:mesh.npoints
#     #         jpa = FLagHada.ialphaFromImesh(domain.ctedet, jp, mesh.nag)
#     #         jpb = FLagHada.ibetaFromImesh(domain.ctedet, jp, mesh.nag)
#     #         if ((ip != jp) 
#     #             && ((abs(ipa-3) == abs(jpa-3)) # symmetric stratae with respect to beta = 0
#     #             && (ipb == jpb))
#     #             || ((abs(ipb-3) == abs(jpb-3))
#     #             && (ipa == jpa))
#     #         )
#     #             vv = Velocity(mesh.points[ip],1.0)
#     #             dist = get_distance(domain, mesh.points[ip], mesh.points[jp])
#     #             geod = [get_exponential(domain, vv, mesh.points[jp], h*dist) for h in LinRange(0.0, 1.0, 40)]
#     #             aas = [geo.alpha for geo in geod]
#     #             bbs = [sign(geo.beta) * get_distance(domain, geo, EPoint(geo.alpha,0.0,-geo.alpha)) for geo in geod]
#     #             plot!(p,aas,bbs,color=:black)
#     #         end
#     #     end
#     # end
#     # display(p)

#     p = plot_EPoints(domain, mesh.points)
#     display(p)

#     # randlooks = 1 .+ round.(Int, (mesh.npoints-1) * rand(5))
#     # for randlook in randlooks 
#     #     # target = mesh.points[randlook]
#     #     # indice = search_in_mesh(Ell, mesh, mesh.points[randlook])
#     #     # result = mesh.points[indice]
#     #     # println("Searching ", target.alpha, ", ", target.beta, ", ", target.gamma)
#     #     # println("Found     ", result.alpha, ", ", result.beta, ", ", result.gamma)
#     #     @test search_in_mesh(Ell, mesh, mesh.points[randlook]) == randlook
#     #     ibeta  = FLagHada.ibetaFromImesh(domain.ctedet, randlook, mesh.nag)
#     #     ialpha = FLagHada.ialphaFromImesh(domain.ctedet, randlook, mesh.nag)
#     #     @test randlook == 1 + (ibeta-1)*mesh.nag + (ialpha-1)
#     #     hatalpha = -domain.alga + 2.0*domain.alga * (ialpha-1)/(mesh.nag-1)
#     #     hatbeta  = -domain.beta + 2.0*domain.beta * (ibeta -1)/(mesh.nb -1)
#     #     # println("ibeta $ibeta for $(mesh.nb) betas, hat ", hatbeta, " for ", mesh.points[randlook].beta)
#     #     # println("ialpha $ialpha for $(mesh.nag) alphas, hat ", hatalpha, " for ", mesh.points[randlook].alpha)
#     #     # println("igamma $igamma for $(mesh.nag) gammas, hat ", hatgamma, " for ", mesh.points[randlook].gamma)
#     #     @test abs(mesh.points[randlook].alpha - hatalpha) <= 1e-8 
#     #     @test abs(mesh.points[randlook].beta  - hatbeta)  <= 1e-8 
#     #     @test abs(mesh.points[randlook].gamma + hatalpha) <= 1e-8 
#     # end

#     # seeds = [
#     #     [mesh.points[1]], # one-point

#     #     [mesh.points[1],mesh.points[2]], # two adjacent points, same strata
#     #     [mesh.points[1],mesh.points[4]], # two non-adjacent points, same strata
#     #     [mesh.points[1],mesh.points[6]], # two vertically aligned points, adjacent stratae
#     #     [mesh.points[1],mesh.points[11]], # two vertically aligned points, non-adjacent stratae
#     #     [mesh.points[1],mesh.points[7]], # diagonal points, adjacent stratae
#     #     [mesh.points[1],mesh.points[12]], # diagonal points, non-adjacent stratae

#     #     [mesh.points[1],mesh.points[2],mesh.points[3]], # three adjacent, same strata
#     #     [mesh.points[1],mesh.points[2],mesh.points[7]], # three adjacent, different strata
#     #     [mesh.points[3],mesh.points[12],mesh.points[14]], # three nonadjacent, different strata

#     #     [mesh.points[1], mesh.points[5], mesh.points[11], mesh.points[21], mesh.points[25], mesh.points[15]], # whole mesh
#     #     [mesh.points[2], mesh.points[5], mesh.points[6], mesh.points[11], mesh.points[21], mesh.points[25], mesh.points[15]] # whole mesh but 1
#     # ]

#     # goals = [
#     #     [1],
        
#     #     [1,2],
#     #     [1,2,3,4],
#     #     [1,6], 
#     #     [1,6,11],
#     #     [1,7],
#     #     [1,7,12],

#     #     [1,2,3],
#     #     [1,2,7],
#     #     [3,8,12,13,14],

#     #     collect(1:mesh.npoints),
#     #     collect(2:mesh.npoints)
#     # ]

#     # for (seed, goal) in zip(seeds, goals)
#     #     convft = get_mesh_convex(domain, mesh, seed)
#     #     println("Convex : ", convft)
#     #     convex = sort(convft)
#     #     @test convex == goal

#     #     if seed == seeds[end]
#     #         println("convex : ", convex)
#     #         p = plot_convex_in_mesh(domain, mesh, convex)
#     #         display(p)

#     #         # vv = Velocity(mesh.points[40],1.0,1.0)
#     #         # geod = [get_exponential(domain, vv, mesh.points[3], h) for h in LinRange(0.0, 1.0, 40)]
#     #         # aas = [geo.alpha for geo in geod]
#     #         # ggs = [geo.gamma for geo in geod]
#     #         # bbs = [sign(geo.beta) * get_distance(domain, geo, EPoint(geo.alpha,0.0,geo.gamma)) for geo in geod]
#     #         # scatter3d!(aas,ggs,bbs,color=:brown, msc=:brown, ms=3)
#     #         # display(p)
#     #     end
#     # end
# end

## Testing the boundary 
mat(A::EPoint) = [A.c11 A.c12; A.c12 A.c22]

function prodscal(V1, V2, Q) 
    # isqrtQ = inv(sqrt(mat(Q)))
    # return sum(diag(isqrtQ * V1 * isqrtQ^2 * V2 * isqrtQ))
    return sum(diag(V1 * V2))
end

domain = Ellipses("test", 1.0, 1.0, false)
XX = EPoint(2.0 .* (rand(3) .- 0.5)...)
sqrtX = sqrt(mat(XX))
isqrtX = inv(sqrtX)
YY = EPoint(2.0 .* (rand(3) .- 0.5)...)
sqrtY = sqrt(mat(YY))
isqrtY = inv(sqrtY)
VV = rand(2,2)
VV .= (VV .+ VV') ./ 2.0

ZZ = exp(10 * VV)
VVx = log(isqrtX * ZZ * isqrtX) / get_distance(domain, XX, EPoint(ZZ))
VVy = log(isqrtY * ZZ * isqrtY) / get_distance(domain, YY, EPoint(ZZ))
vv = Velocity(EPoint(ZZ), 1.0, 1.0)

hh = LinRange(0.0, 1.0, 200)[2:end]
# candidate = - 2.0 * prodscal(VVx, log(isqrtX * mat(YY) * isqrtX), XX)
# quotients = [(get_distance(domain, get_exponential(domain, vv, XX, h), YY)^2 - dxy2)/h for h in hh]
# p = plot(hh, quotients)
# plot!(p, hh, candidate .* ones(length(hh)))
# display(p)

dxy2 = get_distance(domain, XX, YY)^2
derx = - 2.0 * prodscal(VVx, log(isqrtX * mat(YY) * isqrtX), XX)
dery = - 2.0 * prodscal(VVy, log(isqrtY * mat(XX) * isqrtY), YY)
d1 = [((get_distance(domain, get_exponential(domain, vv, XX, h), YY)^2 - dxy2)/h - derx)/h for h in hh]
d2 = [((get_distance(domain, get_exponential(domain, vv, YY, h), XX)^2 - dxy2)/h - dery)/h for h in hh]
p = plot(legend=false)#, xlims=[-0.1,1.1], ylims=[-0.1,0.1])
plot!(p, hh, d1)
plot!(p, hh, d2)
# plot!(p, hh, d1 .- d2)
display(p)