using StatsPlots

T = 350

#   パラメータ設定
α1, α2, α3, α4, α5, α6 = 0.9, 0.1, 0.1, 0.02, 0.1, 0.5
β1, β2 = 0.05, 1.0
γ1, γ2 = 0.015, 0.02
δ1, δ2 = fill(0.02, T), 0.3 #   0.2, 0.3
ϵ1, ϵ2, ϵ3, ϵ4, ϵ5, ϵ6, ϵ7 = 0.5, 1.0, 0.7, 0.05, 0.693, 0.7, 0.3
ζ1, ζ2, ζ3, ζ4, ζ5 = 1.0, fill(150.0, T), fill(150.0, T), 0.02, 1.0  #   1.0, 150.0, 150.0, 0.02, 1.0
for t = 2:T
    ζ2[t] = ζ2[t-1]*(1 + ζ4)    #   ここ元に戻す
    ζ3[t] = ζ3[t-1]*(1 + ζ4)
end
θ1, θ2 = 0.4, 0.1
ι1, ι2, ι3, ι4, ι5 = 0.1, 1.0, 2.0, 0.5, 5.0
κ1, κ2, κ3, κ4, κ5 = 0.5, 0.02, 0.02, 0.9, 0.5
λe = 0.5
μ1, μ2, μ3 = 0.3, 0.1, 0.5
τ1, τ2 = 0.2, 0.2
uT = 0.8
G0 = 100.0
r = 0.01
muc, muk = 0.3, 0.3

#   配列定義
Wc, Wk, Wb, Wg, W, wg = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
Tew, Tei, Tec, Tek, Te = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
p, pk, pec, pek = zeros(T), zeros(T), zeros(T), zeros(T)
G, GD, g, gD = zeros(T), zeros(T), zeros(T), zeros(T)
Cw, Ci, C, cw, ci, c, CwD, CiD = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
Kc, Kk, Kg, kc, kk, kg, K, k = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
Mw, Mi, Mc, Mk, M = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
Eci, Ec, Ecb, Eki, Ek, Ekb, eci, ec, ecb, eki, ek, ekb = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
Lw, Lc, Lk, L, LcD, LkD = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
Hw, Hb, H = zeros(T), zeros(T), zeros(T)
NWw, NWi, NWc, NWk, NWb, NWg = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
NLw, NLi, NLc, NLk, NLb, NLg = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
Tec, Tek, Tei, Tew, Te = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
Tii, Tiw, Ti = zeros(T), zeros(T), zeros(T)
Tfc, Tfk, Tfb, Tf = zeros(T), zeros(T), zeros(T), zeros(T)
Πci, Πc, Πcc, Πcb, Πki, Πk, Πkk, Πkb = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
Ic, Ik, If, Ig, ic, ik, ig, icD, ikD, igD = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
uc, uk, Δuc, Δuk, uce, uke = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
GS, rE = zeros(T), zeros(T)

Δkc, Δkk, Δkg, Δk, ΔKc, ΔKk, ΔKg, ΔK = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
ΔMw, ΔMi, ΔMc, ΔMk, ΔM = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
ΔLw, ΔLc, ΔLk, ΔL = zeros(T), zeros(T), zeros(T), zeros(T)
Δeci, Δecb, Δec, Δeki, Δekb, Δek, ΔEci, ΔEcb, ΔEc, ΔEki, ΔEkb, ΔEk = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
ΔHw, ΔHb, ΔH = zeros(T), zeros(T), zeros(T)
Δp, Δpk, Δpec, Δpek = zeros(T), zeros(T), zeros(T), zeros(T)
Δc, Δg, Δic, Δig = zeros(T), zeros(T), zeros(T), zeros(T)
ΔTfc, ΔTfk, ΔTiw, ΔTew, ΔTii, ΔTei = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
ΔW, ΔCw, ΔCi = zeros(T), zeros(T), zeros(T)
ΔNWi = zeros(T)

EciT, EkiT, EcbT, EkbT = zeros(T), zeros(T), zeros(T), zeros(T)
Tfce, Tfke, Tiwe, Tewe, Tiie, Teie = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
ce, ge, ice, ige = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
We = zeros(T)
Πcie, Πkie = zeros(T), zeros(T)
NWie = zeros(T)

#   初期値設定
Kc[1], Kk[1], Kg[1] = 150, 150, 150
Mw[1], Mi[1], Mc[1], Mk[1], M[1] = 400.0, 50.0, 75.0, 75.0, 600.0
Eci[1], Ec[1], Ecb[1], Eki[1], Ek[1], Ekb[1] = 50.0, 100.0, 50.0, 50.0, 100.0, 50.0
eci[1], ec[1], ecb[1], eki[1], ek[1], ekb[1] = 50.0, 100.0, 50.0, 50.0, 100.0, 50.0
Lw[1], Lc[1], Lk[1], L[1] = 200.0, 0.0, 0.0, 200.0
Hw[1], Hb[1], H[1] = 20.0, 580.0, 600.0
GD[1], gD[1], igD[1] = G0, G0, 0.5*G0
Wc[1], Wk[1], Wg[1], W[1], wg[1] = 0.5*G0, 0.5*G0, 0.5*G0, 1.5*G0, 0.5*G0
p[1], pk[1], pec[1], pek[1] = 1.0, 1.0, 1.0, 1.0
kc[1], kk[1], kg[1] = Kc[1]/pk[1], Kk[1]/pk[1], Kg[1]/pk[1]
K[1], k[1] = Kc[1] + Kk[1] + Kg[1], kc[1] + kk[1] + kg[1]
ic[1], ik[1], ig[1] = 0.15*G0, 0.15*G0, 0.5*G0
ice[1], ige[1] = ic[1], ig[1]
CwD[1], CiD[1] = 1.0*G0, 0.1*G0
Tfc[1], Tfk[1] = 5.0, 5.0
Cw[1], Ci[1] = CwD[1]/max(1, (CwD[1] + CiD[1] + GD[1])/(p[1]*ζ1*kc[1])), CiD[1]/max(1, (CwD[1] + CiD[1] + GD[1])/(p[1]*ζ1*kc[1]))
G[1] = GD[1]/max(1, (CwD[1] + CiD[1] + GD[1])/(p[1]*ζ1*kc[1]))
cw[1], ci[1], g[1] = Cw[1]/p[1], Ci[1]/p[1], G[1]/p[1]
Ic[1], Ik[1], If[1], Ig[1] = ic[1]*p[1], ik[1]*p[1], ic[1]*p[1] + ik[1]*p[1], ig[1]*p[1]
C[1], c[1] = Cw[1] + Ci[1], cw[1] + ci[1]

uc[1], uk[1] = min(1, (C[1] + G[1])/(p[1] * ζ1 * kc[1])), min(1, (Ic[1] + Ig[1])/(pk[1] * ζ1 * kk[1]))
uce[1], uke[1] = uc[1], uk[1]
Πc[1] = C[1] + G[1] - Wc[1] - Tec[1] - Tfc[1] - r*Lc[1]
Πci[1] = max(0, θ1*(Πc[1]-Ic[1])*Eci[1]/Ec[1] + θ2*(Mc[1]-Lc[1]))
Πcb[1] = max(0, θ1*(Πc[1]-Ic[1])*Ecb[1]/Ec[1] + θ2*(Mc[1]-Lc[1]))
Πcc[1] = Πc[1] - Πci[1] - Πcb[1]
Πk[1] = If[1] + Ig[1] - Wk[1] - Tek[1] - Tfk[1] - r*Lk[1]
Πki[1] = max(0, θ1*(Πk[1]-Ik[1])*Eki[1]/Ek[1] + θ2*(Mk[1]-Lk[1]))
Πkb[1] = max(0, θ1*(Πk[1]-Ik[1])*Ekb[1]/Ek[1] + θ2*(Mk[1]-Lk[1]))
Πkk[1] = Πk[1] - Πki[1] - Πkb[1]
Tiw[1], Tii[1], Tei[1], Tew[1] = τ1*W[1], τ1*(Πci[1] + Πki[1]), γ2*(Mi[1] + Eci[1] + Eki[1]), γ2*(Mw[1] + Hw[1])
Tfce[1], Tfke[1], Tiwe[1], Tewe[1], Tiie[1], Teie[1] = Tfc[1], Tfk[1], Tiw[1], Tew[1], Tii[1], Tei[1]
We[1] = W[1]
Πcie[1], Πkie[1] = Πci[1], Πki[1]
rE[1] = (Πci[1] + Πcb[1] + Πki[1] + Πkb[1])/(Ec[1] + Ek[1])

#   会計的整合性を満たすように、残りの変数の値を決める
#   初期値の計算においては、使わない会計恒等式も多い
NLw[1] = -Cw[1] + W[1] - Tiw[1] - Tew[1] - r*Lw[1]
ΔMw[1] = NLw[1] + ΔLw[1] - ΔHw[1]
NLi[1] = Πci[1] + Πki[1] - Ci[1] - Tii[1] - Tei[1]
ΔMi[1] = NLi[1] - pec[1]*Δeci[1] - pek[1]*Δeki[1]
NLc[1] = Πcc[1] - Ic[1]
NLk[1] = Πkk[1] - Ik[1]
ΔMc[1] = NLc[1] + ΔLc[1] + pec[1]*Δec[1]
ΔMk[1] = NLk[1] + ΔLk[1] + pek[1]*Δek[1]
Ti[1] = Tiw[1] + Tii[1]
Te[1] = Tew[1] + Tei[1] + Tec[1] + Tek[1]
Tf[1] = Tfc[1] + Tfk[1] + Tfb[1]
ΔM[1] = ΔMw[1] + ΔMi[1] + ΔMc[1] + ΔMk[1]
ΔL[1] = ΔLw[1] + ΔLc[1] + ΔLk[1]
NLb[1] = -Wb[1] - Tfb[1] + Πcb[1] + Πkb[1] + r*L[1]
ΔHb[1] = NLb[1] + ΔM[1] - ΔL[1] - pec[1]*Δecb[1] - pek[1]*Δekb[1]
ΔH[1] = ΔHw[1] + ΔHb[1]
GS[1] = -G[1] - Wg[1] + Ti[1] + Te[1] + Tf[1]
NLg[1] = -Ig[1] + GS[1]
K[1] = Kc[1] + Kk[1] + Kg[1]

NWw[1] = Mw[1] - Lw[1] + Hw[1]
NWi[1] = Mi[1] + Eci[1] + Eki[1]
NWc[1] = Kc[1] + Mc[1] - Ec[1] - Lc[1]
NWk[1] = Kk[1] + Mk[1] - Ek[1] - Lk[1]
NWb[1] = -M[1] + Ecb[1] + Ekb[1] + L[1] + Hb[1]
NWg[1] = Kg[1] - H[1]

NWie[1] = NWi[1]

#   会計的整合性の確認
println("------initial------")
println("-Ig+GS+ΔH=", -Ig[1]+GS[1]+ΔH[1])
println("NLw+NLi+NLc+NLk+NLb+NLg=", NLw[1]+NLi[1]+NLc[1]+NLk[1]+NLb[1]+NLg[1])
println("-M+Mw+Mi+Mc+Mk=", -M[1]+Mw[1]+Mi[1]+Mc[1]+Mk[1])
println("-Ec+Eci+Ecb=", -Ec[1]+Eci[1]+Ecb[1])
println("-Ek+Eki+Ekb=", -Ek[1]+Eki[1]+Ekb[1])
println("-L+Lw+Lc+Lk=", -L[1]+Lw[1]+Lc[1]+Lk[1])
println("-H+Hw+Hb=", -H[1]+Hw[1]+Hb[1])
println("NWw+NWi+NWc+NWk+NWb+NWg-K=", NWw[1]+NWi[1]+NWc[1]+NWk[1]+NWb[1]+NWg[1]-K[1])

function btw(A, B, C)
    if A > C
        return nothing
    else
        return min(max(A, B), C)
    end
end

function run()
    #   シミュレーション実行
    for t = 2:T
        Tfce[t] = ((1 - λe)*Tfce[t-1] + λe*Tfc[t-1])
        ce[t] = ((1 - λe)*ce[t-1] + λe*c[t-1])
        ge[t] = ((1 - λe)*ge[t-1] + λe*g[t-1])
        Tfke[t] = ((1 - λe)*Tfke[t-1] + λe*Tfk[t-1])
        ice[t] = ((1 - λe)*ice[t-1] + λe*ic[t-1])
        ige[t] = ((1 - λe)*ige[t-1] + λe*ig[t-1])
        We[t] = ((1 - λe)*We[t-1] + λe*W[t-1])
        Tiwe[t] = ((1 - λe)*Tiwe[t-1] + λe*Tiw[t-1])
        Tewe[t] = ((1 - λe)*Tewe[t-1] + λe*Tew[t-1])
        Πcie[t] = ((1 - λe)*Πcie[t-1] + λe*Πci[t-1])
        Πkie[t] = ((1 - λe)*Πkie[t-1] + λe*Πki[t-1])
        Tiie[t] = ((1 - λe)*Tiie[t-1] + λe*Tii[t-1])
        Teie[t] = ((1 - λe)*Teie[t-1] + λe*Tei[t-1])
        NWie[t] = ((1 - λe)*NWie[t-1] + λe*NWi[t-1])

        Ωc, Ωk = ϵ5*p[t-1]*(c[t-1] + g[t-1]), ϵ5*p[t-1]*(ic[t-1] + ig[t-1])
        Wc[t] = (1 - ϵ7)*Wc[t-1] + ϵ7*((1 - ϵ6)*Ωc + ϵ6*Wc[t-1]*exp(ϵ2*((CwD[t-1] + CiD[t-1] + GD[t-1])/(p[t-1]*min(ζ1*(kc[t-1]-Δkc[t-1]), ζ2[t-1])) - uT)))
        Wk[t] = (1 - ϵ7)*Wk[t-1] + ϵ7*((1 - ϵ6)*Ωk + ϵ6*Wk[t-1]*exp(ϵ2*((icD[t-1] + ikD[t-1] + igD[t-1])/(min(ζ5*(kk[t-1]-Δkk[t-1]), ζ3[t-1])) - uT)))
        Wb[t] = (1 - ϵ1)*Wb[t-1] + ϵ1*max(0, ϵ3*(Πcb[t-1] + Πkb[t-1] + r*L[t-1]) + ϵ4*(Hb[t-1] - M[t-1]))
        Tec[t] = γ1*Kc[t-1]
        Tek[t] = γ1*Kk[t-1]
        Tei[t] = γ2*(Mi[t-1] + Eci[t-1] + Eki[t-1])
        ΔTei[t] = Tei[t] - Tei[t-1]
        Tew[t] = γ2*(Mw[t-1] + Hw[t-1])
        ΔTew[t] = Tew[t] - Tew[t-1]
        Te[t] = Tew[t] + Tei[t] + Tec[t] + Tek[t] 
        p[t] = (1 + muc)*(pk[t]*ice[t] + Wc[t] + Tec[t] + Tfce[t] + r*Lc[t-1])/(ce[t] + ge[t])
        pk[t] = (1 + muk)*(Wk[t] + Tek[t] + Tfke[t] + r*Lk[t-1])/(ice[t] + ige[t])
        Δpk[t] = pk[t] - pk[t-1]
        wg[t] = wg[t-1]*(1 + δ1[t] - δ2*(p[t]-p[t-1])/p[t-1])
        Wg[t] = p[t]*wg[t]
        W[t] = Wc[t] + Wk[t] + Wb[t] + Wg[t]
        ΔW[t] = W[t] - W[t-1]
        gD[t] = gD[t-1]*(1 + δ1[t] - δ2*(p[t]-p[t-1])/p[t-1])
        GD[t] = p[t]*gD[t]
        CwD[t] = ((1 - α6)*CwD[t-1] + α6*max(0, α1*(We[t] - Tiwe[t] - Tewe[t] - r*Lw[t-1]) + α2*(Mw[t-1] + Hw[t-1] - Lw[t-1])))*Cw[t-1]/(Cw[t-1]-ΔCw[t-1])
        CiD[t] = ((1 - α6)*CiD[t-1] + α6*max(0, α3*(Πcie[t] + Πkie[t] - Tiie[t] - Teie[t]) + (α4 + α5*Mi[t-1]*(Eci[t-1] + Eki[t-1])/(Ec[t-1] + Ek[t-1]))))*Ci[t-1]/(Ci[t-1] - ΔCi[t-1])
        Cw[t] = CwD[t]/max(1, (CwD[t] + CiD[t] + GD[t])/(p[t]*min(ζ1*kc[t-1], ζ2[t])))
        Ci[t] = CiD[t]/max(1, (CwD[t] + CiD[t] + GD[t])/(p[t]*min(ζ1*kc[t-1], ζ2[t])))
        c[t], C[t] = (Cw[t] + Ci[t])/p[t], Cw[t] + Ci[t]
        Δc[t] = c[t] - c[t-1]
        G[t] = GD[t]/max(1, (CwD[t] + CiD[t] + GD[t])/(p[t]*min(ζ1*kc[t-1], ζ2[t])))
        cw[t], ci[t], g[t] = Cw[t]/p[t], Ci[t]/p[t], G[t]/p[t]
        Δg[t] = g[t] - g[t-1]
        uc[t] = (c[t] + g[t])/(ζ1*kc[t-1])
        icD[t] = btw(0, (uc[t-1] - uT)*kc[t-1] + β1*kc[t-1], max(0, β2*(Mc[t-1] - Lc[t-1])/pk[t]))
        ikD[t] = btw(0, (uk[t-1] - uT)*kk[t-1] + β1*kk[t-1], max(0, β2*(Mk[t-1] - Lk[t-1])/pk[t]))
        igD[t] = igD[t-1]*(1 + δ1[t] - δ2*(pk[t]-pk[t-1])/pk[t-1])
        ic[t] = icD[t]/max(1, (icD[t] + ikD[t] + igD[t])/(min(ζ5*kk[t-1], ζ3[t])))
        ik[t] = ikD[t]/max(1, (icD[t] + ikD[t] + igD[t])/(min(ζ5*kk[t-1], ζ3[t])))
        ig[t] = igD[t]/max(1, (icD[t] + ikD[t] + igD[t])/(min(ζ5*kk[t-1], ζ3[t])))
        Δic[t], Δig[t] = ic[t] - ic[t-1], ig[t] - ig[t-1]
        Ic[t], Ik[t], Ig[t], If[t] = pk[t]*ic[t], pk[t]*ik[t], pk[t]*ig[t], pk[t]*ic[t] + pk[t]*ik[t]
        uk[t] = (If[t]/pk[t] + ig[t])/(ζ5*kk[t-1])
        kc[t], kk[t], kg[t] = (1 - β1)*kc[t-1] + ic[t], (1 - β1)*kk[t-1] + ik[t], (1 - β1)*kg[t-1] + ig[t]
        Kc[t], Kk[t], Kg[t] = pk[t]*kc[t], pk[t]*kk[t], pk[t]*kg[t]
        Δkc[t], Δkk[t], Δkg[t] = kc[t] - kc[t-1], kk[t] - kk[t-1], kg[t] - kg[t-1]
        ΔKc[t], ΔKk[t], ΔKg[t] = Kc[t] - Kc[t-1], Kk[t] - Kk[t-1], Kg[t] - Kg[t-1]
        ΔK[t], K[t], Δk[t], k[t] = ΔKc[t] + ΔKk[t] + ΔKg[t], Kc[t] + Kk[t] + Kg[t], Δkc[t] + Δkk[t] + Δkg[t], kc[t] + kk[t] + kg[t]
        Δec[t] = max(0, κ1*(pk[t]*icD[t] - β1*Kc[t-1])/pec[t-1]) - max(0, κ2*(Mc[t-1] - Lc[t-1])/pec[t-1] - κ3*(c[t]+g[t]))
        #TODO 株式発行部数の増減により、企業の資金調達のインプットーアウトプットの大きさが実質消費に比例する程度にできないか
        ec[t] = ec[t-1] + Δec[t]
        Δek[t] = max(0, κ1*(pk[t]*ikD[t] - β1*Kk[t-1])/pek[t-1]) - max(0, κ2*(Mk[t-1] - Lk[t-1])/pek[t-1] - κ3*(ic[t] + ig[t]))
        ek[t] = ek[t-1] + Δek[t]
        Tii[t] = τ1*(Πci[t-1] + Πki[t-1])
        ΔTii[t] = Tii[t] - Tii[t-1]
        Tiw[t] = τ1*W[t-1]
        ΔTiw[t] = Tiw[t] - Tiw[t-1]
        Ti[t] = Tii[t] + Tiw[t]
        Tfc[t] = τ2*(C[t] + G[t] - Wc[t] - Tec[t] - r*Lc[t-1])
        ΔTfc[t] = Tfc[t] - Tfc[t-1]
        Tfk[t] = τ2*(If[t] + Ig[t] - Wk[t] - Tek[t] - r*Lk[t-1])
        ΔTfk[t] = Tfk[t] - Tfk[t-1]
        Πc[t] = C[t] + G[t] - Wc[t] -Tec[t] - Tfc[t] - r*Lc[t-1]
        Πci[t] = max(0, (θ1*(Πc[t] - Ic[t]) + θ2*(Mc[t-1] - Lc[t-1]))*Eci[t-1]/Ec[t-1])
        Πcb[t] = max(0, (θ1*(Πc[t] - Ic[t]) + θ2*(Mc[t-1] - Lc[t-1]))*Ecb[t-1]/Ec[t-1])
        Πcc[t] = Πc[t] - Πci[t] - Πcb[t]
        Πk[t] = If[t] + Ig[t] - Wk[t] - Tek[t] - Tfk[t] - r*Lk[t-1]
        Πki[t] = max(0, (θ1*(Πk[t] - Ik[t]) + θ2*(Mk[t-1] - Lk[t-1]))*Eki[t-1]/Ek[t-1])
        Πkb[t] = max(0, (θ1*(Πk[t] - Ik[t]) + θ2*(Mk[t-1] - Lk[t-1]))*Ekb[t-1]/Ek[t-1])
        Πkk[t] = Πk[t] - Πki[t] - Πkb[t]
        Tfb[t] = max(0, τ2*(Πcb[t] + Πkb[t] + r*L[t-1] - Wb[t] + Δpec[t]*ecb[t-1] + Δpek[t]*ekb[t-1]))
        Tf[t] = Tfc[t] + Tfk[t] + Tfb[t]
        GS[t] = -G[t] - Wg[t] + Ti[t] + Te[t] + Tf[t]

        NLw[t] = -Cw[t] + W[t] - Tiw[t] - Tew[t] - r*Lw[t-1]
        NLi[t] = -Ci[t] - Tii[t] - Tei[t] + Πci[t] + Πki[t]
        NLc[t] = Πcc[t] - Ic[t]
        NLk[t] = Πkk[t] - Ik[t]
        NLb[t] = -Wb[t] - Tfb[t] + Πcb[t] + Πkb[t] + r*L[t-1]
        NLg[t] = -Ig[t] + GS[t]

        Hw[t] = ι1*Cw[t]
        ΔHw[t] = Hw[t] - Hw[t-1]
        Lw[t] = ι2*W[t]
        ΔLw[t] = Lw[t] - Lw[t-1]
        ΔMw[t] = NLw[t] + ΔLw[t] - ΔHw[t]
        Mw[t] = Mw[t-1] + ΔMw[t]
        EciT[t] = κ4*NWie[t]*(C[t] + G[t])/(C[t] + G[t] + Ic[t] + Ig[t])
        EkiT[t] = κ4*NWie[t]*(Ic[t] + Ig[t])/(C[t] + G[t] + Ic[t] + Ig[t])
        EcbT[t] = ((1 - κ5)*Ecb[t-1] + κ5*rE[t-1]*(L[t-1] + Ecb[t-1] + Ekb[t-1])/(r + rE[t-1]))*Ec[t-1]/(Ec[t-1] + Ek[t-1])
        EkbT[t] = ((1 - κ5)*Ekb[t-1] + κ5*rE[t-1]*(L[t-1] + Ecb[t-1] + Ekb[t-1])/(r + rE[t-1]))*Ek[t-1]/(Ec[t-1] + Ek[t-1])
        pec[t] = (EciT[t] + EcbT[t])/ec[t]
        pek[t] = (EkiT[t] + EkbT[t])/ek[t]
        Δpec[t], Δpek[t] = pec[t] - pec[t-1], pek[t] - pek[t-1]
        Ec[t], Ek[t] = pec[t]*ec[t], pek[t]*ek[t]
        rE[t] = (Πci[t] + Πcb[t] + Πki[t] + Πkb[t])/(Ec[t] + Ek[t])
        Eci[t] = min(κ4*NWie[t]*(C[t] + G[t])/(C[t] + G[t] + Ic[t] + Ig[t]), Ec[t])
        Eki[t] = min(κ4*NWie[t]*(Ic[t] + Ig[t])/(C[t] + G[t] + Ic[t] + Ig[t]), Ek[t])
        Ecb[t], Ekb[t] = Ec[t] - Eci[t], Ek[t] - Eki[t]
        eci[t], ecb[t], ec[t], eki[t], ekb[t], ek[t] = Eci[t]/pec[t], Ecb[t]/pec[t], Ec[t]/pec[t], Eki[t]/pek[t], Ekb[t]/pek[t], Ek[t]/pek[t]
        
        NWw[t] = NWw[t-1] + NLw[t]
        NWi[t] = NWi[t-1] + NLi[t] + Δpec[t]*eci[t-1] + Δpek[t]*eki[t-1]
        NWc[t] = NWc[t-1] + NLc[t] - Δpec[t]*ec[t-1] + ΔKc[t]
        NWk[t] = NWk[t-1] + NLk[t] - Δpek[t]*ek[t-1] + ΔKk[t]
        NWb[t] = NWb[t-1] + NLb[t] + Δpec[t]*ecb[t-1] + Δpek[t]*ekb[t-1]
        NWg[t] = NWg[t-1] + NLg[t] + ΔKg[t]

        Δeci[t], Δecb[t], Δeki[t], Δekb[t] = eci[t] - eci[t-1], ecb[t] - ecb[t-1], eki[t] - eki[t-1], ekb[t] - ekb[t-1]
        ΔEci[t], ΔEcb[t], ΔEc[t], ΔEki[t], ΔEkb[t], ΔEk[t] = Eci[t] - Eci[t-1], Ecb[t] - Ecb[t-1], Ec[t] - Ec[t-1], Eki[t] - Eki[t-1], Ekb[t] - Ekb[t-1], Ek[t] - Ek[t-1]
        Mi[t] = NWi[t] - Eci[t] - Eki[t]
        ΔMi[t] = Mi[t] - Mi[t-1]
        LcD[t] = (1 - ι4)*LcD[t-1] + ι4*max(0, ι3*(Ic[t] + Wc[t] + Tec[t] + Tfc[t] + r*Lc[t-1]) - Mc[t-1])
        Lc[t] = btw(0, LcD[t], max(0, ι5*NLc[t]))
        ΔLc[t] = Lc[t] - Lc[t-1]
        ΔMc[t] = NLc[t] + ΔLc[t] + pec[t]*Δec[t]
        Mc[t] = Mc[t-1] + ΔMc[t]
        LkD[t] = (1 - ι4)*LkD[t-1] + ι4*max(0, ι3*(Wk[t] + Tek[t] + Tfk[t] + r*Lk[t-1]) - Mc[t-1])
        Lk[t] = btw(0, LkD[t], max(0, ι5*NLk[t]))
        ΔLk[t] = Lk[t] - Lk[t-1]
        ΔMk[t] = NLk[t] + ΔLk[t] + pek[t]*Δek[t]
        Mk[t] = Mk[t-1] + ΔMk[t]
        M[t] = Mw[t] + Mi[t] + Mc[t] + Mk[t]
        ΔM[t] = ΔMw[t] + ΔMi[t] + ΔMc[t] + ΔMk[t]
        ΔL[t] = ΔLw[t] + ΔLc[t] + ΔLk[t]
        L[t] = Lw[t] + Lc[t] + Lk[t]
        ΔHb[t] = NLb[t] + ΔM[t] - ΔL[t] - pec[t]*Δecb[t] - pek[t]*Δekb[t]
        Hb[t] = NWb[t] + M[t] - Ecb[t] - Ekb[t] - L[t]
        H[t] = Hw[t] + Hb[t]
        ΔH[t] = ΔHw[t] + ΔHb[t]

    end

    #   会計的整合性の確認
    println("--------final--------")
    println("-Ci-Tii-Tei+Πci+Πki-ΔMi-pecΔeci-pekΔeki=",-Ci[end]-Tii[end]-Tei[end]+Πci[end]+Πki[end]-ΔMi[end]-pec[end]*Δeci[end]-pek[end]*Δeki[end])
    println("-Ig+GS+ΔH=",-Ig[end]+GS[end]+ΔH[end])
    println("NLw+NLi+NLc+NLk+NLb+NLg=",NLw[end]+NLi[end]+NLc[end]+NLk[end]+NLb[end]+NLg[end])
    println("Δeci+Δecb-Δec=",Δeci[end]+Δecb[end]-Δec[end])
    println("Δeki+Δekb-Δek=",Δeki[end]+Δekb[end]-Δek[end])
    println("K-K_{-1}-ΔK=",K[end]-K[end-1]-ΔK[end])
    println("ΔKc-Δpk*k_{c-1}-pk*Δkc=",ΔKc[end]-Δpk[end]*kc[end-1]-pk[end]*Δkc[end])
    println("ΔKk-Δpk*k_{k-1}-pk*Δkk=",ΔKk[end]-Δpk[end]*kk[end-1]-pk[end]*Δkk[end])
    println("ΔKg-Δpk*k_{g-1}-pk*Δkg=",ΔKg[end]-Δpk[end]*kg[end-1]-pk[end]*Δkg[end])
    println("ΔK-Δpk*k_{-1}-pk*Δk=",ΔK[end]-Δpk[end]*k[end-1]-pk[end]*Δk[end])
    println("M-M_{-1}-ΔM=",M[end]-M[end-1]-ΔM[end])    
    println("L-L_{-1}-ΔL=",L[end]-L[end-1]-ΔL[end])
    println("ΔEci-Δpec*e_{ci-1}-pec*Δeci=",ΔEci[end]-Δpec[end]*eci[end-1]-pec[end]*Δeci[end])
    println("ΔEki-Δpek*e_{ki-1}-pek*Δeki=",ΔEki[end]-Δpek[end]*eki[end-1]-pek[end]*Δeki[end])
    println("ΔEcb-Δpec*e_{cb-1}-pec*Δecb=",ΔEcb[end]-Δpec[end]*ecb[end-1]-pec[end]*Δecb[end])
    println("ΔEkb-Δpek*e_{kb-1}-pek*Δekb=",ΔEkb[end]-Δpek[end]*ekb[end-1]-pek[end]*Δekb[end])
    println("ΔEc-Δpec*e_{c-1}-pec*Δec=",ΔEc[end]-Δpec[end]*ec[end-1]-pec[end]*Δec[end])
    println("ΔEk-Δpek*e_{k-1}-pek*Δek=",ΔEk[end]-Δpek[end]*ek[end-1]-pek[end]*Δek[end])
    println("Hb-H_{b-1}-ΔHb=",Hb[end]-Hb[end-1]-ΔHb[end])
    println("H-H_{-1}-ΔH=",H[end]-H[end-1]-ΔH[end])
    println("ec-eci-ecb=",ec[end]-eci[end]-ecb[end])
    println("ek-eki-ekb=",ek[end]-eki[end]-ekb[end])
    println("NWw+NWi+NWc+NWk+NWb+NWg-K=",NWw[end]+NWi[end]+NWc[end]+NWk[end]+NWb[end]+NWg[end]-K[end])
    println("-NWw+Mw-Lw+Hw=",-NWw[end]+Mw[end]-Lw[end]+Hw[end])
    println("-NWc+Kc+Mc-Ec-Lc=",-NWc[end]+Kc[end]+Mc[end]-Ec[end]-Lc[end])
    println("-NWk+Kk+Mk-Ek-Lk=",-NWc[end]+Kc[end]+Mc[end]-Ec[end]-Lc[end])
    println("-NWg+Kg-H=",-NWg[end]+Kg[end]-H[end])
end

run()

function plot_transition()
    plot(p[end-150:end], label="p")
    plot!(pk[end-150:end], label="pk")
    savefig("figs/p.png")

    plot(pec[end-150:end], label="pec")
    plot!(pek[end-150:end], label="pek")
    savefig("figs/pe.png")

    plot(NLw[end-150:end], label="NLw")
    plot!(NLi[end-150:end], label="NLi")
    plot!(NLc[end-150:end], label="NLc")
    plot!(NLk[end-150:end], label="NLk")
    plot!(NLb[end-150:end], label="NLb")
    plot!(NLg[end-150:end], label="NLg")
    savefig("figs/NL.png")

    plot(uc[end-150:end], label="uc")
    plot!(uk[end-150:end], label="uk")
    savefig("figs/u.png")

    plot(W[end-150:end], label="W")
    plot!(Wc[end-150:end], label="Wc")
    plot!(Wk[end-150:end], label="Wk")
    plot!(Wb[end-150:end], label="Wb")
    plot!(Wg[end-150:end], label="Wg")
    savefig("figs/W.png")

    plot(Πc[end-150:end], label="Πc")
    plot!(Πcc[end-150:end], label="Πcc")
    plot!(Πcb[end-150:end], label="Πcb")
    plot!(Πci[end-150:end], label="Πci")
    savefig("figs/Πc.png")

    plot(Πk[end-150:end], label="Πk")
    plot!(Πkk[end-150:end], label="Πkf")
    plot!(Πkb[end-150:end], label="Πkb")
    plot!(Πki[end-150:end], label="Πki")
    savefig("figs/Πk.png")

    plot(kc[end-150:end], label="kc")
    plot!(kk[end-150:end], label="kk")
    plot!(kg[end-150:end], label="kg")
    savefig("figs/k.png")

    plot(ic[end-150:end], label="ic")
    plot!(ik[end-150:end], label="ik")
    plot!(ig[end-150:end], label="ig")
    savefig("figs/i.png")

    plot(C[end-150:end], label="C")
    plot!(Ci[end-150:end], label="Ci")
    plot!(Cw[end-150:end], label="Cw")
    savefig("figs/nominal_C.png")

    plot(ζ2[end-150:end], label="ζ2")
    savefig("figs/ζ2.png")

    plot(C[end-150:end], label="C")
    plot!(Ic[end-150:end], label="Ic")
    plot!(Ik[end-150:end], label="Ik")
    plot!(Ig[end-150:end], label="Ig")
    plot!(G[end-150:end], label="G")
    savefig("figs/nominal_Y.png")

    plot(c[end-150:end], label="c")
    plot!(ic[end-150:end], label="ic")
    plot!(ik[end-150:end], label="ik")
    plot!(ig[end-150:end], label="ig")
    plot!(g[end-150:end], label="g")
    savefig("figs/y.png")

    plot(M[end-150:end], label="M")
    plot!(Mw[end-150:end], label="Mw")
    plot!(Mi[end-150:end], label="Mi")
    plot!(Mc[end-150:end], label="Mc")
    plot!(Mk[end-150:end], label="Mk")
    savefig("figs/M.png")

    plot(L[end-150:end], label="L")
    plot!(Lw[end-150:end], label="Lw")
    plot!(Lc[end-150:end], label="Lc")
    plot!(Lk[end-150:end], label="Lk")
    savefig("figs/L.png")

    plot(NWw[end-150:end], label="NWw")
    plot!(NWi[end-150:end], label="NWi")
    plot!(NWc[end-150:end], label="NWc")
    plot!(NWk[end-150:end], label="NWk")
    plot!(NWb[end-150:end], label="NWb")
    plot!(NWg[end-150:end], label="NWg")
    savefig("figs/NW.png")

    plot(Eci[end-150:end], label="Eci")
    plot!(Ec[end-150:end], label="Ec")
    plot!(Ecb[end-150:end], label="Ecb")
    savefig("figs/nominal_Ec.png")

    plot(Eki[end-150:end], label="Eki")
    plot!(Ek[end-150:end], label="Ek")
    plot!(Ekb[end-150:end], label="Ekb")
    savefig("figs/nominal_Ek.png")

    plot(eci[end-150:end], label="eci")
    plot!(ecb[end-150:end], label="ecb")
    plot!(ec[end-150:end], label="ec")
    savefig("figs/ec.png")

    plot(eki[end-150:end], label="eki")
    plot!(ekb[end-150:end], label="ekb")
    plot!(ek[end-150:end], label="ek")
    savefig("figs/ek.png")

    plot(Hw[end-150:end], label="Hw")
    plot!(Hb[end-150:end], label="Hb")
    plot!(H[end-150:end], label="H")
    savefig("figs/H.png")

    plot(Tiw[end-150:end], label="Tiw")
    plot!(Tii[end-150:end], label="Tii")
    plot!(Ti[end-150:end], label="Ti")
    savefig("figs/Ti.png")

    plot(Tew[end-150:end], label="Tew")
    plot!(Tei[end-150:end], label="Tei")
    plot!(Tec[end-150:end], label="Tec")
    plot!(Tek[end-150:end], label="Tek")
    plot!(Te[end-150:end], label="Te")
    savefig("figs/Te.png")

    plot(Tfb[end-150:end], label="Tfb")
    plot!(Tfc[end-150:end], label="Tfc")
    plot!(Tfk[end-150:end], label="Tfk")
    plot!(Tf[end-150:end], label="Tf")
    savefig("figs/Tf.png")

    plot(Πci[end-150:end], label="Πci")
    plot!(Πc[end-150:end], label="Πc")
    plot!(Πcc[end-150:end], label="Πcc")
    plot!(Πcb[end-150:end], label="Πcb")
    savefig("figs/Πc.png")

    plot(Πki[end-150:end], label="Πki")
    plot!(Πk[end-150:end], label="Πk")
    plot!(Πkk[end-150:end], label="Πkk")
    plot!(Πkb[end-150:end], label="Πkb")
    savefig("figs/Πc.png")

    plot(cw[end-150:end], label="cw")
    plot!(ci[end-150:end], label="ci")
    plot!(c[end-150:end], label="c")
    savefig("figs/c.png")
end
plot_transition()

function save_mock_data()
    mock_C = deepcopy(C)
    mock_c = deepcopy(c)
    mock_I = deepcopy(I)
    mock_i = deepcopy(i)
    mock_p = deepcopy(p)
    mock_W = deepcopy(W)
    mock_G = deepcopy(G)
    mock_g = deepcopy(g)
    mock_pe = deepcopy(pe)
    mock_K = deepcopy(K)
    mock_k = deepcopy(k)
    mock_NLw = deepcopy(NLw)
    mock_NLi = deepcopy(NLi)
    mock_NLf = deepcopy(NLf)
    mock_NLb = deepcopy(NLb)
    mock_NLg = deepcopy(NLg)
    mock_M = deepcopy(M)
    mock_H = deepcopy(H)
    mock_L = deepcopy(L)
    mock_Te = deepcopy(Te)
    mock_Tf = deepcopy(Tf)
    mock_Ti = deepcopy(Ti)
    mock_u = deepcopy(u)
    mock_Π = deepcopy(Π)
    mock_NWw = deepcopy(NWw)
    mock_NWi = deepcopy(NWi)
    mock_NWf = deepcopy(NWf)
    mock_NWb = deepcopy(NWb)
    mock_NWg = deepcopy(NWg)
    return mock_C, mock_c, mock_I, mock_i, mock_p, mock_W, mock_G, mock_g, mock_pe, mock_K, mock_k, mock_NLw, mock_NLi, mock_NLf, mock_NLb, mock_NLg, mock_M, mock_H, mock_L, mock_Te, mock_Tf, mock_Ti, mock_u, mock_Π, mock_NWw, mock_NWi, mock_NWf, mock_NWb, mock_NWg
end
function plot_mock_and_control()
    plot(mock_C[end-150:end]./mock_G[end-150], label="mock")
    plot!(C, label="control")
    savefig("figs/diff_nominal_C.png")
    plot(mock_c[end-150:end]./mock_g[end-150], label="mock")
    plot!(c[end-150:end]./g[end-150], label="control")
    savefig("figs/diff_c.png")
    plot(mock_I[end-150:end]./mock_G[end-150], label="mock")
    plot!(I, label="control")
    savefig("figs/diff_nominal_I.png")
    plot(mock_i[end-150:end]./mock_g[end-150], label="mock")
    plot!(i[end-150:end]./g[end-150], label="control")
    savefig("figs/diff_i.png")
    plot(mock_p[end-150:end]./mock_p[end-150], label="mock")
    plot!(p[end-150:end]./p[end-150], label="control")
    savefig("figs/diff_p.png")
    plot(mock_W[end-150:end]./mock_G[end-150], label="mock")
    plot!(W, label="control")
    savefig("figs/diff_nominal_W.png")
    plot(mock_W[end-150:end]./mock_p[end-150:end]./mock_G[end-150], label="mock")
    plot!(W[end-150:end]./p, label="control")
    savefig("figs/diff_real_W.png")
    plot(mock_G[end-150:end]./mock_G[end-150], label="mock")
    plot!(G, label="control")
    savefig("figs/diff_nominal_G.png")
    plot(mock_g[end-150:end]./mock_g[end-150], label="mock")
    plot!(g[end-150:end]./g[end-150], label="control")
    savefig("figs/diff_g.png")
    plot(mock_pe[end-150:end]./mock_pe[end-150], label="mock")
    plot!(pe[end-150:end]./pe[end-150], label="control")
    savefig("figs/diff_pe.png")
    plot(mock_K[end-150:end]./mock_G[end-150], label="mock")
    plot!(K, label="control")
    savefig("figs/diff_nominal_K.png")
    plot(mock_k[end-150:end]./mock_g[end-150], label="mock")
    plot!(k[end-150:end]./g[end-150], label="control")
    savefig("figs/diff_k.png")
    plot(mock_NLw[end-150:end]./mock_G[end-150], label="mock")
    plot!(NLw, label="control")
    savefig("figs/diff_NLw.png")
    plot(mock_NLi[end-150:end]./mock_G[end-150], label="mock")
    plot!(NLi, label="control")
    savefig("figs/diff_NLi.png")
    plot(mock_NLf[end-150:end]./mock_G[end-150], label="mock")
    plot!(NLf, label="control")
    savefig("figs/diff_NLf.png")
    plot(mock_NLb[end-150:end]./mock_G[end-150], label="mock")
    plot!(NLb, label="control")
    savefig("figs/diff_NLb.png")
    plot(mock_NLg[end-150:end]./mock_G[end-150], label="mock")
    plot!(NLg, label="control")
    savefig("figs/diff_NLg.png")
    plot(mock_M[end-150:end]./mock_G[end-150], label="mock")
    plot!(M, label="control")
    savefig("figs/diff_M.png")
    plot(mock_H[end-150:end]./mock_G[end-150], label="mock")
    plot!(H, label="control")
    savefig("figs/diff_H.png")
    plot(mock_L[end-150:end]./mock_G[end-150], label="mock")
    plot!(L, label="control")
    savefig("figs/diff_L.png")
    plot(mock_Te[end-150:end]./mock_G[end-150], label="mock")
    plot!(Te, label="control")
    savefig("figs/diff_Te.png")
    plot(mock_Tf[end-150:end]./mock_G[end-150], label="mock")
    plot!(Tf, label="control")
    savefig("figs/diff_Tf.png")
    plot(mock_Ti[end-150:end]./mock_G[end-150], label="mock")
    plot!(Ti, label="control")
    savefig("figs/diff_Ti.png")
    plot(mock_u, label="mock")
    plot!(u, label="control")
    savefig("figs/diff_u.png")
    plot(mock_Π[end-150:end]./mock_G[end-150], label="mock")
    plot!(Π, label="control")
    savefig("figs/diff_Π.png")
    plot(mock_NWw[end-150:end]./mock_G[end-150], label="mock")
    plot!(NWw, label="control")
    savefig("figs/diff_NWw.png")
    plot(mock_NWi[end-150:end]./mock_G[end-150], label="mock")
    plot!(NWi, label="control")
    savefig("figs/diff_NWi.png")
    plot(mock_NWf[end-150:end]./mock_G[end-150], label="mock")
    plot!(NWf, label="control")
    savefig("figs/diff_NWf.png")
    plot(mock_NWb[end-150:end]./mock_G[end-150], label="mock")
    plot!(NWb, label="control")
    savefig("figs/diff_NWb.png")
    plot(mock_NWg[end-150:end]./mock_G[end-150], label="mock")
    plot!(NWg, label="control")
    savefig("figs/diff_NWg.png")
end

#mock_C, mock_c, mock_I, mock_i, mock_p, mock_W, mock_G, mock_g, mock_pe, mock_K, mock_k, mock_NLw, mock_NLi, mock_NLf, mock_NLb, mock_NLg, mock_M, mock_H, mock_L, mock_Te, mock_Tf, mock_Ti, mock_u, mock_Π, mock_NWw, mock_NWi, mock_NWf, mock_NWb, mock_NWg = save_mock_data()
#   初期化

#   パラメータ変更

#   controlのシミュレーションと、mockとcontrolの比較のプロット
#run()
#plot_mock_and_control()
