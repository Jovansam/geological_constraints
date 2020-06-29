"""
Geological constraints model of resource exploitation formulated as an equilibrium problem
28-06-2020
s.j.okullo@outlook.com
"""

using JuMP,Plots
using Ipopt
#using KNITRO
#using GAMS

function geological_constraints()

    Α = 50.0
    Β = 0.1
    δ = 0.05
    γ = 0.1
    ψ₁ = 1.0
    ψ₂ = 2.0
    λ₁ = 5.0
    λ₂ = 4.0
    R₀ = 0.0
    S₀ = 10.0
    T = 40

    #gc = Model(KNITRO.Optimizer)
    gc = Model(Ipopt.Optimizer)
    #gc = Model(GAMS.Optimizer); set_optimizer_attribute(gc,"Solver","Ipopt")  # https://github.com/GAMS-dev/gams.jl

    @variable(gc,q[t=1:T]>=0)
    @variable(gc,x[t=1:T]>=0)
    @variable(gc,S[t=1:T]>=0)
    @variable(gc,R[t=1:T]>=0)

    @constraint(gc, [t=1:(T-1)], R[t+1] - R[t] == x[t] - q[t])
    @constraint(gc, [t=1:(T-1)], S[t+1] - S[t] == -x[t])
    @constraint(gc, [t=1:T], q[t] <= γ*R[t])
    @constraint(gc, [t=1], R[t] == R₀)
    @constraint(gc, [t=1], S[t] == S₀)

    @NLexpression(gc, Q[t=1:T], q[t])
    @NLexpression(gc, P[t=1:T], Α - Β*Q[t])
    @NLexpression(gc, C[t=1:T], ψ₁*(1/ψ₂)*q[t]^ψ₂)
    @NLexpression(gc, X[t=1:T], λ₁*(1/λ₂)*x[t]^λ₂)

    @NLobjective(gc, Max, sum( ((1+δ)^(-t))*(P[t]*Q[t] - C[t] - X[t]) for t ∈ 1:T))

    JuMP.optimize!(gc)

    return JuMP.value.(x),JuMP.value.(q),JuMP.value.(R),JuMP.value.(S),JuMP.value.(P),T
end

function plot_figures()
    @time x,q,R,S,P,T = geological_constraints()
    p1=plot(1:T,q,label="q")
    p2=plot(1:T,R,label="R")
    p3=plot(1:T,S,label="S")
    p4=plot(1:T,x,label="x")
    p5=plot(R,q,label="q vs R")
    p6=plot(P,label="P")
    plot(p1,p2,p3,p4,p5,p6,lw=2)
end

plot_figures()
