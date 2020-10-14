using ModelingToolkit
using OrdinaryDiffEq, DiffEqCallbacks
using BenchmarkTools
using LinearAlgebra
using DataFrames

@parameters t, I[1:3, 1:3], Iinv[1:3, 1:3], g[1:3], m, α
@variables q[1:4](t) x[1:3](t) ω[1:3](t) v[1:3](t)
@derivatives D'~t

function transition_matrix(q)
    q0, q1, q2, q3 = q
    P = [
        2*(q0^2 + q1^2) - 1 2*(q1*q2 - q0*q3)   2*(q1*q3 + q2*q0)
        2*(q1*q2 + q0*q3)   2*(q0^2 + q2^2) - 1 2*(q2*q3 - q1*q0)
        2*(q1*q3 - q2*q0)   2*(q2*q3 + q1*q0)   2*(q0^2 + q3^2) - 1
    ]
end

function quaternion_derivative(q, ω)
    a, b, c = ω
    Mₛ = [
        0 -a -b -c
        a  0  c -b
        b -c  0  a
        c  b -a  0
    ]
    1/2 * Mₛ * q
end

function dynamics(q, x, ω, v)
    P = transition_matrix(q)
    vcat(
        D.(q) .~ quaternion_derivative(q, ω),
        D.(x) .~ P * v, # world velocity
        D.(ω) .~ Iinv * cross(I * ω, ω),
        D.(v) .~ 1/m * (m * P' * g - α .* v) # body acceleration
    )
end

function postprocess(u) # dynamics_f does not work
    q, x, ω, v = u[1:4, end], u[5:7, end], u[8:10, end], u[11:13, end]
    dynamics_f(q, x, ω, v)
end

get_values(eqs) = getproperty.(getproperty.(eqs, :rhs), :value)
function extract_result(u)
    syms = [:q1, :q2, :q3, :q4, :x, :y, :z, :ωx, :ωy, :ωz, :vx, :vy, :vz]
    df = DataFrame(repeat([Float64], 13), syms)
    for i ∈ 1:size(u)[2]
        push!(df, Dict(syms .=> u[:, i]))
    end
    df
end

function quaternion_norm_callback(resid, u, p, t)
    resid[1] = @views norm(u[1:4])^2 - 1
    resid[2:13] .= 0
end

quaternion_manifold_proj = ManifoldProjection(quaternion_norm_callback)

u0 = vcat(
    q .=> [1., 0, 0, 0],
    x .=> 0.,
    ω .=> [0.1, 0.15, 0.2],
    v .=> [1., 0, 0],
)

p = vcat(
    (I .=> Diagonal([0.25, 0.2, 0.1]))...,
    # (Iinv .=> inv(I))...
    (Iinv .=> Diagonal([4, 5., 10.]))...,
    g .=> [0., 0., -9.81],
    m => 3.,
    α => 0.01
)

tspan = (0.0, 10000.0)
sys = ODESystem(dynamics(q, x, ω, v))
prob = ODEProblem(sys, u0, tspan, p)

u = solve(prob, Vern7(), save_everystep=false, callback=quaternion_manifold_proj)
df = extract_result(u)

# dynamics_f = build_function(dynamics(q, x, ω, v), q, x, ω, v; expression=Val{false})
# postprocess(u)

