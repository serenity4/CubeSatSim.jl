using ModelingToolkit
using OrdinaryDiffEq

@parameters t
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
        D.(x) .~ v,
        D.(ω) .~ 0.,
        D.(v) .~ 0.
    )
end

u0 = vcat(
    q .=> [1., 0, 0, 0],
    x .=> 0.,
    ω .=> [0., 0., 0.2],
    v .=> [1., 0, 0]
)
tspan = (0.0, 100.0)
sys = ODESystem(dynamics(q, x, ω, v))
prob = ODEProblem(sys, u0, tspan)
solve(prob, Tsit5())