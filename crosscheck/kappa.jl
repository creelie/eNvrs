# Independent non-rigorous computation of the front growth factor kappa(Lambda) (Julia).
# Inviscid front coordinates, window k = -KB..KA with Y_0 = 1, Y_1 = c at the start section,
# end section Y_2 = c Y_1; the step map renormalizes by G = Y_1 and shifts by one shell.
# The fixed point is found by iterating the step map; kappa = G at the fixed point.
# Arithmetic: Float64 by default, BigFloat (256 bits) with the argument "big".
function field!(dy, y, Rk)
    n = length(y)
    @inbounds for m in 1:n
        inc = m > 1 ? y[m-1]^2 : zero(eltype(y))
        out = m < n ? y[m]*y[m+1] : zero(eltype(y))
        dy[m] = Rk[m]*(inc - out)
    end
end
function rk4!(y, h, Rk, k1, k2, k3, k4, t)
    field!(k1, y, Rk); @. t = y + h/2*k1
    field!(k2, t, Rk); @. t = y + h/2*k2
    field!(k3, t, Rk); @. t = y + h*k3
    field!(k4, t, Rk); @. y = y + h/6*(k1 + 2k2 + 2k3 + k4)
end
function stepmap(w, Rk, i0, c, h)
    T = eltype(w); y = copy(w); n = length(y)
    k1, k2, k3, k4, t, prev = (similar(y) for _ in 1:6)
    hv(y) = y[i0+2] - c*y[i0+1]
    s = zero(T); hp = hv(y)
    while true
        prev .= y; rk4!(y, h, Rk, k1, k2, k3, k4, t)
        hn = hv(y)
        if hn >= 0
            a, b, fa, fb = zero(T), h, hp, hn
            for it in 1:80
                tau = a - fa*(b - a)/(fb - fa)
                y .= prev; rk4!(y, tau, Rk, k1, k2, k3, k4, t)
                ft = hv(y)
                ft < 0 ? (a = tau; fa = ft) : (b = tau; fb = ft)
                abs(ft) < eps(T)*1e-3 && break
            end
            s += a; y .= prev; rk4!(y, a, Rk, k1, k2, k3, k4, t)
            break
        end
        hp = hn; s += h
    end
    G = y[i0+1]
    wn = vcat(y[2:end] ./ G, zero(T))
    return wn, G, s
end
function kappa(Lam, KB, KA; T=Float64, h=1e-3, iters=600)
    Lam = T(Lam); r = Lam^(T(2)/3); c = T(1)/2
    ks = -KB:KA; Rk = [r^k for k in ks]; i0 = KB + 1
    w = [k <= 0 ? T(1.02)^k : zero(T) for k in ks]; w[i0+1] = c
    G = one(T); s = zero(T); d = one(T)
    for it in 1:iters
        wn, G, s = stepmap(w, Rk, i0, c, T(h))
        d = maximum(abs.(wn .- w)); w = wn
        d < 1e-14 && break
    end
    return G, s, d
end
big = length(ARGS) > 0 && ARGS[1] == "big"
for (Lam, KB) in [(1.2, 140), (1.3, 140), (1.5, 106), (1.7, 81), (1.8, 73), (1.87, 68), (1.88, 68)]
    t0 = time()
    G, s, d = kappa(Lam, KB, 8; h = 5e-4)
    a = 1/3 + log(G)/(2log(Lam))
    println("Lambda=$(Lam) kappa=$(round(G, digits=12)) s*=$(round(s, digits=8)) alpha*=$(round(a, digits=6)) diff=$(d) [$(round(time()-t0, digits=1))s]")
end
