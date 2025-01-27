using LinearAlgebra

function gen_rand(n, b, a)
    return rand(n, n) * (b-a) .+ a
end

function compute_svd(A)
    AtA = A' * A
    eigvals, V = eigen(AtA)
    idx = sortperm(eigvals, rev=true)
    eigvals = eigvals[idx]
    V = V[:, idx]
    sigma = sqrt.(eigvals)
    Σ = zeros(size(A, 1), size(A, 2))
    for i in 1:minimum(size(A))
        Σ[i, i] = sigma[i]
    end
    U = zeros(size(A,1), size(A,1))
    for i in eachindex(sigma)
        U[:, i] = A * V[:, i] / sigma[i]
    end
    return U, Σ, V
end

A = gen_rand(3, 10, -10)

U, Σ, Vt = compute_svd(A)
U_jl, neΣ_jl, Vt_jl = svd(A)
Σ_jl = zeros(size(A, 1), size(A, 2))
for i in eachindex(neΣ_jl)
    Σ_jl[i, i] = neΣ_jl[i]
end

b = IOBuffer()

println("My U")
show(b, "text/plain", U)
println(String(take!(b)))

println("Jl U")
show(b, "text/plain", U_jl)
println(String(take!(b)))

println("My Σ")
show(b, "text/plain", Σ)
println(String(take!(b)))

println("Jl Σ")
show(b, "text/plain", Σ_jl)
println(String(take!(b)))

println("My Vt")
show(b, "text/plain", Vt)
println(String(take!(b)))

println("Jl Vt")
show(b, "text/plain", Vt_jl)
println(String(take!(b)))

println("----------------------------------------------------------------------------------------")
println("A")
show(b, "text/plain", A)
println(String(take!(b)))

println("Mul my")
show(b, "text/plain", U*Σ*Vt')
println(String(take!(b)))

println("Mul jl")
show(b, "text/plain", U_jl*Σ_jl*Vt_jl')
println(String(take!(b)))

println("Is U ort?")
show(b, "text/plain", U'*U)
println(String(take!(b)))

println("Is Vt ort?")
show(b, "text/plain", Vt'*Vt)
println(String(take!(b)))