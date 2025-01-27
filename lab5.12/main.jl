using LinearAlgebra
using Plots
function gen_rand_simetric(n, b, a)
    m = rand(n, n) * (b-a) .+ a
    return Symmetric(m)
end
function norm(v)
    res = 0
    for i in eachindex(v)
        res += v[i]^2
    end
    return sqrt(res)
end
function gen_identiti(n)
    return Matrix(1.0*I, n, n)
end
function Gershgorian(m)
    g = copy(m)
    n = size(m, 1)
    diag = [m[i, i] for i in 1:n]
    r = [sum(abs, g[i, :])-diag[i] for i in 1:n]
    inter = [[diag[i]-r[i], diag[i] + r[i]] for i in 1:n]
    return inter
end
function Krilov(m)
    n = size(m, 1)
    y_s = []
    push!(y_s, ones(n, 1))
    for _ in 1:n
        push!(y_s, A * y_s[end])
    end
    B = hcat(y_s[n:-1:1]...)
    p = B\y_s[end]
    return p, y_s
end
function poli(p, x)
    n = size(p, 1)
    return sum([-p[i]*x^(n-i) for i in 1:n])+x^n
end
function Danilevskiy(m)
    n = size(m, 1)
    B = gen_identiti(n)
    D = copy(m)
    for i in n:-1:2
        B_i = gen_identiti(n)
        B_i[i-1, :] = -D[i, :] / D[i, i -1]
        B_i[i-1,i-1] = 1 / D[i, i -1]
        D = inv(B_i) * D * B_i
        B *= B_i
    end
    return B, D
end
function Durand_Kerner(p, start, eps)
    new = []
    while true
        for i in 1:length(start)
            znam = 1
            for j in 1:length(new)
                znam *= (start[i]-new[j])
            end
            for j in length(new)+1+1:length(start)
                znam *= (start[i]-start[j])
            end
            r = start[i] - poli(p, start[i]) / znam
            push!(new, r)
        end
        if maximum(abs, new - start) < eps
            return new
        else
            start = new
            new = []
        end
    end
end
function find_vectors_D(B, spectr)
    vects = []
    n = size(B, 1)
    for lambda in spectr
        p_vec = [lambda^i for i in n-1:-1:0]
        push!(vects, B*p_vec)
    end
    return vects
end
function find_vectors(y_s, p, spectr)
    n = size(p, 1)
    vects = []
    for i in 1:n
        v = y_s[n]
        q = [1.0]
    for j in 1:n-1
        push!(q, spectr[i]*q[j]-p[j])
        v+= q[j+1]* y_s[n-j]
    end
        push!(vects, v)
    end
    return vects
end
println("Dan Matrix Krilov =======================")
A = [2.2 1. 0.5 2.; 1 1.3 2. 1.; 0.5 2 0.5 1.6; 2. 1. 1.6 2]
P, y_s = Krilov(A)
inter = Gershgorian(A)
maxa = minimum([i[1] for i in inter])
maxb = maximum([i[2] for i in inter])
x = range(maxa, maxb, length=10000)
savefig(plot(x, [poli(P, i) for i in x]), "myplot.pdf")
spec = Durand_Kerner(P, [(i[2]+i[1])/2 for i in inter], 10^-1)
println("Spec ",spec)
vectors = find_vectors(y_s,P, spec)
norm_vect = []
for v in vectors
    push!(norm_vect, v / norm(v))
end
println("Ortogonal")
for v in 1:length(norm_vect)-1
    for w in v+1:length(norm_vect)
        if v != w
            println(dot(norm_vect[v], norm_vect[w]))
        end
    end
end
println("|Sum(spec)-tr(A)| ", abs(sum(spec)- tr(A)))
println("Dan Matrix Danilevskiy =======================")
A = [2.2 1. 0.5 2.; 1 1.3 2. 1.; 0.5 2 0.5 1.6; 2. 1. 1.6 2]
B, D = Danilevskiy(A)
P = D[1, :]
inter = Gershgorian(A)
maxa = minimum([i[1] for i in inter])
maxb = maximum([i[2] for i in inter])
x = range(maxa, maxb, length=10000)
savefig(plot(x, [poli(P, i) for i in x]), "myplot1.pdf")
spec = Durand_Kerner(P, [(i[2]+i[1])/2 for i in inter], 10^-1)
println("Spec ",spec)
vectors = find_vectors_D(B, spec)
norm_vect = []
for v in vectors
    push!(norm_vect, v / norm(v))
end
println("Ortogonal")
for v in 1:length(norm_vect)-1
    for w in v+1:length(norm_vect)
        if v != w
            println(dot(norm_vect[v], norm_vect[w]))
        end
    end
end
println("|Sum(spec)-tr(A)| ", abs(sum(spec)- tr(A)))
println("Rand Matrix Krilov =======================")
A = gen_rand_simetric(5, 10, -10)
P, y_s = Krilov(A)
inter = Gershgorian(A)
maxa = minimum([i[1] for i in inter])
maxb = maximum([i[2] for i in inter])
x = range(maxa, maxb, length=10000)
savefig(plot(x, [poli(P, i) for i in x]), "myplot3.pdf")
spec = Durand_Kerner(P, [(i[2]+i[1])/2 for i in inter], 10^-1)
println("Spec ",spec)
vectors = find_vectors(y_s,P, spec)
norm_vect = []
for v in vectors
    push!(norm_vect, v / norm(v))
end
println("Ortogonal")
for v in 1:length(norm_vect)-1
    for w in v+1:length(norm_vect)
        if v != w
            println(dot(norm_vect[v], norm_vect[w]))
        end
    end
end
println("|Sum(spec)-tr(A)| ", abs(sum(spec)- tr(A)))
println("Rand Matrix Danilevskiy =======================")
B, D = Danilevskiy(A)
P = D[1, :]
inter = Gershgorian(A)
maxa = minimum([i[1] for i in inter])
maxb = maximum([i[2] for i in inter])
x = range(maxa, maxb, length=10000)
savefig(plot(x, [poli(P, i) for i in x]), "myplot4.pdf")
spec = Durand_Kerner(P, [(i[2]+i[1])/2 for i in inter], 10^-1)
println("Spec ",spec)
vectors = find_vectors_D(B, spec)
norm_vect = []
for v in vectors
    push!(norm_vect, v / norm(v))
end
println("Ortogonal")
for v in 1:length(norm_vect)-1
    for w in v+1:length(norm_vect)
        if v != w
            println(dot(norm_vect[v], norm_vect[w]))
        end
    end
end
println("|Sum(spec)-tr(A)| ", abs(sum(spec)- tr(A)))