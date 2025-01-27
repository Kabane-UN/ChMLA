using Base.Threads
using Plots

function gen_rand(n, b, a)
    return rand(n, n) * (b-a) .+ a
end

function normal_mat_mul(A, B)
    R = []
    for A_row in eachrow(A)
        row = Float64[]
        for B_col in eachcol(B)
            push!(row, sum(A_row .* B_col))
        end
        push!(R, row)
    end
    return hcat(R...)'
end

function strass_mat_mul(A, B, n_min)
    n = size(A, 1)
    if n <= n_min
        return normal_mat_mul(A,B)
    else
        C = zeros(n, n)
        m = div(n,  2)
        u = 1:m 
        v = m+1:n 
        P₁ = strass_mat_mul(A[u, u]+A[v, v], B[u, u]+B[v, v], n_min)
        P₂ = strass_mat_mul(A[v, u]+A[v, v], B[u, u], n_min)
        P₃ = strass_mat_mul(A[u, u], B[u, v]-B[v, v], n_min)
        P₄ = strass_mat_mul(A[v, v], B[v, u]-B[u, u], n_min)
        P₅ = strass_mat_mul(A[u, u]+A[u, v], B[v, v], n_min)
        P₆ = strass_mat_mul(A[v, u]-A[u, u], B[u, u]+B[u, v], n_min)
        P₇ = strass_mat_mul(A[u, v]-A[v, v], B[v, u]+B[v, v], n_min)
        C[u, u] = P₁+P₄-P₅+P₇
        C[u, v] = P₃+P₅
        C[v, u] = P₂+P₄
        C[v, v] = P₁+P₃-P₂+P₆
        return C
    end
end

function multy_strass_mat_mul(A, B, n_min)
    n = size(A, 1)
    if n <= n_min
        return normal_mat_mul(A,B)
    else
        C = zeros(n, n)
        m = div(n,  2)
        u = 1:m
        v = m+1:n
        tasks = []
        push!(tasks, @spawn multy_strass_mat_mul(A[u, u]+A[v, v], B[u, u]+B[v, v], n_min))
        push!(tasks, @spawn multy_strass_mat_mul(A[v, u]+A[v, v], B[u, u], n_min))
        push!(tasks, @spawn multy_strass_mat_mul(A[u, u], B[u, v]-B[v, v], n_min))
        push!(tasks, @spawn multy_strass_mat_mul(A[v, v], B[v, u]-B[u, u], n_min))
        push!(tasks, @spawn multy_strass_mat_mul(A[u, u]+A[u, v], B[v, v], n_min))
        push!(tasks, @spawn multy_strass_mat_mul(A[v, u]-A[u, u], B[u, u]+B[u, v], n_min))
        push!(tasks, @spawn multy_strass_mat_mul(A[u, v]-A[v, v], B[v, u]+B[v, v], n_min))
        P₁ ,P₂,P₃,P₄,P₅, P₆,P₇ = fetch.(tasks)
        
        C[u, u] = P₁+P₄-P₅+P₇
        C[u, v] = P₃+P₅
        C[v, u] = P₂+P₄
        C[v, v] = P₁+P₃-P₂+P₆
        return C
    end
end
function wine_mut_mul(A, B)
    n = size(A, 1)
    C = zeros(n, n)
    m = n ÷ 2
    rowFactor = [0.0 for _ in 1:n]
    for i in 1:n 
        rowFactor[i] = A[i, 1]*A[i, 2]
        for j in 2:m
            rowFactor[i] = rowFactor[i]+A[i, 2*j-1]*A[i, 2*j]
        end
    end
    colFactor = [0.0 for _ in 1:n]
    for i in 1:n 
        colFactor[i] = B[1, i]*B[2, i]
        for j in 2:m
            colFactor[i] = colFactor[i]+B[2*j-1, i]*B[2*j, i]
        end
    end
    for i ∈ 1:n 
        for j ∈ 1:n 
            C[i, j] += sum([(A[i, 2*k-1]+B[2*k, j])*(B[2*k-1, j]+A[i, 2*k]) for k in 1:m])-rowFactor[i]-colFactor[j]
        end
    end
    return C
end
function wine_strass_mat_mul(A, B, n_min)
    n = size(A, 1)
    if n <= n_min
        return normal_mat_mul(A,B)
    else
        C = zeros(n, n)
        m = div(n,  2)
        u = 1:m 
        v = m+1:n 
        S₁ = (A[v, u]+A[v, v])
        S₂ = (S₁ - A[u, u])
        S₃ = (A[u, u]-A[v, u])
        S₄ = (A[u, v] - S₂)
        S₅ = (B[u, v]-B[u, u])
        S₆ = (B[v, v] - S₅)
        S₇ = (B[v, v]-B[u, v])
        S₈ = (S₆-B[v, u])
        P₁ = wine_strass_mat_mul(S₂, S₆, n_min)
        P₂ = wine_strass_mat_mul(A[u, u], B[u, u], n_min)
        P₃ = wine_strass_mat_mul(A[u, v], B[v, u], n_min)
        P₄ = wine_strass_mat_mul(S₃, S₇, n_min)
        P₅ = wine_strass_mat_mul(S₁, S₅, n_min)
        P₆ = wine_strass_mat_mul(S₄, B[v, v], n_min)
        P₇ = wine_strass_mat_mul(A[v, v], S₈, n_min)
        T₁ = P₁+P₂
        T₂ = T₁+P₄
        C[u, u] = P₂+P₃
        C[u, v] = T₁+P₅+P₆
        C[v, u] = T₂-P₇
        C[v, v] = T₂+P₅
        return C
    end
end
function multy_wine_strass_mat_mul(A, B, n_min)
    n = size(A, 1)
    if n <= n_min
        return normal_mat_mul(A,B)
    else
        C = zeros(n, n)
        m = div(n,  2)
        u = 1:m 
        v = m+1:n 
        S₁ = (A[v, u]+A[v, v])
        S₂ = (S₁ - A[u, u])
        S₃ = (A[u, u]-A[v, u])
        S₄ = (A[u, v] - S₂)
        S₅ = (B[u, v]-B[u, u])
        S₆ = (B[v, v] - S₅)
        S₇ = (B[v, v]-B[u, v])
        S₈ = (S₆-B[v, u])

        tasks = []
        push!(tasks, @spawn multy_wine_strass_mat_mul(S₂, S₆, n_min))
        push!(tasks, @spawn multy_wine_strass_mat_mul(A[u, u], B[u, u], n_min))
        push!(tasks, @spawn multy_wine_strass_mat_mul(A[u, v], B[v, u], n_min))
        push!(tasks, @spawn multy_wine_strass_mat_mul(S₃, S₇, n_min))
        push!(tasks, @spawn multy_wine_strass_mat_mul(S₁, S₅, n_min))
        push!(tasks, @spawn multy_wine_strass_mat_mul(S₄, B[v, v], n_min))
        push!(tasks, @spawn multy_wine_strass_mat_mul(A[v, v], S₈, n_min))
        P₁ ,P₂,P₃,P₄,P₅, P₆,P₇ = fetch.(tasks)

        T₁ = P₁+P₂
        T₂ = T₁+P₄
        C[u, u] = P₂+P₃
        C[u, v] = T₁+P₅+P₆
        C[v, u] = T₂-P₇
        C[v, v] = T₂+P₅
        return C
    end
end
begin
    x = []
    
    norm = []
    wine = []
    straus = []
    multy_straus = []
    wine_straus = []
    multy_wine_straus = []
    for i in 1:10
        push!(x, 2^i)
        n = 2^i
        A = gen_rand(n, 10, -10)
        B = gen_rand(n, 10, -10)
        start = time()
        norm_res = @time normal_mat_mul(A, B)
        push!(norm, time()-start)

        start = time()
        wine_res = @time wine_mut_mul(A, B)
        push!(wine, time()-start)

        start = time()
        strass_res = @time strass_mat_mul(A, B, 8)
        push!(straus, time()-start)
        
        start = time()
        wine_straus_res = @time wine_strass_mat_mul(A, B, 8)
        push!(wine_straus, time()-start)

        start = time()
        multy_res = @time multy_strass_mat_mul(A, B, 8)
        push!(multy_straus, time()-start)


        start = time()
        multy_wine_strass_res = @time multy_wine_strass_mat_mul(A, B, 8)
        push!(multy_wine_straus, time()-start)
    end
    savefig(plot(x, [norm, wine, straus, multy_straus, wine_straus, multy_wine_straus], label=["Normal" "Wine" "Straus" "Straus Multy" "Wine-Straus" "Wine-Straus Multy"]), "myplot.pdf")
end