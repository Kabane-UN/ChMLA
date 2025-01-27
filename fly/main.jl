using LinearAlgebra
function LU(A)
    n = size(A, 1)
    L = zeros(n, n)
    U = zeros(n, n)
    for i in 1:n 
        L[i, i] = 1.
    end
    for i in 1:n 
        for j in 1:n
            s = sum([L[i, k]*U[k, j] for k in 1:i])
            if i <= j               
                U[i, j] = A[i, j] - s
            elseif  i > j 
                L[i, j] = (A[i, j]-s)/U[j, j]
            end
        end
    end
    return L, U 
end
begin
    

    A = [1 0 0 ; 0 1 0; 0 0 1;]
    L, U = LU(A)
    println(L, U)
    de = 1.0
    for i in 1:size(A, 1)
        global de = de*U[i, i]
    end
    println("My: ", de)
    println("Lib: ", det(A))
end