using Random
using Plots

function random_number_array(n)
    return rand(0:9, n)
end


function standard_multiplication(x, y)
    n = length(x)
    result = zeros(Int, 2n)
    for i in n:-1:1
        carry = 0
        for j in n:-1:1
            result[i + j] += x[i] * y[j] + carry
            carry = div(result[i + j], 10)
            result[i + j] = mod(result[i + j], 10)
        end
        result[i] += carry
    end
    return result
end

function karatsuba(x, y)
    if x < 10 || y < 10
        return x * y
    end
    n = max(floor(Int, log10(x) + 1), floor(Int, log10(y) + 1))
    m = div(n, 2)

    high1, low1 = div(x, 10^m), x % 10^m
    high2, low2 = div(y, 10^m), y % 10^m

    z0 = karatsuba(low1, low2)
    z1 = karatsuba(low1 + high1, low2 + high2)
    z2 = karatsuba(high1, high2)

    return z2 * 10^(2m) + (z1 - z2 - z0) * 10^m + z0
end

function vec2int(v)
    res = 0
    for i in 1:length(v)
        res += v[i]*10^(length(v)-i)
    end
    return res
    
end

begin
    x = []
    norm = []
    kara = []
    for i in 1:13
        push!(x, 2^i)
        n = 2^i
        A = random_number_array(n)
        B = random_number_array(n)
        println(vec2int(standard_multiplication(A, B)))
        norm_res = @elapsed standard_multiplication(A, B)
        push!(norm, norm_res)
        a = vec2int(A)
        b = vec2int(B)
        println(a, " ", b, " ", a*b)
        println(karatsuba(a, b))
        kara_res = @elapsed karatsuba(a, b)
        push!(kara, kara_res)
    end
    savefig(plot(x, [norm, kara], label=["Normal" "karatsuba"]), "myplot.pdf")
end