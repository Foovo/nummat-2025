module dn01

export RedkaMatrika, getindex, setindex!, firstindex, lastindex, sor

"""
    RedkaMatrika je predstavljena z matrikami V, ki hrani neničelne elemente po vrsticah
    in matriko I, ki hrani indekse elementov po vrsticah
"""
mutable struct RedkaMatrika
    V
    I
end

import Base: getindex, setindex!, firstindex, lastindex, *

"""
    v = getindex(rm::RedkaMatrika, i::Int64, j::Int64)
Pridobi element iz redke matrike na mestu (i,j)
"""
function getindex(rm::RedkaMatrika, i::Int64, j::Int64)

    n, m = size(rm.I)

    if i > n || j > n || i < 1 || j < 1
        throw("Indeksi so preveliki")
    end

    #Preveri če obstaja element j na vrstici i v matriki I -> če ja vrni element na istem indeksu v V

    I = rm.I
    V = rm.V

    for k in 1:m
        if I[i, k] == j 
            return V[i, k]
        end
    end

    return 0
end


"""
    setindex!(rm::RedkaMatrika, v::Float64, i::Int64, j::Int64)
Nastavi element v redki matriki na mestu (i,j) na vrednost v
"""
function setindex!(rm::RedkaMatrika, v::Float64, i::Int64, j::Int64)
    n, m = size(rm.I)
    I = rm.I
    V = rm.V

    if i > n || j > n || i < 1 || j < 1
        throw("Indeksi so preveliki")
    end

    if rm[i, j] != 0
        for k in 1:m
            if I[i, k] == j 
                V[i, k] = v
                return
            end
        end
    end

    for k in 1:m
        #we can directly insert element
        if V[i, k] == 0
            V[i,k] = v
            I[i,k] = j 
            return
        end
    end

    #dodaj prazen stolpec ničel
    rm.V = hcat(V, zeros(size(V, 1)))
    rm.I = hcat(I, zeros(size(I, 1)))

    rm.V[i, m+1] = v
    rm.I[i, m+1] = j
    
end

"""
Vrni prvi indeks Redke matrike
"""
function firstindex(rm::RedkaMatrika)
    return 1
end

"""
Vrni zadnji indeks matrike rm -> kar število elementov
"""
function lastindex(rm::RedkaMatrika)
    n, _ = size(rm.I)
    return n*n
end

"""
Vrni zadnji indeks matrike rm po dimenziji i
"""
function lastindex(rm::RedkaMatrika, i::Int64)
    n, _ = size(rm.I)

    if i > 2
        throw("Indeks je prevelik za dimenzije")
    end

    return n
end

"""
    res = *(rm::RedkaMatrika, v::Vector)
Zmnoži matriko rm z vektorjem v in vrni rezultat res
"""
function *(rm::RedkaMatrika, v::Vector)
    res = zeros(size(v))
    
    for i in firstindex(rm):lastindex(rm, 1)
        for j in firstindex(rm):lastindex(rm, 2)
            res[i] += rm[i, j]*v[j]
        end
    end

    return res
end

#w je ponavadi od 0-2, za w=1 dobimo kar gauss seidel
"""
    x = sor_korak(rm::RedkaMatrika, x0::Vector, b::Vector, w::Float64)
Izvedi 1 korak SOR iteracije za redko matriko rm, začetni približek x0 
desno stran vektor b in parameter omega
"""
function sor_korak(rm::RedkaMatrika, x0::Vector, b::Vector, w::Float64)
    x = copy(x0)

    for i = firstindex(rm):lastindex(rm, 1)
        x_sum = 0
        for j = 1:length(x0)
            if i != j
                x_sum += rm[i, j]*x[j]
            end
        end
        x[i] = (1-w)*x[i] + w/(rm[i,i]) * (b[i] - x_sum) 
    end

    return x
end

"""
    x, i = sor(rm::RedkaMatrika, b, x0, omega; tol=1e-8)

Izvedi max 1000 sor korakov za matriko rm, desno stran vektor b, začetni 
približek x0 in omego omega. 
"""
function sor(rm::RedkaMatrika, b, x0, omega; tol=1e-8)

    maxit = 1000

    for i = 1:maxit
        x = sor_korak(rm, x0, b, omega)
        if maximum(abs, rm*x - b) < tol
            return x, i
        end
        x0 = x
    end
    throw("Iteracija ne konvergira po $maxit korakih!")
  end
    
end