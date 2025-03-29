using Graphs
using dn01

"""
b = desne_strani(G::AbstractGraph, sprem, koordinate)
Poišči desne strani sistema linearnih enačb za eno koordinato vložitve grafa `G`
s fizikalno metodo. Argument `sprem` je vektor vozlišč grafa, ki nimajo
določenih koordinat. Argument `koordinate` vsebuje eno koordinato za vsa
vozlišča grafa. Metoda uporabi le koordinato vozlišč, ki so pritrjena.
Indeksi v vektorju `b` ustrezajo vozliščem v istem vrstnem redu,
kot nastopajo v argumentu `sprem`.
"""
function desne_strani(G::AbstractGraph, sprem, koordinate)
    set = Set(sprem)
    m = length(sprem)
    b = zeros(m)
    for i = 1:m
        v = sprem[i]
        for v2 in neighbors(G, v)
            if !(v2 in set) # dodamo le točke, ki so fiksirane
                b[i] -= koordinate[v2]
            end
        end
    end
    return b
end


using SparseArrays
"""
A = matrika(G::AbstractGraph, sprem)
Poišči matriko sistema linearnih enačb za vložitev grafa `G` s fizikalno metodo.
Argument `sprem` je vektor vozlišč grafa, ki nimajo določenih koordinat.
Indeksi v matriki `A` ustrezajo vozliščem v istem vrstnem redu,
kot nastopajo v argumentu `sprem`.
"""
function matrika(G::AbstractGraph, sprem)
    # preslikava med vozlišči in indeksi v matriki
    v_to_i = Dict([sprem[i] => i for i in eachindex(sprem)])
    m = length(sprem)
    A = RedkaMatrika(zeros(m, 1), zeros(m, 1))
    for i = 1:m
        vertex = sprem[i]
        sosedi = neighbors(G, vertex)
        for vertex2 in sosedi
            if haskey(v_to_i, vertex2)
                j = v_to_i[vertex2]
                A[i, j] = 1.0
            end
        end
        A[i, i] = float(-length(sosedi))
    end
    return A
end


"""
G = krožna_lestev(n)
Ustvari graf krožna lestev z `2n` točkami.
"""
function krožna_lestev(n)
    G = SimpleGraph(2 * n)
    # prvi cikel
    for i = 1:n-1
        add_edge!(G, i, i + 1)
    end
    add_edge!(G, 1, n)
    # drugi cikel
    for i = n+1:2n-1
        add_edge!(G, i, i + 1)
    end
    add_edge!(G, n + 1, 2n)
    # povezave med obema cikloma
    for i = 1:n
        add_edge!(G, i, i + n)
    end
    return G
end

using Logging
"""
x = cg(A, b; atol=1e-10)
Metoda konjugiranih gradientov za reševanje sistema enačb `Ax = b`
s pozitivno definitno matriko `A`. Argument `A` ni nujno matrika, lahko je
tudi drugega tipa, če ima implementirano množenje z vektorjem `b`.
Metoda ne preverja, ali je argument `A` pozitivno definiten.
"""
function cg(A, b; atol=1e-8)
    # za začetni približek vzamemo kar desne strani
    x = copy(b)
    r = b - A * b
    p = r
    res0 = sum(r .* r)
    for i = 1:length(b)
        Ap = A * p
        alpha = res0 / sum(p .* Ap)
        x = x + alpha * p
        r = r - alpha * Ap
        res1 = sum(r .* r)
        if sqrt(res1) < atol
            @info "Metoda KG konvergira po $i korakih."
            break
        end
        p = r + (res1 / res0) * p
        res0 = res1
    end
    return x
end


#Vložitev 
function vlozitev()
    n = 10
    G = krožna_lestev(n)
    sprem = setdiff(vertices(G), 1:5) #[1,2,3....50]
    A = matrika(G, sprem)

    b = desne_strani(G, sprem, rand(5))

    x,i = sor(A, b, rand(length(b)), 0.5)

    return x
end


using Plots
function najdi_najboljsi_omega(rm::RedkaMatrika, b, x0) 
    segments = 20
    omegas = LinRange(0, 2, segments)

    is = []
    ws = []
    best_i = 1000
    best_w = -1
    for k = 1:segments
        w = omegas[k]

        try
            x, i = sor(rm, b, x0, w)
            push!(is, i)
            push!(ws, w)

            if i < best_i
                best_i = i
                best_w = w
            end
        catch e
            push!(is, 1000)
            push!(ws, w)
        end
   
    end

    plot(ws, is, xlabel="Omege", ylabel="Št. iteracij", title="Število iteracij do konvergence glede na omego")

    return best_i, best_w, is, ws
end

V = [3 1 0; 2 1 0; 1 0 0; 2 1 4]
I = [1 2 0; 2 4 0; 3 0 0; 1 2 4]

rm = RedkaMatrika(copy(V), copy(I))

b=[5.0;8;3;20]
x0 = [0.0;0;0;0]

i, w, is, ws = najdi_najboljsi_omega(rm, b, x0) 
plot(ws, is, xlabel="Omege", ylabel="Št. iteracij", title="Število iteracij do konvergence glede na omego", legend=false)
savefig("./plot.png")
#best i: 18
#best w: 0.947368

