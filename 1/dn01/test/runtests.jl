
using Test
using dn01

A = [3 1 0 0; 0 2 0 1; 0 0 1 0; 2 1 0 4]

V = [3 1 0; 2 1 0; 1 0 0; 2 1 4]
I = [1 2 0; 2 4 0; 3 0 0; 1 2 4]

rm = RedkaMatrika(copy(V), copy(I))
M = RedkaMatrika(copy(V), copy(I))

@testset "Definiran je tip za razpršene matrike" begin
    @test rm.I == I
end

@testset "Elementi matrike tipa RedkaMatrika so dostopni kot rm(i, j)" begin
    @test rm[1, 1] == A[1,1]
    @test rm[1, 2] == A[1,2]
    @test rm[1, 4] == A[1,4]
    @test rm[4, 1] == A[4,1]
end

@testset "Elementi matrike tipa RedkaMatrika so nastavljivi kot rm(i, j) = v" begin
    M[1,1] = 5.0
    @test M[1, 1] == 5

    M[3, 1] = 1.0
    @test M[3,1] == 1


    @test M[1, 4] == A[1,4]
    print(M)

    
    @test rm[3, 2] == A[3, 2]
end

@testset "Firstindex" begin
    @test firstindex(rm) == 1
end

@testset "Lastindex" begin
    @test lastindex(rm) == 16
    @test lastindex(rm, 1) == 4
    @test lastindex(rm, 2) == 4
end

@testset "Množenje z vektorjem" begin
    v = [1; 2; 3; 4]
    @test rm*v == [5; 8; 3; 20]
end

@testset "sor" begin
    b=[5.0;8;3;20]
    x0 = [0.0;0;0;0]

    x, i = sor(rm, b, x0, 0.5)
    
    @test isapprox(x,[1; 2; 3; 4])
end
