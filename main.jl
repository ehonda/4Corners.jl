module 4Corners

# Function definitions

using Plots
using LaTeXStrings
using CoordinateTransformations
using Distances

# Mathematical functions

# Rectangle points for our rectangle with side length √(2) centered at the origin
A = [-1 / √2, -1 / √2]
B = [1 / √2, -1 / √2]
C = [1 / √2, 1 / √2]
D = [-1 / √2, 1 / √2]

Δ_cart(P) = (euclidean(P, A) + euclidean(P, B) + euclidean(P, C) + euclidean(P, D))
Δ_cart(x, y) = Δ_cart([x, y])

Δ_cartfrompolar(r, ϕ) = Δ_cart(CartesianFromPolar()(Polar(r, ϕ)))

Δ_polarfromcart(x, y) =
    let polar = PolarFromCartesian()([x, y])
        Δ_cartfrompolar(polar.r, polar.θ)
    end

function plot_single_cart(f, res=100)
    x = range(-1 / √2, 1 / √2, length=res)
    y = range(-1 / √2, 1 / √2, length=res)

    zf = @. f(x', y)

    plotf = contour(x, y, zf, color=:turbo, fill=true)
    circles(plotf)

    plot(plotf)
end

function circles(plot)
    for r in 0.1:0.1:1
        ϕ = range(-π, π, length=100)
        x = @. r * cos(ϕ)
        y = @. r * sin(ϕ)

        if r == 0
            plot!(plot, x, y, color=:black, linewidth=2, label="circles")
        else
            plot!(plot, x, y, color=:black, linewidth=2, label=false)
        end

        # plot!(x, y, color=:black, linewidth=2, label=false)
    end
end

#diff_contour_cart(Δ_cart, Δ_polarfromcart)
function diff_contour_cart(f₁, f₂, res=100)
    x = range(-1 / √2, 1 / √2, length=res)
    y = range(-1 / √2, 1 / √2, length=res)

    δ(x, y) = f₁(x, y) - f₂(x, y)

    zf₁ = @. f₁(x', y)
    zf₂ = @. f₂(x', y)
    zδ = @. δ(x', y)

    plotf₁ = contour(x, y, zf₁, color=:turbo, fill=true)
    circles(plotf₁)

    plotf₂ = contour(x, y, zf₂, color=:turbo, fill=true)
    circles(plotf₂)

    plotδ = contour(x, y, zδ, color=:turbo, fill=true)

    plot(plotf₁, plotf₂, plotδ, layout=(1, 3), size=(1200, 400))
end

"""
    Δ(r, ϕ)

The original distance function in polar coordinates. Represents the distance of a point P in the plane to the corner
points of the square with side length √(2) centered at the origin. 

# Examples

```jldoctest
julia> Δ(0, 0)
4.0
```

"""
Δ(r, ϕ) = (√(r^2 + √(2) * r * (cos(ϕ) + sin(ϕ)) + 1)
           + √(r^2 + √(2) * r * (-cos(ϕ) + sin(ϕ)) + 1)
           + √(r^2 + √(2) * r * (-cos(ϕ) - sin(ϕ)) + 1)
           + √(r^2 + √(2) * r * (cos(ϕ) - sin(ϕ)) + 1))

"""
    Δ₁(r, ϕ)

Δ rewritten by using the [harmonic addition formula](https://en.wikipedia.org/wiki/List_of_trigonometric_identities#Sine_and_cosine).

# Examples

```jldoctest
julia> Δ₁(0, 0)
4.0
```

"""
Δ₁(r, ϕ) = (√(r^2 + 2 * r * (cos(ϕ - π / 4)) + 1)
            + √(r^2 - 2 * r * (cos(ϕ + π / 4)) + 1)
            + √(r^2 - 2 * r * (cos(ϕ - π / 4)) + 1)
            + √(r^2 + 2 * r * (cos(ϕ + π / 4)) + 1))

∂ϕΔ(r, ϕ) = ((2 * r * sin(-(ϕ - π / 4)) / √(r^2 + 2 * r * (cos(ϕ - π / 4)) + 1))
             + (2 * r * sin(ϕ + π / 4) / √(r^2 - 2 * r * (cos(ϕ + π / 4)) + 1))
             + (2 * r * sin(ϕ - π / 4) / √(r^2 - 2 * r * (cos(ϕ - π / 4)) + 1))
             + (2 * r * sin(-(ϕ + π / 4)) / √(r^2 + 2 * r * (cos(ϕ + π / 4)) + 1)))

∂rΔ(r, ϕ) = (((2r + 2(cos(ϕ - π / 4))) / √(r^2 + 2 * r * (cos(ϕ - π / 4)) + 1))
             + ((2r - 2(cos(ϕ + π / 4))) / √(r^2 - 2 * r * (cos(ϕ + π / 4)) + 1))
             + ((2r - 2(cos(ϕ - π / 4))) / √(r^2 - 2 * r * (cos(ϕ - π / 4)) + 1))
             + ((2r + 2(cos(ϕ + π / 4))) / √(r^2 + 2 * r * (cos(ϕ + π / 4)) + 1)))

gradient_Δ(r, ϕ) = [∂rΔ(r, ϕ), ∂ϕΔ(r, ϕ)]

gradient_Δ_norm(r, ϕ) = euclidean(gradient_Δ(r, ϕ), [0, 0])

δ_digits(i, r, ϕ) =
    let ops = (-1) .^ digits(i, base=2, pad=2)
        √(r^2 + ops[1] * 2r * cos(ϕ + ops[2] * π / 4) + 1)
    end

Δ_digits(r, ϕ) = sum(δ_digits(i, r, ϕ) for i in 0:3)
# This is the same order as in the original definition of Δ₁
# Δ_digits(r, ϕ) = δ_digits(2, r, ϕ) + δ_digits(1, r, ϕ) + δ_digits(3, r, ϕ) + δ_digits(0, r, ϕ)

# Utility functions

function gridlines(plot, rmax=1, labels=true)
    # odds
    for k in [-3 // 4, -1 // 4, 1 // 4, 3 // 4]
        r = range(0, rmax, length=2)
        ϕ = fill(k * π, 2)

        if k == -3 // 4 && labels
            plot!(plot, r, ϕ, color=:black, linewidth=2, label=L"\frac{(2k + 1)\pi}{4}")
        else
            plot!(plot, r, ϕ, color=:black, linewidth=2, label=false)
        end
    end

    # evens
    for k in [-1, -1 // 2, 0, 1 // 2]
        r = range(0, rmax, length=2)
        ϕ = fill(k * π, 2)

        if k == 0 && labels
            plot!(plot, r, ϕ, color=:grey, linewidth=2, label=L"\frac{2k\pi}{4}")
        else
            plot!(plot, r, ϕ, color=:grey, linewidth=2, label=false)
        end
    end
end

function plot_single(f, res=100, levels=15, rmax=1)
    r = range(0, rmax, length=res)
    ϕ = range(-π, π, length=res)

    zf = @. f(r', ϕ)

    plotf = contour(r, ϕ, zf, color=:turbo, fill=true, levels=levels)
    gridlines(plotf, rmax)

    plot(plotf)
end

function diff_contour(f₁, f₂, res=100)
    r = range(0, 1, length=res)
    ϕ = range(-π, π, length=res)

    δ(r, ϕ) = f₁(r, ϕ) - f₂(r, ϕ)

    zf₁ = @. f₁(r', ϕ)
    zf₂ = @. f₂(r', ϕ)
    zδ = @. δ(r', ϕ)

    plotf₁ = contour(r, ϕ, zf₁, color=:turbo, fill=true)
    gridlines(plotf₁)

    plotf₂ = contour(r, ϕ, zf₂, color=:turbo, fill=true)
    gridlines(plotf₂)

    plotδ = contour(r, ϕ, zδ, color=:turbo, fill=true)
    plot(plotf₁, plotf₂, plotδ, layout=(1, 3), size=(1200, 400))
end

function plot_partials_ϕ(rmax_count=4, levels=15, res=100, size=(3000, 3000))
    plots = []
    for rmax in range(0.1, 1, length=rmax_count)
        r = range(0, rmax, length=res)
        ϕ = range(-π, π, length=res)

        zf = @. ∂ϕΔ(r', ϕ)

        plotf = contour(r, ϕ, zf, color=:turbo, fill=true, levels=levels)
        gridlines(plotf, rmax, false)

        push!(plots, plotf)
    end

    #plot(plots..., layout=(1, rmax_count), size=(1200, 400))
    plot(plots..., size=size)
end

# Tests

using Test

# TODO: Improvements to be made
#           - Not as verbose
#           - Better description
#           - Generator function
#           - Use `@testset let`
#           - They all use the same random numbers (see println output)
@testset "test $i" for i in 1:10
    r = rand()
    ϕ = π * (rand() * 2 - 1)
    println("r = $r, ϕ = $ϕ")
    @test Δ(r, ϕ) ≈ Δ₁(r, ϕ)
end


end # module 4Corners
