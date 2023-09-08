# Function definitions

using Plots

# Mathematical functions

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

# Utility functions

function diff_contour(f₁, f₂, res=100)
    r = range(0, 1, length=res)
    ϕ = range(-π, π, length=res)

    δ(r, ϕ) = f₁(r, ϕ) - f₂(r, ϕ)

    zf₁ = @. f₁(r', ϕ)
    zf₂ = @. f₂(r', ϕ)
    zδ = @. δ(r', ϕ)

    plotf₁ = contour(r, ϕ, zf₁, color=:turbo, fill=true)
    plotf₂ = contour(r, ϕ, zf₂, color=:turbo, fill=true)
    plotδ = contour(r, ϕ, zδ, color=:turbo, fill=true)
    plot(plotf₁, plotf₂, plotδ, layout=(1, 3), size=(1200, 400))
end
