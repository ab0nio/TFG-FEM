module NedelecHex
include("Geometry/ReferenceHex.jl")        # sin dependencias
include("Math/GaussQuadrature.jl")         # sin dependencias
include("Math/CoordinateTransforms.jl")    # sin dependencias
include("Bases/Monomials.jl")              # sin dependencias
include("Bases/DOFs.jl")                   # usa los cuatro anteriores
include("Bases/CoefficientGenerator.jl")   # usa DOFs y Monomials
end
