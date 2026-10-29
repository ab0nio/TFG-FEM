module NedelecHex
include("Geometry/ReferenceHex.jl")            # sin dependencias
include("Geometry/Mesh.jl")                    # sin dependencias del paquete
include("Math/GaussQuadrature.jl")             # sin dependencias
include("Math/CoordinateTransforms.jl")        # sin dependencias
include("Bases/Monomials.jl")                  # sin dependencias
include("Bases/DOFs.jl")                       # usa los cuatro anteriores
include("Bases/CoefficientGenerator.jl")       # usa DOFs y Monomials
include("Bases/BasisEval.jl")                  # usa Monomials; eval en ref y físico
include("Assembly/GlobalDOFs.jl")              # usa Mesh + ReferenceHex
include("Assembly/Assembly.jl")                # usa GaussQuadrature + BasisEval
include("Benchmark/HalfFilledCavity.jl")       # usa todo lo anterior
end
