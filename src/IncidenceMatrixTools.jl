
module IncidenceMatrixTools

using FromFile

@from "MakeIncidenceMatrix.jl" using MakeIncidenceMatrix
@from "MakeWeightMatrices.jl" using MakeWeightMatrices

export ij_To_k
export k_To_ij
# export makeIncidenceMatrix2D
export makeIncidenceMatrix3D
# export solChecks
export makeGhostVertexMask
export makeGhostEdgeMask
export vertexVolumeWeightsMatrix
export vertexVolumeWeightsInverseMatrix
export edgeLengthMatrix
export edgeLengthInverseMatrix
export setBoundaryEdgesToZero

end