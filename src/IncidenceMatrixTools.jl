#
#  IncidenceMatrixTools.jl
#  IncidenceMatrixTools
#
#  Created by Christopher Revell on 31/07/2024.

module IncidenceMatrixTools

    using FromFile

    @from "MakeIncidenceMatrix.jl" using MakeIncidenceMatrix
    @from "MakeWeightMatrices.jl" using MakeWeightMatrices

    export ij_To_k
    export k_To_ij
    export makeIncidenceMatrix3D
    export vertexVolumeWeightsMatrix
    export vertexVolumeWeightsInverseMatrix
    export edgeLengthMatrix
    export edgeLengthInverseMatrix
    export edgePerpendicularAreaMatrix
    export makeGhostVertexMask
    export makeGhostEdgeMask

end