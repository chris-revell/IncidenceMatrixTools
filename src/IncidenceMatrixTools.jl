
module IncidenceMatrixTools

    using FromFile

    @from "MakeIncidenceMatrix.jl" using MakeIncidenceMatrix
    @from "MakeWeightMatrices.jl" using MakeWeightMatrices

    export ij_To_k
    export k_To_ij
    export makeIncidenceMatrix3D
    export makeGhostVertexMask
    export makeGhostEdgeMask
    export makeGhostEdgeMaskNew
    export vertexVolumeWeightsMatrix
    export vertexVolumeWeightsInverseMatrix
    export edgeLengthMatrix
    export edgeLengthInverseMatrix
    export setBoundaryEdgesToZero

end