#
#  MakeWeightMatrices.jl
#  IncidenceMatrixTools
#
#  Created by Christopher Revell on 04/06/2024.

module MakeWeightMatrices

using SparseArrays
using LinearAlgebra
using InvertedIndices

# Vertex weights
# Forming a diagonal matrix of volumes around each vertex, 
# divided by 2 at the periphery
function vertexVolumeWeightsMatrix(dims, spacing)
    matW = fill(prod(spacing), dims...)
    for i=1:length(dims)
        selectdim(matW, i, 1) ./= 2.0
        selectdim(matW, i, dims[i]) ./= 2.0
    end
    vecW = reshape(matW, prod(dims))
    W = spdiagm(vecW)
    return W
end

# Vertex weights inverse
# Forming a diagonal matrix of 1/volumes around each vertex, 
# with volumes divided by 2 at the periphery
function vertexVolumeWeightsInverseMatrix(dims, spacing)
    matW = fill(prod(spacing), dims...)
    for i=1:length(dims)
        selectdim(matW, i, 1) ./= 2.0
        selectdim(matW, i, dims[i]) ./= 2.0
    end
    vecW = reshape(matW, prod(dims))
    W⁻¹ = spdiagm(1.0./vecW)
    return W⁻¹
end

# Diagonal matrix of edge lengths
function edgeLengthMatrix(dims, spacing)
    lvec = Float64[]
    for i=1:length(dims)
        nEdgesi = (dims[i]-1)*prod(dims[Not(i)])
        l_i = fill(spacing[i], nEdgesi)
        append!(lvec, l_i)
    end
    L = spdiagm(lvec)
    return L
end

# Diagonal inverse matrix of edge lengths
function edgeLengthInverseMatrix(dims, spacing)
    lvec = Float64[]
    for i=1:length(dims)
        nEdgesi = (dims[i]-1)*prod(dims[Not(i)])
        l_i = fill(spacing[i], nEdgesi)
        append!(lvec, l_i)
    end
    l⁻¹ = spdiagm(1.0./lvec)
    return l⁻¹
end

# Diagonal matrix of areas perpendicular to each edge, 
# meaning the area through which diffusive 
# flux in the direction of a given edge passes
# Factor of 1/2 applied to peripheral edges, assuming peripheral 
# vertices lie on the edge of the solution domain so that the surrounding
# grid points are halved.
function edgePerpendicularAreaMatrix(dims, spacing)  
    Aperpvec = Float64[]
    for i=1:length(dims)
        edgeDims = copy(dims)
        edgeDims[i] -= 1
        nEdgesi = prod(edgeDims)
        AMat = fill(prod(spacing[Not(i)]), edgeDims...)
        for j = [jj for jj in 1:length(edgeDims) if jj!=i]
            selectdim(AMat, j, 1) ./= 2.0
            selectdim(AMat, j, edgeDims[j]) ./= 2.0
        end
        append!(Aperpvec, reshape(AMat, nEdgesi))
    end
    Aperpₑ   = spdiagm(Aperpvec) 
    return Aperpₑ
end

# Ghost point mask is a 1D vector in which the value at component i 
# is true if vertex i in the flattened vector of vertices is an internal vertex
# but false if vertex i is a ghost vertex in the flattened vector of vertices 
# This can be used to exclude ghost points from calculations over the whole state vector
function makeGhostVertexMask(dims)
    ghostMaskVertex = fill(true, dims...)
    for i=1:length(dims)
        selectdim(ghostMaskVertex, i, 1) .= false
        selectdim(ghostMaskVertex, i, dims[i]) .= false
    end
    # ghostMaskVertex = spdiagm(reshape(ghostMaskVertex, prod(dims)))
    return reshape(ghostMaskVertex, prod(dims))
end

# New faster version
function makeGhostEdgeMask(dims)
    dimEdgeCount = Int64[]
    for i=1:length(dims)
        push!(dimEdgeCount, (dims[i]-1)*prod(dims[Not(i)]))
    end
    nEdges  = sum(dimEdgeCount)

    ghostEdgeMaskVec = fill(true, nEdges)

    dimEdgeArrayStrides = (1, (dims[1]-1), (dims[1]-1)*dims[2])
    for kk=1:dims[3]
        for jj=1:dims[2]
            for ii in [1,(dims[1]-1)]
                edgeIndex = 1 + (ii-1)*dimEdgeArrayStrides[1] + (jj-1)*dimEdgeArrayStrides[2] + (kk-1)*dimEdgeArrayStrides[3]
                ghostEdgeMaskVec[edgeIndex] = false
            end
        end
    end  
    dimEdgeArrayStrides = (1, dims[1], dims[1]*(dims[2]-1))
    for kk=1:dims[3]
        for jj in [1,(dims[2]-1)]
            for ii=1:dims[1]
                edgeIndex = dimEdgeCount[1] + 1 + (ii-1)*dimEdgeArrayStrides[1] + (jj-1)*dimEdgeArrayStrides[2] + (kk-1)*dimEdgeArrayStrides[3]
                ghostEdgeMaskVec[edgeIndex] = false
            end
        end
    end  
    dimEdgeArrayStrides = (1, dims[1], dims[1]*dims[2])
    for kk in [1,(dims[3]-1)]
        for jj=1:dims[2]
            for ii=1:dims[1]
                edgeIndex = dimEdgeCount[1] + dimEdgeCount[2] + 1 + (ii-1)*dimEdgeArrayStrides[1] + (jj-1)*dimEdgeArrayStrides[2] + (kk-1)*dimEdgeArrayStrides[3]
                ghostEdgeMaskVec[edgeIndex] = false
            end
        end
    end  
    return ghostEdgeMaskVec
end


# function setBoundaryEdgesToZero(dims)
#     ghostMaskArray = fill(true, dims...)
#     for i=1:length(dims)
#         selectdim(ghostMaskArray1, i, 1) .= false
#         selectdim(ghostMaskArray1, i, dims[i]) .= false
#     end
#     ghostMask = reshape(ghostMaskArray, prod(dims))
#     return ghostMask
# end

export solChecks
export makeGhostVertexMask
export makeGhostEdgeMask
# export makeGhostEdgeMaskNew
export vertexVolumeWeightsMatrix
export vertexVolumeWeightsInverseMatrix
export edgeLengthMatrix
export edgeLengthInverseMatrix
export edgePerpendicularAreaMatrix
# export setBoundaryEdgesToZero

end



# function makeGhostEdgeMask(dims)
#     ghostMaskEdge = Bool[]
#     for i=1:length(dims)
#         dimsVec_i = copy(dims)
#         dimsVec_i[i] = dimsVec_i[i]-1
#         nEdgesi = prod(dimsVec_i)
#         l_i = fill(true, dimsVec_i...)
#         l_iSize = size(l_i)
#         for (j,s) in enumerate(size(l_i))
#             selectdim(l_i, j, 1) .= false
#             selectdim(l_i, j, s) .= false
#         end
#         append!(ghostMaskEdge, reshape(l_i, nEdgesi))
#     end
#     return ghostMaskEdge
# end