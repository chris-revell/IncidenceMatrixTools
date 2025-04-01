#
#  MakeIncidenceMatrix.jl
#  IncidenceMatrixTools
#
#  Created by Christopher Revell on 02/05/2024.

module MakeIncidenceMatrix

using LinearAlgebra
using SparseArrays
using Base.Threads

# State matrix uMat has dimensions Nplus x Nplus. 
# Column major order, so when flattened to a state vector u, u[1:Nplus] == uMat[:,1], u[1+Nplus:2*Nplus]==uMat[2,:] etc
# State vector has size Nplus*Nplus
# Let us say that in the state matrix, the first dimension is ν and the second is x.
# So in the state vector, the first Nplus components correspond to 0<=ν<=1.0 for x=0.0; second Nplus components correspond to 0<=ν<=1.0 for x=dx, and so on 

# Function to convert index (i,j) to access uMat[i,j] where uMat is of size N1xN2 to index k, 
# which is the corresponding index of the same data point when uMat is flattened to a vector 
ij_To_k(i, j, N1) = (j-1)*N1+i 
# Function to convert index k to access u[k] where u is a vector of size N1*N2 to index (i,j) to access the same data 
# in uMat, where uMat is u reshaped to be a 2D matrix of size N1xN2 
function k_To_ij(k, N1)
    j = floor(Int64, (k-1)/N1)
    i = k-j*N1
    return (i,j)
end

# A is a matrix of nVerts x nEdges
# A[1:nEdgesi, :] corresponds to i directed edges
# A[1:nEdgesj, :] corresponds to i directed edges
# A[1:nEdgesk, :] corresponds to i directed edges
# Vertices and edges in each direction ordered in 1D vector as per julia system data layout
function makeIncidenceMatrix3D(dims)

    if length(dims) == 3
        ni, nj, nk = dims
    elseif length(dims) == 2
        ni, nj= dims
        nk = 1
    end

    nVerts = ni*nj*nk
    nEdgesi = (ni-1)*nj*nk
    nEdgesj = ni*(nj-1)*nk
    nEdgesk = ni*nj*(nk-1)
    nEdges = nEdgesi+nEdgesj+nEdgesk
    is = Int64[]
    js = Int64[]
    ks = Int64[]
    vertArrayStrides = (1, ni, ni*nj)

    # Link i-directed edges to corresponding vertices
    dimEdgeArrayStrides = (1, (ni-1), (ni-1)*nj)
    for kk=1:nk
        for jj=1:nj
            for ii=1:(ni-1)
                edgeIndex = 1 + (ii-1)*dimEdgeArrayStrides[1] + (jj-1)*dimEdgeArrayStrides[2] + (kk-1)*dimEdgeArrayStrides[3]
                v1 = 1 + (ii-1)*vertArrayStrides[1] + (jj-1)*vertArrayStrides[2] + (kk-1)*vertArrayStrides[3]
                v2 = 1 + (ii)*vertArrayStrides[1] + (jj-1)*vertArrayStrides[2] + (kk-1)*vertArrayStrides[3]
                # A[edgeIndex,v1] = -1
                push!(is, edgeIndex)
                push!(js, v1)
                push!(ks, -1)
                # A[edgeIndex,v2] = 1
                push!(is, edgeIndex)
                push!(js, v2)
                push!(ks, 1)
            end
        end
    end  
    # Link j-directed edges to corresponding vertices
    dimEdgeArrayStrides = (1, ni, ni*(nj-1))
    for kk=1:nk
        for jj=1:(nj-1)
            for ii=1:ni
                edgeIndex = nEdgesi + 1 + (ii-1)*dimEdgeArrayStrides[1] + (jj-1)*dimEdgeArrayStrides[2] + (kk-1)*dimEdgeArrayStrides[3]
                v1 = 1 + (ii-1)*vertArrayStrides[1] + (jj-1)*vertArrayStrides[2] + (kk-1)*vertArrayStrides[3]
                v2 = 1 + (ii-1)*vertArrayStrides[1] + (jj)*vertArrayStrides[2] + (kk-1)*vertArrayStrides[3]
                # A[edgeIndex,v1] = -1
                push!(is, edgeIndex)
                push!(js, v1)
                push!(ks, -1)
                # A[edgeIndex,v2] = 1
                push!(is, edgeIndex)
                push!(js, v2)
                push!(ks, 1)
            end
        end
    end  
    # Link k-directed edges to corresponding vertices
    dimEdgeArrayStrides = (1, ni, ni*nj)
    for kk=1:(nk-1)
        for jj=1:nj
            for ii=1:ni
                edgeIndex = nEdgesi + nEdgesj + 1 + (ii-1)*dimEdgeArrayStrides[1] + (jj-1)*dimEdgeArrayStrides[2] + (kk-1)*dimEdgeArrayStrides[3]
                v1 = 1 + (ii-1)*vertArrayStrides[1] + (jj-1)*vertArrayStrides[2] + (kk-1)*vertArrayStrides[3]
                v2 = 1 + (ii-1)*vertArrayStrides[1] + (jj-1)*vertArrayStrides[2] + (kk)*vertArrayStrides[3]
                # A[edgeIndex,v1] = -1
                push!(is, edgeIndex)
                push!(js, v1)
                push!(ks, -1)
                # A[edgeIndex,v2] = 1
                push!(is, edgeIndex)
                push!(js, v2)
                push!(ks, 1)
            end
        end
    end  

    A = sparse(is, js, ks)

    return A
end

export ij_To_k
export k_To_ij
export makeIncidenceMatrix3D

end
