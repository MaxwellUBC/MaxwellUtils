import jInv.Mesh.getEdgeIntegralOfPolygonalChain

export getEdgeIntegralOfPolygonalChain

function getEdgeIntegralOfPolygonalChain(
  mesh::AbstractMesh,
  trx::Array{TrxRcv,1};
  normalize=false)
  
  nEx,nEy,nEz = getEdgeNumbering(mesh)
  
  m = sum(mesh.ne)
  n = length(trx)
  
  I = Int64[]
  J = Int64[]
  A = Float64[]
  for i = 1:n
    U = getEdgeIntegralOfPolygonalChain(mesh,trx[i].trxpts',nEx,nEy,nEz,normalize=normalize)
    append!(I, find(U))
    append!(J, fill(i, countnz(U)))
    append!(A, nonzeros(U))
  end
  return sparse(I,J,A,m,n)
  
end
