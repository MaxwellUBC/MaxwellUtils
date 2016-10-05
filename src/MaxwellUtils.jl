module MaxwellUtils

using jInv.Utils
using jInv.Mesh
using jInv.ForwardShare
hasJOcTree = false
try
  using JOcTree
  hasJOcTree = true
catch
  hasJOcTree = false
end

# ----------------------------------------------------------
# Old UBC EM frequency domain data format

type Transmitter
   trxpts      # (3,npts) points making up the transmitter
   omega       # 2*pi*frequency
   rcvpts      # (3,nrcv) receiver locations
   # The following is used only when only_loc=false
   data        # (12,nrcv)
   sd          # (12,nrcv) standard deviation
   data_exist  # (12,nrcv) true if there is data
   ndata       # number of data
end # type Transmitter
export Transmitter


# The 12 data columns are ordered:
#  1   2    3   4    5   6    7   8    9   10   11  12
# Exr Exi  Eyr Eyi  Ezr Ezi  Hxr Hxi  Hyr Hyi  Hzr Hzi

# ---------------------------------------------------------

# ---------------------------------------------------------
# New UBC EM frequency domain data format

type TrxRcv
   idx::Int         # unique integer index value
   trxtype::Int     # type of transmitter or receiver
   trxpts::Array    # (3,npts) points making up the transmitter or receiver
end # type TrxRcv

type freqinfo
   idx::Int          # unique integer frequency index
   omega::Float64    # 2*pi*frequency
end  # type freqinfo

type datainfo
   trx_idx::Int
   frq_idx::Int
  # omega  # 2*pi*frequency
   rcv_idx::Int
   dataid::Int
   dobs::Array{Float64}   # (2)  observed data
   sd::Array{Float64}     # (2)  standard deviation
end  # type datainfo

export TrxRcv, freqinfo, datainfo

# ---------------------------------------------------------

include("readUBCData.jl")
include("writeUBCData.jl")
include("readDataFiles.jl")
include("writeDataFiles.jl")
include("getTrxOmega.jl")
include("getSxRxFromData.jl")
include("getInitialModel.jl")
include("getDobsWdFromTrx.jl")
include("arrayUtils.jl")
include("diagM.jl")
include("QuickHull.jl")
include("getEdgeIntegralOfPolygonalChain.jl")

if hasJOcTree
  include("setupMeshParam.jl")
  include("readTopo.jl")
  include("createSmallMeshFromTX.jl")
  include("exportMesh.jl")
  include("importMesh.jl")
  include("setupBigOctreeMeshPolygon.jl")
  include("prepareMesh2MeshOT.jl") 
  include("createOcTreeFromTRX.jl")
end

end
