export importOcTreeMeshRoman, importOcTreeModelRoman

using jInv.Utils


function importOcTreeMeshRoman(meshfile::ASCIIString)
	
	# open file (throws error if file doesn't exist)
    f    = open(meshfile,"r")
	
	# number of cells of underlying tensor mesh along dimension 1, 2, 3
	line = split(readline(f))
    m1 = parse(Int64,line[1])
    m2 = parse(Int64,line[2])
    m3 = parse(Int64,line[3])
	
	# top corner coordinates
	line = split(readline(f))
    x1 = parse(Float64,line[1])
    x2 = parse(Float64,line[2])
    x3 = parse(Float64,line[3])
	
	# cell size
	line = split(readline(f))
    h1 = parse(Float64,line[1])
    h2 = parse(Float64,line[2])
    h3 = parse(Float64,line[3])
	
	line = split(readline(f))
    n = parse(Int64,line[1])
	
	# read rest of file at ones
	lines = readlines(f)
	
	# close file
	close(f)
	
	# check correct number of lines read
	if n != length(lines)
		error("Invalid number of (i,j,k,bsz) lines in file $meshfile.")
	end
	
	# allocate space
	i1  = zeros(Int64, n)
	i2  = zeros(Int64, n)
	i3  = zeros(Int64, n)
	bsz = zeros(Int64, n)
	
	# convert string array to numbers
	for i = 1:n
		line   = split(lines[i])
	   i1[i]  = parse(Int64,line[1])
	   i2[i]  = parse(Int64,line[2])
	   i3[i]  = parse(Int64,line[3])
	   bsz[i] = parse(Int64,line[4])
	end
	
	# Roman's code starts the OcTree at the top corner. Change to bottom
	# corner.
	i3 = m3 + 2 .- i3 - bsz
	x3 = x3 - m3 * h3

  # S   = sortrows([i3 i2 i1 bsz])
  # i1  = S[:,3]
  # i2  = S[:,2]
  # i3  = S[:,1]
  # bsz = S[:,4]
  
  S   = sub2ind( (m1,m2,m3), i1,i2,i3 )
  p   = sortpermFast(S)[1]
  i1  = i1[p]
  i2  = i2[p]
  i3  = i3[p]
  bsz = bsz[p]
	
	# create mesh object
	S = sparse3(i1,i2,i3,bsz,[m1,m2,m3])
 	M = getOcTreeMeshFV(S,[h1,h2,h3];x0=[x1,x2,x3])
	
end

function importOcTreeModelRoman(modelfile::ASCIIString, mesh::OcTreeMesh, T::DataType=Float64)
	
	# open file (throws error if file doesn't exist)
  f = open(modelfile,"r")
	
	# read everything
	s = readlines(f)
	
	# close
	close(f)
	
	m = length(s)
	
	# check if we have the correct number of cell values
	if m != mesh.nc
		error("Incorrect number of cell values")
	end

	# Roman's code starts the OcTree at the top corner. Here, we start with the bottom corner. Therefore, we need to permute the cells values.
	m1,m2,m3     = mesh.n
	i1,i2,i3,bsz = find3(mesh.S)
	i3           = m3 + 2 .- i3 - bsz

	n = nnz(mesh.S)
  
  # S = cell(n)
  # for i=1:n
  #   S[i] = (i3[i],i2[i],i1[i])
  # end
  # p = sortperm(S)
  
  S = sub2ind( (m1,m2,m3), i1,i2,i3 )
  p = sortpermFast(S)[1]
	
	d = split(s[1])
	n = length(d)
	
	# convert to numbers
  if n == 1
    u = Array{T}(m)
    for i = 1:m
      u[p[i]] = parse(T,s[i])
    end
  else
  	u = Array{T}(m,n)
  	for i = 1:m
  		pi = p[i]
  		d  = split(s[i])
  		for j = 1:n
        u[pi,j] = parse(T,d[j])
  		end
  	end
  end
	
	return u

end
