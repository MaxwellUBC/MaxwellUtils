export writeAllFiles, writeTrxRcvFile, writeFrqFile, writeDataFile

function writeAllFiles( datafile::ASCIIString,
                        trxfile::ASCIIString,
                        rcvfile::ASCIIString,
                        frqfile::ASCIIString,
												data::Array{datainfo},
												trx::Array{TrxRcv},
												rcv::Array{TrxRcv},
												frq::Array{freqinfo},
                        only_loc::Bool )  # true to only write locations, false for data and sd

writeTrxRcvFile(trxfile,trx)
writeTrxRcvFile(rcvfile,rcv)
writeFrqFile(frqfile,frq)
writeDataFile(datafile, data, only_loc)

end

function writeTrxRcvFile( trxfile::ASCIIString, trx::Array{TrxRcv} )
# Read the transmitter or receiver information into
# an array of type TrxRcv.

   f = open(trxfile,"w")
	 
	 ntrx = length(trx)
	 for itrx = 1:ntrx
		 idx     = trx[itrx].idx
		 trxtype = trx[itrx].trxtype
		 trxpts  = trx[itrx].trxpts
		 npts    = size(trxpts, 2)
		 @printf(f, "%d %d %d\n", idx, npts, trxtype)
		 for ipts = 1:npts
			 @printf(f, "%.15g %.15g %.15g\n", trxpts[1,ipts], trxpts[2,ipts], trxpts[3,ipts])
		 end
		 @printf(f, "\n")
	 end
	 
	 close(f)
	 
end

function writeFrqFile( frqfile::ASCIIString, frq::Array{freqinfo})

	f = open(frqfile,"w")

	nfrq = length(frq)
	for ifrq = 1:nfrq
		idx       = frq[ifrq].idx
		frequency = frq[ifrq].omega / (2 * pi)
		@printf(f, "%d %.15g\n", idx, frequency)
	end

	close(f)

end

function writeDataFile( datafile::ASCIIString, data::Array{datainfo}, only_loc::Bool )

	f = open(datafile,"w")
	
	ndata = length(data)
	for idata = 1:ndata
		
		trx_idx = data[idata].trx_idx
		frq_idx = data[idata].frq_idx
		rcv_idx = data[idata].rcv_idx
		dataid  = data[idata].dataid
		@printf(f, "%d %d %d %d", trx_idx, frq_idx, rcv_idx, dataid)
		
		if !only_loc
			
			ncomp = 2
			
			dobs  = data[idata].dobs
			sd    = data[idata].sd
			
			for icomp = 1:ncomp
				@printf(f, " %.15g %.15g", dobs[icomp], sd[icomp])
			end
			
		end
		
		@printf(f, "\n")
		
	end

	close(f)
	
end
