# tools/grapham_read.r
#
# Copyright (c) Matti Vihola and Ashley Ford 2008-2013
#
# This file is part of Grapham.
# 
# Grapham is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Grapham is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Grapham.  If not, see <http://www.gnu.org/licenses/>.

grapham_read <- function(fname, nthin=1, nmax=100e6, nblock=1e6, start=0, coda=FALSE) {
  
  stopifnot(file.exists(fname))
  
  ifh <- file(fname, "rb")
   
  # Read the header line
  hdr <- readLines(ifh, n=1)
  hdr_ <- strsplit(substr(hdr, 2, nchar(hdr)-1), "\",\"", fixed=TRUE)
  Nvars <- length(hdr_[[1]])
  
  # Make sure that we read complete records, and that thinning
  # is done at equal steps
  nblock <- max(Nvars*nthin, nblock - (nblock %% (Nvars*nthin)))
  
  # Initialise the output variable
  data = matrix(nrow=0,ncol=Nvars)
  colnames(data) <- hdr_[[1]]
  
  # Omit records from the start, if requested
  readBin(ifh, "double", n=Nvars*start)
  
  # Read data block by block
  repeat {
    dd <-readBin(ifh, "double", n=nblock)
    # Check whether there was enough data
    Ndd <- length(dd)
    if (Ndd==0) {
      break 
    } else if ((Ndd %% Nvars) != 0) {
      print("Warning: last record was broken -- data omitted!")
      dd <- dd[seq(1,Ndd-(Ndd %% Nvars))]
    }
    ddm <- matrix(dd, ncol=Nvars, byrow=TRUE);
    data <- rbind(data, ddm[seq(1,nrow(ddm),nthin),drop=FALSE])
    if (nrow(data)>=nmax) {
      data <- data[seq(1,nmax),]
      break
    }
  }
  
  close(ifh)
  
  if (coda) {
    require(coda)
    return(mcmc(data, start=1+start, thin=nthin))
  } else {
    return(as.data.frame(data))
  }
}
