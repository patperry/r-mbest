unpivotRp <- function(qr,nvars){
  rank <- qr$rank
  Rp <- matrix(0, rank, nvars)
  if (rank > 0L) {
    R <- qr.R(qr)
    pivot <- qr$pivot
    Rp[,pivot] <- R[seq_len(rank),]
  }
  return(Rp)
}
