
# noen contiglengder 

A <- c(156, 105,141,116,95,173,86,128)
B <- c(69,31,230,47,194,344,21,154) 
C <- c(7,241,40,16,594,3,1,98) 
D <- c(62,505,8,15,125,31,251,4) 

# lag kode som beregner N50 

N50 <- function(contig.lengder){
  x.sortert <- sort(contig.lengder, decreasing = T) 
  c.sum <- cumsum(x.sortert) 
  contig.nummer <- min(which(c.sum >= max(c.sum)*0.5))
  N50 <- x.sortert[contig.nummer]
  return(N50) 
}

# lag kode som beregner N25
N25 <- function(contig.lengder){
  x.sortert <- sort(contig.lengder, decreasing = T) 
  c.sum <- cumsum(x.sortert) 
  contig.nummer <- min(which(c.sum >= max(c.sum)*0.25))
  N25 <- x.sortert[contig.nummer]
  return(N25) 
}

# lag kode som beregner N75
N75 <- function(contig.lengder){
  x.sortert <- sort(contig.lengder, decreasing = T) 
  c.sum <- cumsum(x.sortert) 
  contig.nummer <- min(which(c.sum >= max(c.sum)*0.75))
  N75 <- x.sortert[contig.nummer]
  return(N75) 
}

