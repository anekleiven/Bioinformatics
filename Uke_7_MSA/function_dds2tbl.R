# BIN210 UKE 7
# laste inn pakker 

library(msa)
library(microseq)

# kode fra R-lab uke 7 

ecoli.ss <- readDNAStringSet(filepath = "data/upstream30_Ecoli.fasta")

ecoli.tbl <- readFasta('data/upstream30_Ecoli.fasta') 

ecoli.msa <- msa(ecoli.ss, method = "Muscle")

ecoli.ss


# funksjon for å konvertere DNAStringSrt til tabell  

dss2tbl <- function(dss){
  tbl <- tibble(Header = names(dss), 
                Sequence = as.character(dss))
  return(tbl)
}

ecoli.tbl2 <- dss2tbl(ecoli.ss) 
ecoli.tbl2


# funksjon for å konvertere tabell til DNAStringSet 
tbl2dss <- function(tbl){
  dss <- DNAStringSet(tbl$Sequence)
  names(dss) <- tbl$Header
  return(dss)
}


# lage msa-objekt med funksjonen msa() og tbl2dss() 
ecoli.ss2 <- tbl2dss(ecoli.tbl2)
ecoli.ss2 

ecoli.msa2 <- ecoli.tbl |> 
  tbl2dss() |> 
  msa(method = "Muscle")

ecoli.msa2
