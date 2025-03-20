# Implementation of the Needleman-Wunsch algorithm for global sequence alignment

needleman_wunsch <- function(seq1, seq2, match = 1, mismatch = -1, gap = -1) {
  n <- nchar(seq1) + 1
  m <- nchar(seq2) + 1
  
  # Initialize scoring matrix
  score_matrix <- matrix(0, nrow = n, ncol = m)
  
  # Initialize gap penalties for first row and column
  score_matrix[1, ] <- seq(0, (m - 1) * gap, by = gap)
  score_matrix[, 1] <- seq(0, (n - 1) * gap, by = gap)
  
  # Fill in the scoring matrix
  for (i in 2:n) {
    for (j in 2:m) {
      char1 <- substr(seq1, i - 1, i - 1)
      char2 <- substr(seq2, j - 1, j - 1)
      
      match_score <- ifelse(char1 == char2, match, mismatch)
      diagonal <- score_matrix[i - 1, j - 1] + match_score
      left <- score_matrix[i, j - 1] + gap
      up <- score_matrix[i - 1, j] + gap
      
      score_matrix[i, j] <- max(diagonal, left, up)
    }
  }
  
  print("Scoring Matrix:")
  print(score_matrix)
  
  # Traceback to get alignment
  alignment1 <- ""
  alignment2 <- ""
  i <- n
  j <- m
  
  while (i > 1 || j > 1) {
    if (i > 1 && j > 1 && score_matrix[i, j] == score_matrix[i - 1, j - 1] + 
        ifelse(substr(seq1, i - 1, i - 1) == substr(seq2, j - 1, j - 1), match, mismatch)) {
      alignment1 <- paste0(substr(seq1, i - 1, i - 1), alignment1)
      alignment2 <- paste0(substr(seq2, j - 1, j - 1), alignment2)
      i <- i - 1
      j <- j - 1
    } else if (i > 1 && score_matrix[i, j] == score_matrix[i - 1, j] + gap) {
      alignment1 <- paste0(substr(seq1, i - 1, i - 1), alignment1)
      alignment2 <- paste0("-", alignment2)
      i <- i - 1
    } else {
      alignment1 <- paste0("-", alignment1)
      alignment2 <- paste0(substr(seq2, j - 1, j - 1), alignment2)
      j <- j - 1
    }
  }
  
  cat("Alignment:", "\n")
  cat(alignment1, "\n")
  cat(alignment2, "\n")
}

# Example usage
needleman_wunsch("ATTCT", "ATCT")
needleman_wunsch("AGTACGCA", "TATGC")
needleman_wunsch("ATCT--", "ATTCT-")

