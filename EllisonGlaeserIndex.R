

weighted_EG <- function(s, x, Psi, HI) {
  
  # --- Input checks ---
  s <- as.numeric(s)
  x <- as.numeric(x)
  
  if(length(s) != length(x)) {
    stop("s and x must have the same length")
  }
  
  if(!all(dim(Psi) == c(length(s), length(s)))) {
    stop("Psi must be a square matrix with dimensions length(s) x length(s)")
  }
  
  if(HI < 0 || HI > 1) {
    stop("HI must be between 0 and 1")
  }
  
  # --- Compute GS ---
  d <- s - x
  GS <- as.numeric(t(d) %*% Psi %*% d)
  
  # --- Compute x' Psi x ---
  xPx <- as.numeric(t(x) %*% Psi %*% x)
  
  # --- Weighted Ellison-Glaeser index ---
  gamma_S <- (GS - HI * (1 - xPx)) /
    ((1 - HI) * (1 - xPx))
  
  return(list(
    gamma_S = gamma_S,
    GS = GS,
    xPsiX = xPx
  ))
}


# Regional employment shares for an industry
s <- c(0.40, 0.30, 0.20, 0.10)

# Aggregate employment shares
x <- c(0.25, 0.25, 0.25, 0.25)

# Spatial weights matrix: Psi = I + W
W <- matrix(c(
  0,1,0,0,
  1,0,0,0,
  0,0,1,0,
  0,1,0,1
), nrow = 4, byrow = TRUE)

# Row-standardize W
W <- W / rowSums(W)

# Replace NaN rows if isolated regions exist
W[is.na(W)] <- 0

Psi <- diag(4) + W

# Herfindahl index of plant sizes
HI <- 0.18 #

weighted_EG(s, x, Psi, HI)


nfirm = 4
si   <- c(0.45, 0.3, 0.20, 0.05)

smax <- rep(1/nfirm, nfirm)
Wmax <- matrix(rep(1/(nfirm-1), nfirm^2), nrow = nfirm, byrow = TRUE)
diag(Wmax)<- 0

t(smax)%*%Wmax%*%smax

W <- matrix(c(
  0,1,0,1,
  1,0,0,0,
  1,0,0,1,
  0,1,1,0
), nrow = 4, byrow = TRUE)

W <- matrix(c(
  0,1,0,0,
  0,0,0,0,
  0,0,0,0,
  0,0,0,0
), nrow = 4, byrow = TRUE)

# Row-standardize W
W <- W / rowSums(W)
W[is.na(W)] <- 0

# W <-  W - maxW  #maxW W/(nfirm-1)
# Replace NaN rows if isolated regions exist
diag(W) <- 0

W
HI <- sum(si^2)

HI - (t(si)%*%W%*%si)*(1/nfirm)



