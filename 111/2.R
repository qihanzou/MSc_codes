bures_distance <- function(A, B) {
  A_half <- sqrtm(A)
  inner <- A_half %*% B %*% A_half
  inner_sqrt <- sqrtm(inner)
  trace_A <- sum(diag(A))
  trace_B <- sum(diag(B))
  trace_middle <- sum(diag(inner_sqrt))
  val <- trace_A + trace_B - 2 * trace_middle
  # For Numerical stable
  val <- max(val, 0)
  return(sqrt(val))
}

find_T_matrix <- function(Sigma_x, Sigma_z) {
  S_z_sqrt <- sqrtm(Sigma_z)
  middle <- S_z_sqrt %*% Sigma_x %*% S_z_sqrt
  middle_sqrt_inv <- solve(sqrtm(middle))
  Tm <- S_z_sqrt %*% middle_sqrt_inv %*% S_z_sqrt
  return(Tm)
}

is_positive_definite <- function(A) {
  if (!isSymmetric(A)) return(FALSE)
  result <- tryCatch({
    chol(A)
    TRUE
  }, error = function(e) FALSE)
  return(result)
}

