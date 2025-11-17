#' Computing the product of two matrices with symbolic coefficients
#'
#' This function computes the product of two squares matrices of same size with character strings as coefficients and simplifies product with "0" or "1"
#'
#' @param A Matrix of character strings
#' @param B Matrix of character strings
#' @return The product AB cleaned of product with "0" and "1"
matrix_product_lit <- function(A, B) {
  n_row <- nrow(A)
  n_col <- ncol(B)
  n_common <- ncol(A)

  result <- matrix("", nrow = n_row, ncol = n_col)


  for (i in 1:n_row) {
    for (j in 1:n_col) {
      terms <- c()
      #Computation of each term of the product matrix AB
      for (k in 1:n_common) {
        a <- A[i, k]
        b <- B[k, j]

        # Step 1 : Product of terms
        prod <- ""
        # if one of the two terms is "0", the product is "0"
        if (a == "0" || b == "0") {
          prod <- "0"
          # if one of the two terms is "1", the product is the other term
        } else if (a == "1") {
          prod <- b
        } else if (b == "1") {
          prod <- a
          # in all other cases, the product is the concatenation of the two terms with "*" between them
        } else {
          prod <- paste0(a, " * ", b)
        }

        terms <- c(terms, prod)
      }

      # Step 2 : Sum of all the products computed before
      # All zeros are removed and the sum is the concatenation of remaining terms with "+" between each of them
      non_zero_terms <- terms[terms != "0"]
      if (length(non_zero_terms) == 0) {
        result[i, j] <- "0"
      } else {
        result[i, j] <- paste(non_zero_terms, collapse = " + ")
      }
    }
  }

  return(result)
}

#' General constructor for projection matrices (pre/post, numeric/literal)
#'
#' Builds a projection matrix according to breeding timing and mode (numeric or symbolic).
#'
#' @param n Number of classes.
#' @param t List of transition rates, format: list(c(i, j, value)) — `value` is either a float or a character string.
#' @param r List of reproduction rates, format: list(c(i, j, value)) — `value` is either a float or a character string.
#' @param mode "num" for numeric mode, "lit" for symbolic/literal mode.
#' @param timing "pre" or "post" to select breeding timing.
#' @return If mode is "num": list(matrix, growth_rate, stable_distribution).
#'         If mode is "lit": symbolic matrix (character matrix).
#' @export
#' @examples
#' t <- list(c(2,1, 0.5), c(3,2, 0.6))
#' r <- list(c(1,3, 1.2))
#' projection_matrix(3, t, r, mode = "num", timing = "pre")
#'
#' t <- list(c(2,1, 0.5), c(3,2, 0.6))
#' r <- list(c(1,3, 1.2))
#' projection_matrix(3, t, r, mode = "num", timing = "post")
#'
#' t <- list(c(2,1, "t12"), c(3,2, "t23"))
#' r <- list(c(1,3, "r31"))
#' projection_matrix(3, t, r, mode = "lit", timing = "post")
#'
#' t <- list(c(2,1, "t12"), c(3,2, "t23"))
#' r <- list(c(1,3, "r31"))
#' projection_matrix(3, t, r, mode = "lit", timing = "pre")
projection_matrix <- function(n, t, r, mode = "num", timing = "pre") {

  # Matrices T et R
  T <- matrix(0, nrow = n, ncol = n)
  R <- diag(n)

  # Convert identity into "1" if symbolic option
  if (mode == "lit") mode(R) <- "character"

  # Filling T
  for (trans in t) {
    i <- as.integer(trans[1])
    j <- as.integer(trans[2])
    T[j, i] <- trans[3]
  }

  # Filling R
  for (repro in r) {
    i <- as.integer(repro[1])
    j <- as.integer(repro[2])
    if (mode == "num") {
      R[j, i] <- R[j, i] + as.numeric(repro[3])
    } else {
      # literal option : concatenation with "+" if value other than "0"
      if (R[j, i] == "1") {
        R[j, i] <- paste0("(", R[j, i], " + ", repro[3], ")")
      } else if (R[j, i] == "0") {
        R[j, i] <- repro[3]
      } else {
        R[j, i] <- paste0("(", R[j, i], " + ", repro[3], ")")
      }
    }
  }

  # Product of matrices
  if (mode == "num") {
    if (timing == "pre") {
      M <- T %*% R
    } else {
      M <- R %*% T
    }
    eig <- eigen(M)
    lambda <- max(Re(eig$values))
    vect <- Re(eig$vectors[, which.max(Re(eig$values))])
    return(list(
      matrix = M,
      growth_rate = lambda,
      stable_distribution = vect / sum(vect)
    ))

  } else {
    # literal option → use function matrix_product_lit
    M <- if (timing == "pre") matrix_product_lit(T, R) else matrix_product_lit(R, T)
    return(M)
  }
}
