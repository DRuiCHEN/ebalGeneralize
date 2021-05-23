#' Compute the (extended) entropy balancing weights proposed in the paper.
#' 
#' @param x covariate matrix, arranged as (H, G).
#' @param trt logical vector of treatment assignment.
#' @param target_moments mean of first K_h covariates of the target sample.
#' @param H_add_intercept whether to include 1 in H, default as TRUE.
ebal_generalize <- function(x, trt, 
                            target_moments = NULL, 
                            H_add_intercept = TRUE) 
{
   stopifnot(length(target_moments) <= NCOL(x))
   
   n <- NROW(x)
   
   # If H_add_intercept is TRUE, include 1 in H, so that the sums of the weights
   # in the treated and control group will be kept as 1 in the computation
   if (H_add_intercept) 
   {
      x <- cbind(1, x)
      target_moments <- c(1, target_moments)
   }
   
   # Divide the covariates into the H part and G part. The mean of H part of the
   # weighted treated/control group will be equalized to target_moments. The
   # means of the G part of the treated and control groups will be equalized to
   # each other.
   K_h <- length(target_moments)
   K_g <- NCOL(x) - K_h
   H1 <- x[trt == 1, (0:K_h)[-1], drop = FALSE]
   H0 <- x[trt != 1, (0:K_h)[-1], drop = FALSE]
   G1 <- x[trt == 1, K_h + (0:K_g)[-1], drop = FALSE]
   G0 <- x[trt != 1, K_h + (0:K_g)[-1], drop = FALSE]

   # An internal function to compute the primal and dual arguments from theta
   primal_dual_args <- function(theta) {
      lambda1 <- theta[(0:K_h)[-1]]
      lambda0 <- theta[K_h + (0:K_h)[-1]]
      gamma <- theta[2 * K_h + (0:K_g)[-1]]
      list(lambda1 = lambda1, lambda0 = lambda0, gamma = gamma,
           w1 = drop(exp(H1 %*% lambda1 + G1 %*% gamma - 1)),
           w0 = drop(exp(H0 %*% lambda0 - G0 %*% gamma - 1)))
   }
   
   # Set up the dual problem: functions for computing the dual objective and its
   # gradient
   dual_fn <- function(theta) {
      args <- primal_dual_args(theta)
      (sum(args$w1) + sum(args$w0)) / n - sum((args$lambda1 + args$lambda0) * target_moments)
   }
   
   dual_grad <- function(theta) {
      args <- primal_dual_args(theta)
      c(colSums(args$w1 * H1) / n - target_moments,
        colSums(args$w0 * H0) / n - target_moments,
        colSums(args$w1 * G1) / n - colSums(args$w0 * G0) / n)
   }
   
   # Optimization is done with the built-in optim function
   opt <- optim(rep(0, 2 * K_h + K_g), fn = dual_fn, gr = dual_grad, method = "BFGS")
   
   # Output
   args <- primal_dual_args(opt$par)
   w <- numeric(n)
   w[trt == 1] <- args$w1
   w[trt != 1] <- args$w0
   list(w = w,
        opt_res = opt,
        args = args)
}


#' Simple version of entropy balancing
#'
#' The resulting weights calibrate the whole x sample to the target moments (not
#' distinguishing treated and control), which is equivalent to a expential
#' tilting calibration.
#'
#' @param x covariate matrix.
#' @param target_moments mean of covariates of the target sample. The length of
#'   this vector should equal NCOL(x).
#' @param H_add_intercept whether to include 1 in H, default as TRUE.
ebal_simple <- function(x,
                        target_moments,
                        add_intercept = TRUE)
{
   stopifnot(length(target_moments) == NCOL(x))
   
   # If H_add_intercept is TRUE, include 1 in H, so that the sums of the weights
   # in the treated and control group will be kept as 1 in the computation
   if (add_intercept) 
   {
      x <- cbind(1, x)
      target_moments <- c(1, target_moments)
   }
   
   # Function for computing the primal arguments from the dual arguments
   primal_w <- function(lambda) drop(exp(x %*% lambda - 1))
   
   # Set up the dual problem: functions for computing the dual objective and its
   # gradient
   dual_fn <- function(lambda) mean(primal_w(lambda)) - sum(target_moments * lambda)
   
   dual_grad <- function(lambda) colMeans(primal_w(lambda) * x) - target_moments

   # Optimization is done with the built-in optim function
   opt <- optim(rep(0, NCOL(x)), fn = dual_fn, gr = dual_grad, method = "BFGS")
   
   # Output
   list(w = primal_w(opt$par),
        lambda = opt$par,
        opt_res = opt)
}
