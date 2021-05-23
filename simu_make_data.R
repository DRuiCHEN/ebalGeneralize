make_data <- function(n,
                      p = 5,
                      setting_src = 1,
                      setting_trt = 1, 
                      # setting_cate = 1,
                      setting_main = 1,
                      err_sd = 1)
{
   
   x_fn <- function(n) {
      matrix(4 * runif(n * p) - 2, n, p)
   }
   
   prob_src_fn <- function(x) {
      switch(setting_src,
             `1` = plogis(.4 * x[,1] + .3 * x[,2] - .2 * x[,4]),
             `2` = plogis(.3 * x[,1] + .2 * x[,2] * x[,4] - .2 * x[,4]))
   }
   
   prob_trt_fn <- function(x) {
      z <- switch(setting_trt,
                  `1` = (.7 * x[,2] + .5 * x[,3]),
                  `2` = (.35 * x[,2] + .25 * x[,3] + .2 * x[,4] - .4 * x[,5]),
                  `3` = (.35 * x[,2] - .4 * pmax(x[,3], x[,4]) - .4 * x[,5]))
      plogis(z)
      # plogis(2 * z)
   }
   
   eff_cate_fn <- function(x) {
      x[,1] - .6 * x[,2] - .4 * x[,3]
   }
   
   eff_main_fn <- function(x) {
      switch(setting_main,
             `1` = .5 * x[,1] + .3 * x[,2] + .3 * x[,3] - .4 * x[,4] - .7 * x[,5],
             `2` = .5 * x[,1] + .3 * x[,2]^2 + .2 * exp(x[,3] - x[,4] - 1) - .7 * x[,5])
   }
   
   err_fn <- function(n) rnorm(n, 0, err_sd)
   
   .make_data_inner(n, x_fn, err_fn, prob_src_fn, prob_trt_fn, eff_main_fn, eff_cate_fn)
   
}

.make_data_inner <- function(n, x_fn, err_fn, prob_src_fn, prob_trt_fn, eff_main_fn, eff_cate_fn)
{
   x <- x_fn(n)
   prob_src <- prob_src_fn(x)
   prob_trt <- prob_trt_fn(x)
   eff_main <- eff_main_fn(x)
   eff_cate <- eff_cate_fn(x)
   err <- err_fn(n)
   
   s <- runif(n) <= prob_src
   trt <- runif(n) <= prob_trt
   
   y0 <- eff_main - eff_cate / 2
   y1 <- eff_main + eff_cate / 2
   
   y <- s * (trt * y1 + (1 - trt) * y0 + err)
   
   list(x   = x,
        s   = as.integer(s),
        trt = s * trt,
        y   = y,
        prob_src = prob_src,
        prob_trt = prob_trt,
        y0  = y0,
        y1  = y1)
}
