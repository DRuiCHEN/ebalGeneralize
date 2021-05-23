library(tidyverse)
library(foreach)
library(doParallel)
registerDoParallel()

source("method.R")
source("utils.R")
source("simu_make_data.R")

n <- 800
p <- 5
H_vars <- 1:3
n_simu <- 400

setting_src = 1
link_trt = 1

system.time(
   ATET <- foreach(1:200, .combine = "mean") %dopar% {
      with(make_data(5e4, p, setting_src = setting_src), mean(y1[s != 1] - y0[s != 1]))
   }
)

result <- foreach(setting_trt = 1:3, .combine = bind_rows) %:%
   foreach(setting_main = 1:2, .combine = bind_rows) %do% {
      tmp <- foreach(1:n_simu, .combine = bind_rows) %dopar% {
         dat <- make_data(n, p, 
                          setting_src = setting_src,
                          setting_trt = setting_trt, 
                          setting_main)
         with(
            dat,
            {
               xs <- x[s == 1,]
               xt <- x[s != 1,]
               trts <- trt[s == 1]
               ys <- y[s == 1]
               
               xs_sub <- xs[, H_vars]
               target_moments <- colMeans(xt)[H_vars]
               
               w_src <- ebal_simple(xs_sub, target_moments)$w
               prob_trt <- get_prop_score(xs, trts)$s
               prob_trt_sub <- get_prop_score(xs_sub, trts)$s
               
               wts <-
                  tibble(
                     IPW = trts / prob_trt + (1 - trts) / (1 - prob_trt),
                     EBAL = ebal_generalize(xs_sub, trts, target_moments)$w,
                     "IPW + ET" = IPW * w_src,
                     proposed = ebal_generalize(xs, trts, target_moments)$w
                  )
               
               wts <- wts %>%
                  mutate(trts = trts) %>%
                  group_by(trts) %>%
                  mutate_at(vars(-group_cols()), function(z) z / sum(z)) %>%
                  ungroup() %>%
                  select(-trts)
               
               colSums((2 * trts - 1) * ys * wts) - ATET
            }
         )
      }
      tmp %>% mutate(trt = setting_trt, main = setting_main)
   }
