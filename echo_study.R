library(tidyverse)
library(foreach)
library(doParallel)
registerDoParallel()

source("method.R")
source("utils.R")


# prepare data ------------------------------------------------------------

load("echo_data.rdata")

vars_H <- c("age", "weight", "genderF", "icd_chf1", "icd_afib1", "icd_respfail1", "icd_malignancy1")

dat <- dat %>% 
   mutate_at(vars(starts_with("lab_") & !ends_with("_flag")), function(v) log(1 + v)) %>% 
   mutate_at(all_of(features_cont), function(fea) {fea <- scale(fea); attributes(fea) <- NULL; fea}) 


fml_subgrp <- c(features_cont, features_disc) %>% paste(collapse = " + ") %>% sprintf(" ~ %s", .)
x <- model.matrix(as.formula(fml_subgrp), data = dat)[, -1]
x <- cbind(x[, vars_H], x[, setdiff(colnames(x), vars_H)])
trt <- as.integer(dat$echo == 1)
vars_disc <- setdiff(colnames(x), features_cont)
y <- dat$mort_28_day == 0


# run estimation ----------------------------------------------------------

set.seed(22)

# takes about 80s
system.time(
   tib_err <- 
      foreach(p_src_scale = c(1, 5), .combine = bind_rows) %:%
      foreach(p_trt_case = 0:1, .combine = bind_rows) %do%
      {
         src_coefs <- p_src_scale * c(age = -.3,
                                      weight = 0,
                                      genderF = 0,
                                      icd_chf1 = .3,
                                      icd_afib1 = .4,
                                      icd_respfail1 = .3,
                                      icd_malignancy1 = .4,
                                      vs_temp_first = 0,
                                      vs_map_first = 0,
                                      vs_heart_rate_first = 0,
                                      lab_platelet_first = 0,
                                      lab_po2_first = 0,
                                      lab_lactate_first = 0,
                                      lab_bun_first = 0,
                                      saps = 0,
                                      sofa = 0,
                                      elix_score = 0)
         
         foreach(1:800, .combine = bind_rows) %dopar% 
            {
               # split into source, target, test
               probs <- pnorm(drop(x %*% src_coefs - .5 * p_src_scale))
               probs <- .5 + .8 * (probs - .5)
               s_idx <- sample(1:NROW(x), round(.4 * NROW(x)), prob = probs)
               s <- rep(0, NROW(x))
               s[s_idx] <- 1
               xs <- x[s == 1,]
               trts <- trt[s == 1]
               ys <- y[s == 1]
               xt <- x[s != 1,]
               trtt <- trt[s != 1]
               yt <- y[s != 1]
               
               in_test <- sample(c(rep(TRUE, round(NROW(xt) * .5)), rep(FALSE, NROW(xt) - round(NROW(xt) * .5))))
               x_test <- xt[in_test,]
               trt_test <- trtt[in_test]
               y_test <- yt[in_test]
               xt <- xt[!in_test,]
               trtt <- trtt[!in_test]
               yt <- yt[!in_test]
               
               # introduce more confounding by subsampling trt and ctrl separately
               trt_idx1 <- which(trts == 1)
               trt_idx0 <- which(trts == 0)
               if (p_trt_case == 0) {
                  fz <- rep(0, NROW(xs))
               }
               if (p_trt_case == 1) {
                  z <- xs[, c("saps", "sofa", "elix_score")] #, "vs_temp_first", "vs_map_first", "vs_heart_rate_first")]
                  fz <- drop(.3 * z[, 1] + .4 * z[,2] - .5 * z[,3])
               }
               if (p_trt_case == 2) {
                  z <- xs[, c("saps", "sofa", "elix_score")]
                  fz <- .3 * z[,1]^2 + .2 * z[,2]^2 + .2 * z[,3]^2 + .4 * z[,1] * z[,2] - .3 * z[,1] + .4* z[,2] + .2 * z[,3] - .6
               }
               prob1 <- pnorm(fz[trt_idx1])
               prob0 <- pnorm(-fz[trt_idx0])
               prob1 <- .5 + .8 * (prob1 - .5)
               prob0 <- .5 + .8 * (prob0 - .5)
               sub_idx1 <- sample(trt_idx1, round(.5 * length(trt_idx1)), prob = prob1)
               sub_idx0 <- sample(trt_idx0, round(.5 * length(trt_idx0)), prob = prob0)
               sub_idx <- sort(c(sub_idx1, sub_idx0))
               xs <- xs[sub_idx, ]
               trts <- trts[sub_idx]
               ys <- ys[sub_idx]
               
               # compute the "true" ATET with treatment/outcome information on the target sample
               
              # prob_trt_test <- get_prop_score(x_test, trt_test)$s
              # ATET <- weighted.mean(y_test[trt_test == 1], 1 / prob_trt_test[trt_test == 1]) - 
              #    weighted.mean(y_test[trt_test != 1], 1 / (1 - prob_trt_test[trt_test != 1]))
               
               prob_trt_test <- weightit(trt_test ~ x_test, method = "ebal", estimand = "ATE")$weights
               ATET <- weighted.mean(y_test[trt_test == 1], prob_trt_test[trt_test == 1]) - 
                  weighted.mean(y_test[trt_test != 1], prob_trt_test[trt_test != 1])
               
               # construct weights
               xs_sub <- xs[, vars_H]
               target_moments <- colMeans(xt[, vars_H])
               
               w_src <- ebal_simple(xs_sub, target_moments)$w
               prob_trt <- get_prop_score(xs, trts)$s
               
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
               
               # compute the difference between the weighting estimates and the truth
               colSums((2 * trts - 1) * ys * wts) - ATET
            } %>% 
            mutate(src = p_src_scale,
                   trt = p_trt_case)
      }
)


# summarize the results (as a latex table) --------------------------------

tib_err %>% 
   gather(method, err, -src, -trt) %>% 
   mutate(method = as_factor(method),
          setting = paste(src, trt), .keep = "unused") %>% 
   group_by(setting, method) %>% 
   summarise(bias = mean(err) * 100,
             rmse = sqrt(mean(err^2)) * 100) %>% 
   mutate(nums = sprintf("%.2f & %.2f", bias, rmse), .keep = "unused") %>%
   spread(setting, nums) %>% 
   unite(out, -method, sep = " && ") %>% 
   mutate(out = paste(out, "\\\\")) %>% 
   unite(out, everything(), sep = ' & ') %>% 
   `$`(out) 
