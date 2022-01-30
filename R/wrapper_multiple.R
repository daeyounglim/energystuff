library(energystuff)

simfun <- function(start, end, seednumber=18007, type="size", save.file=TRUE, ncores=1, verbose=FALSE) {
    if (type == "power") {
        load("multipleAccount_nsim_1000_1p2_Tpower_v3.RData")
    } else if (type == "size") {
        load("multipleAccount_nsim_10000_1p2_Talpha_v3.RData")
    } else {
        stop("Unknown type of analysis. Please check the type argument.")
    }
    df <- lapply(start:end, function(i) df[[i]])
    n.sim <- end - start + 1
    # percentage_increment <- 1.2
    nyear.grid <- c(3, 5, 10, 15)
    alpha.grid <- c(0.025, 0.05, 0.1)
    dnames <- list(c("T=3", "T=5", "T=10", "T=15"),
                   c("a=0.025", "a=0.05", "a=0.1"),
                   c(1:n.sim))
    dnames_v2 <- list(c("T=3", "T=5", "T=10", "T=15"),
                   c("a=0.025", "a=0.05", "a=0.1"),
                   c(1:n.sim),
                   c("holm", "hochberg", "hommel", "BH", "BY"))
    accept_reject <- array(0, dim = c(4, 3, n.sim), dimnames = dnames)
    accept_reject_multiple <- array(0, dim = c(4, 3, n.sim, 5), dimnames = dnames_v2)
  accept_reject_indiv <- array(0, dim = c(4, 3, 5, 3), dimnames = list(c("T=3", "T=5", "T=10", "T=15"), c("a=0.025", "a=0.05", "a=0.1"),  c("holm","hochberg","hommel","BH","BY"), c("k=10","k=11","k=12")))

  tdr_indiv <- array(0, dim = c(4, 3, 5, 3), dimnames = list(c("T=3", "T=5", "T=10", "T=15"), c("a=0.025", "a=0.05", "a=0.1"),  c("holm","hochberg","hommel","BH","BY"), c("k=10","k=11","k=12")))

  accept_reject_indiv_all <- array(0, dim = c(4, 3, n.sim, 5, 3), dimnames = list(c("T=3", "T=5", "T=10", "T=15"), c("a=0.025", "a=0.05", "a=0.1"), c(1:n.sim),  c("holm","hochberg","hommel","BH","BY"), c("k=10","k=11","k=12")))

    
    elapsed_time <- array(0, dim = c(4, 3, n.sim))
    pb <- txtProgressBar(style=3)
    set.seed(seednumber)
    icount <- 1
    for (i.sim in 1:n.sim) {
        dff <- df[[i.sim]]
        for (i in 1:length(nyear.grid)) {
            dfff <- dff[[i]]
            for (j in 1:length(alpha.grid)) {
                alpha <- alpha.grid[j]
                elapsed_time[i, j, i.sim] <- system.time(out <- hybridmonitor.multiple(dfff$y, dfff$y_Kp1, method=c("holm", "hochberg", "hommel", "BH", "BY"), control = list(sig_level = alpha, B1=1000, B2=2000, B3=2000), ncores=ncores, verbose=verbose))[3]

                accept_reject[i, j, i.sim] <- out$test
                accept_reject_multiple[i, j, i.sim, ] <- out$test_grand
                accept_reject_indiv_all[i,j,i.sim,,] <- (out$obs_phat_indiv[,10:12] > out$cr_val_indiv[,10:12])
                accept_reject_indiv[i,j,,] <- accept_reject_indiv[i,j,,] + (out$obs_phat_indiv[,10:12] > out$cr_val_indiv[,10:12])
                tdr_indiv[i,j,,] <- tdr_indiv[i,j,,] + as.numeric(out$obs_phat_indiv[,10:12] > out$cr_val_indiv[,10:12]) * as.numeric(out$test_grand)
                setTxtProgressBar(pb, icount / (n.sim*length(nyear.grid)*length(alpha.grid)))
                icount <- icount + 1
            }
        }
    }
    close(pb)
    if (save.file) {
        save(accept_reject, accept_reject_multiple, accept_reject_indiv, elapsed_time, accept_reject_indiv_all, tdr_indiv, file = paste0("multiple_sim_start_",start,"_end_",end,"_seednumber_",seednumber,"_Talpha_v3.RData"))
    } else {
        return(list(accept_reject=accept_reject, accept_reject_multiple = accept_reject_multiple, accept_reject_indiv = accept_reject_indiv))
    }
}
