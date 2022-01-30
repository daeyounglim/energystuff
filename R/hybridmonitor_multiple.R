#' Run hybrid monitoring procedure for multiple accounts
#' 
#' This is a function that run hybrid monitoring algorithm for multiple accounts
#' 
#' @author Daeyoung Lim, \email{daeyoung.lim@uconn.edu}
#' @param y an array (12 x K x T) of measured energy usage to be used as the response variable.
#' @param y_Kp1 a matrix of the last cycle of energy usage to be monitored.
#' @param m_idx (Optional) a vector of integers (1-12) that indicate the months to be tested. Note that 1 does not represent January. This completely depends on the user's definition of a *cycle*. If not provided, it will by default run for all 12 months (i.e., full cycle).
#' @param method (Optional) a vector of characters ("bonferroni", "holm", "hochberg", "hommel", "BH", "BY") indicating the selection of multiplicity corrections. If not provided, the multiplicity correction will default to Bonferroni's method.
#' @param control (Optional) a list of parameters for the [Monte Carlo sampling](https://en.wikipedia.org/wiki/Monte_Carlo_method) of the posterior predictive distribution of *p*-values, as well as the *p*-value of the test data.
#' @param ncores (Optional)
#' @return a dataframe with the result of the hypothesis testing, empirical distribution of the *p*-values, observed *p*-value
#' @examples
#' \dontrun{
#' fit <- hybridmonitor.multiple(y, y_Kp1)
#' }
#' @importFrom parallel detectCores
#' @seealso \code{\link{hybridmonitor.single}} for single-account monitoring
#' @md
#' @export
'hybridmonitor.multiple' <- function(y, y_Kp1, m_idx=NULL, method, control=list(), ncores=NULL, verbose=FALSE) {
	if (missing(method)) {
		stop("Method argument should be specified.\nPlease select from p.adjust.methods.")
	}
	cntrl <- list(B1=2000, B2=2000, B3 = 5000, sig_level=0.05)
	cntrl[names(control)] <- control
	K <- ncol(y_Kp1)
	z_alpha_K <- qnorm((1-cntrl$sig_level)^(1/K))
	if (is.null(m_idx)) {
		m_idx <- c(0L:11L)
	} else {
		if (max(m_idx) > 12) {
			stop("m_idx cannot contain indices bigger than 12.")
		} else if (min(m_idx) < 1) {
			stop("m_idx cannot contain indices smaller than 1.")
		}
		m_idx <- m_idx - 1
	}
	if (!is.null(ncores)) {
		ncores_ <- parallel::detectCores()
		if (ncores > ncores_) {
			## parallel computing with # of cores greater than existing cores
			## undermines performance (worse than single-thread computing)
			stop(paste0("The number of cores must not exceed ", ncores_))
		}
	} else {
		## if not provided, use 2 cores or the number of cores, whichever is smaller
		ncores <- min(2, parallel::detectCores())
	}	
	T <- dim(y)[3]

	methods <- match.arg(tolower(method), c("bonferroni", "holm", "hochberg", "hommel", "bh", "by"), several.ok = TRUE)
	do_bonferroni <- "bonferroni" %in% methods
	do_holm <- "holm" %in% methods
	do_hochberg <- "hochberg" %in% methods
	do_hommel <- "hommel" %in% methods
	do_BH <- "bh" %in% methods
	do_BY <- "by" %in% methods

	rnames <- c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY")[c(do_bonferroni, do_holm, do_hochberg, do_hommel, do_BH, do_BY)]

	out <- .Call(`_energystuff_hybridmonitor_multiple_select_months`,
		as.array(y),
		as.matrix(y_Kp1),
		as.double(cntrl$sig_level),
		as.double(z_alpha_K),
		as.integer(cntrl$B1),
		as.integer(cntrl$B2),
		as.integer(cntrl$B3),
		as.integer(T),
		as.integer(m_idx),
		as.logical(do_bonferroni),
		as.logical(do_holm),
		as.logical(do_hochberg),
		as.logical(do_hommel),
		as.logical(do_BH),
		as.logical(do_BY),
		as.logical(verbose),
		as.integer(ncores)
		)
	out$idx <- out$idx + 1L # C++ zero-indexing -> R indexing
	out$test_grand <- c(out$test_grand)
	names(out$test_grand) <- rnames
	rownames(out$cr_val_grand) <- rnames
	rownames(out$obs_phat_grand) <- rnames
	rownames(out$cr_val_indiv) <- rnames
	rownames(out$obs_phat_indiv) <- rnames
	out$test_indiv <- apply(out$cr_val_indiv < out$obs_phat_indiv, 1, which)

	out
}
