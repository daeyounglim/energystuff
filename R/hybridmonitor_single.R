#' Run hybrid monitoring procedure for a single account
#' 
#' This is a function that run hybrid monitoring algorithm for a single account
#' 
#' @author Daeyoung Lim, \email{daeyoung.lim@uconn.edu}
#' @param y a vector of measured energy usage to be used as the response variable
#' @return a dataframe with the result of the hypothesis testing, empirical distribution of the p-values, observed p-value
#' @examples
#' \dontrun{
#' fit <- hybridmonitor.single(y, control=list(C=rep(5,12)))
#' }
#' @importFrom parallel detectCores
#' @export
'hybridmonitor.single' <- function(y, y_Kp1, control=list(), ncores=NULL, verbose=FALSE) {
	cntrl <- list(C=rep(1.5, 12), B1=10000, B2=5000, B3 = 10000, sig_level=0.05)
	cntrl[names(control)] <- control
	z_alpha <- qnorm(cntrl$sig_level, lower.tail=FALSE)
	if (!is.null(ncores)) {
		ncores_ <- parallel::detectCores()
		if (ncores > ncores_) {
			stop(paste0("The number of cores must not exceed ", ncores_))
		}
	} else {
		ncores <- parallel::detectCores()
	}
	out <- .Call(`_energystuff_hybridmonitor_single`,
		as.matrix(y),
		as.double(y_Kp1),
		as.double(cntrl$sig_level),
		as.double(z_alpha),
		as.integer(cntrl$B1),
		as.integer(cntrl$B2),
		as.integer(cntrl$B3),
		as.logical(verbose),
		as.integer(ncores)
		)
	out
}