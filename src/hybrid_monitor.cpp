#include <iostream>
#include <cmath>
#include <RcppArmadillo.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <algorithm>
#include <iterator>
#include <progress.hpp>
#include <progress_bar.hpp>
#include "quantile.h"
#ifdef _OPENMP
  #include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo,RcppProgress)]]

// [[Rcpp::export]]
Rcpp::List hybridmonitor_single(const arma::mat& y,
						        const arma::vec& y_Kp1, // = y_{K+1})
						        const double& sig_level,
					 	        const double& z_alpha,
					 	        const int& B1,
					 	        const int& B2,
					 	        const int& B3,
					 	        const bool verbose,
					 	        const int& ncores) {
	using namespace arma;
	const int T = y.n_cols; // # of years

	/********************************************
	Get beta0 and sigma2_hat with historical data
	********************************************/
	vec betahat = arma::mean(y, 1);
	double sig2hat = 0.0;
	for (int j = 0; j < T; ++j) {
		sig2hat += arma::accu(arma::square(y.col(j) - betahat));
	}
	sig2hat /= 12.0 * static_cast<double>(T-1);

	/****************************************************************
	Construct an empirical distribution of acceptance probability and
	compute the critical value with simulated future data
	****************************************************************/
	vec P(B2, fill::zeros);
	if (verbose) {
		Rcpp::Rcout << "Constructing empirical distribution of acceptance probabilities" << std::endl;
	}
	{	
		Progress prog(B2, verbose);
		#ifdef _OPENMP
		#pragma omp parallel for schedule(static) num_threads(ncores)
		#endif
		for (int b2 = 0; b2 < B2; ++b2) {
			vec delta(B1, fill::zeros);
			vec betaf(12, fill::randn);
			double sig2f = 6.0 * static_cast<double>(T-1) * sig2hat / ::Rf_rgamma(6.0 * static_cast<double>(T-1), 1.0);
			betaf *= std::sqrt(sig2f / static_cast<double>(T));
			betaf += betahat;
			vec yf(12, fill::randn);
			yf *= std::sqrt(sig2f);
			yf += betaf;
			for (int b1 = 0; b1 < B1; ++b1) {
				vec beta(12, fill::randn);
				double sig2 = 6.0 * static_cast<double>(T-1) * sig2hat / ::Rf_rgamma(6.0 * static_cast<double>(T-1), 1.0);
				beta *= std::sqrt(sig2 / static_cast<double>(T));
				beta += betahat;
				delta(b1) = (arma::dot(beta, yf - beta) / (std::sqrt(sig2) * arma::norm(beta)) < z_alpha) ? 1.0 : 0.0;
			}
			P(b2) = arma::mean(delta);
			prog.increment();
		}
	}
	double gamma_alpha = quantile_cpp(P, sig_level);

	/*******************************************************
	Estimate the acceptance rate with observed (K+1)-th data
	*******************************************************/
	vec delta(B3, fill::zeros);
	if (verbose) {
		Rcpp::Rcout << "Computing observed acceptance rate" << std::endl;
	}
	{
		Progress prog(B3, verbose);
		#ifdef _OPENMP
		#pragma omp parallel for schedule(static) num_threads(ncores)
		#endif
		for (int b3 = 0; b3 < B3; ++b3) {
			vec beta(12, fill::randn);
			double sig2 = 6.0 * static_cast<double>(T-1) * sig2hat / ::Rf_rgamma(6.0 * static_cast<double>(T-1), 1.0);
			beta *= std::sqrt(sig2 / static_cast<double>(T));
			beta += betahat;
			delta(b3) = (arma::dot(beta, y_Kp1 - beta) / (std::sqrt(sig2) * arma::norm(beta)) < z_alpha) ? 1.0 : 0.0;	
		}
		prog.increment();
	}
	double phat = arma::mean(delta);
	return Rcpp::List::create(Rcpp::Named("test") = (phat >= gamma_alpha) ? 0 : 1, // test function (= decision whether to reject H0)
							  Rcpp::Named("emp_dist") = P, // empirical distribution of p-values
							  Rcpp::Named("cr_val") = gamma_alpha, // critical value from empirical distribution
							  Rcpp::Named("obs_phat") = phat, // observed p-value
							  Rcpp::Named("betahat") = betahat,
							  Rcpp::Named("sig2hat") = sig2hat); 

}

// [[Rcpp::export]]
Rcpp::List hybridmonitor_multiple(const arma::cube& y, // 12 x K x T (one slice per account)
						          const arma::mat& y_Kp1, // = y_{K+1}) (12 x K)
						          const double& sig_level,
					 	          const double& z_alpha,
					 	          const int& B1,
					 	          const int& B2,
					 	          const int& B3,
					 	          const int& T,
					 	          const bool verbose,
					 	          const int& ncores) {
	using namespace arma;
	const int K = y_Kp1.n_cols;

	/********************************************************
	Get widehat{beta} and widehat{Sigma} with historical data
	********************************************************/
	mat betahat = arma::mean(y, 2);
	mat S(12, 12, fill::zeros);
	for (int t = 0; t < T; ++t) {
		mat y_k = y.slice(t);
		mat resid_k = y_k - betahat;
		S += resid_k * resid_k.t();
	}
	mat Sinv = arma::inv_sympd(S);
	mat Schol = arma::chol(Sinv); // To sample from iwishrnd; faster with Cholesky factor
	const double df = static_cast<double>(K * (T - 1) - 13);
	/****************************************************************
	Construct an empirical distribution of acceptance probability and
	compute the critical value with simulated future data
	****************************************************************/
	vec P(B2, fill::zeros);
	if (verbose) {
		Rcpp::Rcout << "Constructing empirical distribution of acceptance probabilities" << std::endl;
	}
	{	
		Progress prog(B2, verbose);
		#ifdef _OPENMP
		#pragma omp parallel for schedule(static) num_threads(ncores)
		#endif
		for (int b2 = 0; b2 < B2; ++b2) {
			// Generate Sig_f ~ Inverse-Wishart(K*(T-1)-13, S)
			mat Sigf = arma::iwishrnd(S, df, Schol);
			mat Sigfchol = arma::chol(Sigf);
			// Generate beta_f ~ Normal(betahat, Sig_f/T)
			mat betaf(12, K, fill::randn);
			betaf = (Sigfchol.t() * betaf) / std::sqrt(static_cast<double>(T));
			betaf += betahat; 
			// Generate y_f ~ Normal(beta_f, Sig_f)
			mat yf(12, K, fill::randn);
			yf = Sigfchol.t() * yf + betaf;
			vec delta(B1, fill::zeros);
			for (int b1 = 0; b1 < B1; ++b1) {
				// Generate Sig ~ Inverse-Wishart(K*(T-1)-13, S)
				mat Sig = arma::iwishrnd(S, df, Schol);
				mat Sigchol = arma::chol(Sig);
				// Generate beta ~ Normal(betahat, Sig/T)
				mat beta(12, K, fill::randn);
				beta = (Sigchol.t() * beta) / std::sqrt(static_cast<double>(T)) + betahat;
				// Compute 
				rowvec num = arma::sum(beta % arma::solve(arma::trimatu(Sigchol), arma::solve(arma::trimatl(Sigchol.t()), yf - beta)));
				rowvec denom = arma::sqrt(arma::sum(beta % arma::solve(arma::trimatu(Sigchol),arma::solve(arma::trimatl(Sigchol.t()), beta))));
				rowvec Ws = num / denom;
				delta(b1) = Ws.max() < z_alpha ? 1.0 : 0.0;
			}
			P(b2) = arma::mean(delta);
			prog.increment();
		}
	}
	double gamma_alpha = quantile_cpp(P, sig_level);

	/*******************************************************
	Estimate the acceptance rate with observed (K+1)-th data
	*******************************************************/
	vec delta(B3, fill::zeros);
	uvec idx(B3, fill::zeros);
	vec Wss(B3, fill::zeros);
	if (verbose) {
		Rcpp::Rcout << "Computing observed acceptance rate" << std::endl;
	}
	{
		Progress prog(B3, verbose);
		#ifdef _OPENMP
		#pragma omp parallel for schedule(static) num_threads(ncores)
		#endif
		for (int b3 = 0; b3 < B3; ++b3) {
			mat Sig = arma::iwishrnd(S, df, Schol);
			mat Sigchol = arma::chol(Sig);
			mat beta(12, K, fill::randn);
			beta = (Sigchol.t() * beta) / std::sqrt(static_cast<double>(T)) + betahat;
			rowvec num = arma::sum(beta % arma::solve(arma::trimatu(Sigchol), arma::solve(arma::trimatl(Sigchol.t()), y_Kp1 - beta)));
			rowvec denom = arma::sqrt(arma::sum(beta % arma::solve(arma::trimatu(Sigchol), arma::solve(arma::trimatl(Sigchol.t()), beta))));
			rowvec Ws = num / denom;
			uword i = Ws.index_max();
			idx(b3) = i;
			delta(b3) = Ws(i) < z_alpha ? 1.0 : 0.0;
			Wss(b3) = Ws(i);
		}
		prog.increment();
	}
	double phat = arma::mean(delta);
	return Rcpp::List::create(Rcpp::Named("test") = (phat >= gamma_alpha) ? 0 : 1, // test function (= decision whether to reject H0)
							  Rcpp::Named("idx") = idx,
							  Rcpp::Named("Wmax") = Wss,
							  Rcpp::Named("emp_dist") = P, // empirical distribution of acceptance probabilities
							  Rcpp::Named("cr_val") = gamma_alpha, // critical value from empirical distribution
							  Rcpp::Named("obs_phat") = phat, // observed acceptance probability
							  Rcpp::Named("betahat") = betahat);
}

template <typename T>
T cummax(const T& x) {
	using namespace arma;
	int n = x.n_elem;
	T out(n, fill::zeros);
	out(0) = x(0);
	for (int i = 1; i < n; ++i) {
		out(i) = (x(i) >= out(i-1)) ? x(i) : out(i-1);
	}
	return out;
}

template <typename T>
T cummin(const T& x) {
	using namespace arma;
	int n = x.n_elem;
	T out(n, fill::zeros);
	out(0) = x(0);
	for (int i = 1; i < n; ++i) {
		out(i) = (x(i) <= out(i-1)) ? x(i) : out(i-1);
	}
	return out;
}

template <typename T>
T bonferroni(const T& P) {
	using namespace arma;
	uvec o = sort_index(P);
	uvec ro = sort_index(o);
	int n = P.n_elem;
	T ones(n, fill::ones);
	T cm = static_cast<double>(n) * conv_to<T>::from(P(o));
	T out = arma::min(ones, cm);
	return conv_to<T>::from(out(ro));
}

template <typename T>
T holm(const T& P) {
	using namespace arma;
	int n = P.n_elem;
	uvec o = sort_index(P);
	uvec ro = sort_index(o);
	T p_ordered = arma::conv_to<T>::from(P(o));
	T out(n, fill::zeros);
	out(0) = 1.0 - std::pow(1.0 - p_ordered(0), static_cast<double>(n));
	for (int k = 1; k < n; ++k) {
		double tmp = 1.0 - std::pow(1.0 - p_ordered(k), static_cast<double>(n));
		out(k) = (out(k-1) >= tmp) ? out(k-1) : tmp;
	}
	return conv_to<T>::from(out(ro));
}

template <typename T>
T hochberg(const T& P) {
	using namespace arma;
	int n = P.n_elem;
	T i = arma::linspace<T>(static_cast<double>(n), 1.0, n);
	uvec o = sort_index(P, "descend");
	uvec ro = sort_index(o);
	T adj_p = (static_cast<double>(n) - i + 1.0) % conv_to<T>::from(P(o));
	T ones(n, fill::ones);
	T cm = cummin(adj_p);
	T out = arma::min(ones, cm);
	return conv_to<T>::from(out(ro));
}

template <typename T>
T hommel(const T& P_) {
	using namespace arma;
	int n = P_.n_elem;
	T i = arma::linspace<T>(1.0, static_cast<double>(n), n);
	uvec o = sort_index(P_);
	T P = conv_to<T>::from(P_(o));
	uvec ro = sort_index(o);
	double min_val = arma::min(static_cast<double>(n) * P / i);
	T q(n, fill::zeros);
	q.fill(min_val);
	T pa(n, fill::zeros);
	pa.fill(min_val);
	for (int j = n-1; j > 1; --j) {
		int ij = n - j + 1;
		int i2 = j - 1;
		T denom = linspace<T>(2.0, static_cast<double>(j), i2);
		T p_tail = P.tail(i2) * static_cast<double>(j) / denom;
		double q1 = p_tail.min();
		T p_head = P.head(ij) * static_cast<double>(j);
		p_head.elem( find(p_head > q1) ).fill(q1);
		q.head(ij) = p_head;
		double qtmp = q(n-j);
		q.tail(i2).fill(qtmp);
		pa = arma::max(pa, q);
	}
	T out = arma::max(pa, P);
	return conv_to<T>::from(out(ro));
}

template <typename T>
T BH(const T& P) {
	using namespace arma;
	int n = P.n_elem;
	T i = arma::linspace<T>(static_cast<double>(n), 1.0, n);
	uvec o = sort_index(P, "descend");
	uvec ro = sort_index(o);

	T adj_p = (static_cast<double>(n) / i) % conv_to<T>::from(P(o));
	T ones(n, fill::ones);
	T cm = cummin(adj_p);
	T out = arma::min(ones, cm);
	return conv_to<T>::from(out(ro));
}

template <typename T>
T BY(const T& P) {
	using namespace arma;
	int n = P.n_elem;
	T i = arma::linspace<T>(static_cast<double>(n), 1.0, n);
	uvec o = sort_index(P, "descend");
	uvec ro = sort_index(o);

	double q = arma::accu(1.0 / i);
	T adj_p = (q * static_cast<double>(n) / i) % conv_to<T>::from(P(o));
	T ones(n, fill::ones);
	T cm = cummin(adj_p);
	T out = arma::min(ones, cm);
	return conv_to<T>::from(out(ro));
}

// [[Rcpp::export]]
Rcpp::List hybridmonitor_multiple_select_months(const arma::cube& y, // 12 x K x T (one slice per account)
										           const arma::mat& y_Kp1_, // = y_{K+1}) (12 x K)
										           const double& sig_level,
										           const double& z_alpha_K,
										           const int& B1,
										           const int& B2,
										           const int& B3,
										           const int& T,
										           const arma::uvec& m_idx, // the indices of months to test
										           const bool do_bonferroni,
										           const bool do_holm,
										           const bool do_hochberg,
										           const bool do_hommel,
										           const bool do_BH,
										           const bool do_BY,
										           const bool verbose,
										           const int& ncores) {
	using namespace arma;
	const int K = y_Kp1_.n_cols;
	const int J = m_idx.n_elem;
	const double upper_sig_level = 1.0 - sig_level;
	const int n_multiple_testing = static_cast<int>(do_bonferroni) + static_cast<int>(do_holm) + static_cast<int>(do_hochberg) +
									static_cast<int>(do_hommel) + static_cast<int>(do_BH) + static_cast<int>(do_BY);
	mat y_Kp1 = y_Kp1_.rows(m_idx);
	/********************************************************
	Get \widehat{beta} and \widehat{Sigma} with historical data
	********************************************************/
	mat betahat_ = arma::mean(y, 2);
	mat betahat = betahat_.rows(m_idx);
	mat S_(12, 12, fill::zeros);
	for (int t = 0; t < T; ++t) {
		mat y_k = y.slice(t);
		mat resid_k = y_k - betahat_;
		S_ += resid_k * resid_k.t();
	}
	mat S = S_(m_idx, m_idx);
	mat Sinv = arma::inv_sympd(S);

	mat Schol = arma::chol(Sinv); // To sample from iwishrnd; faster with Cholesky factor
	const double df = static_cast<double>(K * (T - 1) - J - 1);
	/****************************************************************
	Construct an empirical distribution of p-value and
	compute the critical value with simulated future data
	****************************************************************/
	vec P(B2, fill::zeros);
	cube P_multiple_testing(n_multiple_testing, K, B2, fill::zeros);
	mat P_(n_multiple_testing, B2, fill::zeros);
	if (verbose) {
		Rcpp::Rcout << "Constructing empirical distribution of p-values" << std::endl;
	}
	{	
		Progress prog(B2, verbose);
		#ifdef _OPENMP
		#pragma omp parallel for schedule(static) num_threads(ncores)
		#endif
		for (int b2 = 0; b2 < B2; ++b2) {
			if (!Progress::check_abort()) {
				// Generate Sig_f ~ Inverse-Wishart(K*(T-1)-13, S)
				mat Sigf = arma::iwishrnd(S, df, Schol);
				mat Sigfchol = arma::chol(Sigf);
				// Generate beta_f ~ Normal(betahat, Sig_f/T)
				mat betaf(J, K, fill::randn);
				betaf = (Sigfchol.t() * betaf) / std::sqrt(static_cast<double>(T));
				betaf += betahat; 
				// Generate y_f ~ Normal(beta_f, Sig_f)
				mat yf(J, K, fill::randn);
				yf = Sigfchol.t() * yf + betaf;
				double delta = 0.0;
				mat delta_multiple_testing(n_multiple_testing, K, fill::zeros);
				vec delta_(n_multiple_testing, fill::zeros);
				for (int b1 = 0; b1 < B1; ++b1) {
					// Generate Sig ~ Inverse-Wishart(K*(T-1)-13, S)
					mat Sig = arma::iwishrnd(S, df, Schol);
					mat Sigchol = arma::chol(Sig);
					// Generate beta ~ Normal(betahat, Sig/T)
					mat beta(J, K, fill::randn);
					beta = (Sigchol.t() * beta) / std::sqrt(static_cast<double>(T)) + betahat;
					// Compute test statistic and adjusted p-values
					rowvec num = arma::sum(beta % arma::solve(arma::trimatu(Sigchol), arma::solve(arma::trimatl(Sigchol.t()), yf - beta)));
					rowvec denom = arma::sqrt(arma::sum(beta % arma::solve(arma::trimatu(Sigchol),arma::solve(arma::trimatl(Sigchol.t()), beta))));
					rowvec Ws = num / denom;
					delta += (Ws.max() < z_alpha_K) ? 1.0 : 0.0;
					rowvec Us = 1.0 - arma::normcdf(Ws);
					int counter = 0;
					if (do_bonferroni) {
						rowvec p_adj = bonferroni(Us);
						urowvec test_multiple_testing = p_adj <= sig_level;
						delta_(counter) += arma::any(test_multiple_testing) ? 1.0 : 0.0;
						delta_multiple_testing.row(counter) += arma::conv_to<rowvec>::from(test_multiple_testing);
						++counter;
					}
					if (do_holm) {
						rowvec p_adj = holm(Us);
						urowvec test_multiple_testing = p_adj <= sig_level;
						delta_(counter) += arma::any(test_multiple_testing) ? 1.0 : 0.0;
						delta_multiple_testing.row(counter) += arma::conv_to<rowvec>::from(test_multiple_testing);
						++counter;
					}
					if (do_hochberg) {
						rowvec p_adj = hochberg(Us);
						urowvec test_multiple_testing = p_adj <= sig_level;
						delta_(counter) += arma::any(test_multiple_testing) ? 1.0 : 0.0;
						delta_multiple_testing.row(counter) += arma::conv_to<rowvec>::from(test_multiple_testing);
						++counter;
					}
					if (do_hommel) {
						rowvec p_adj = hommel(Us);
						urowvec test_multiple_testing = p_adj <= sig_level;
						delta_(counter) += arma::any(test_multiple_testing) ? 1.0 : 0.0;
						delta_multiple_testing.row(counter) += arma::conv_to<rowvec>::from(test_multiple_testing);
						++counter;
					}
					if (do_BH) {
						rowvec p_adj = BH(Us);
						urowvec test_multiple_testing = p_adj <= sig_level;
						delta_(counter) += arma::any(test_multiple_testing) ? 1.0 : 0.0;
						delta_multiple_testing.row(counter) += arma::conv_to<rowvec>::from(test_multiple_testing);
						++counter;
					}
					if (do_BY) {
						rowvec p_adj = BY(Us);
						urowvec test_multiple_testing = p_adj <= sig_level;
						delta_(counter) += arma::any(test_multiple_testing) ? 1.0 : 0.0;
						delta_multiple_testing.row(counter) += arma::conv_to<rowvec>::from(test_multiple_testing);
						++counter;
					}
				}
				P(b2) = delta / static_cast<double>(B1);
				P_.col(b2) = delta_ / static_cast<double>(B1);
				P_multiple_testing.slice(b2) = delta_multiple_testing / static_cast<double>(B1);
				prog.increment();
			}
		}
	}
	double gamma_alpha = quantile_cpp(P, upper_sig_level);
	vec gamma_alpha_ = quantile_v_cpp(P_, upper_sig_level);
	mat gamma_alpha_multiple_testing = quantile_m_cpp(P_multiple_testing, upper_sig_level);
	

	/*******************************************************
	Estimate the p-value with observed (K+1)-th data
	*******************************************************/
	double delta = 0.0;
	vec delta_(n_multiple_testing, fill::zeros);
	mat delta_multiple_testing(n_multiple_testing, K, fill::zeros);
	uvec idx(B3, fill::zeros);
	vec Wss(B3, fill::zeros);
	if (verbose) {
		Rcpp::Rcout << "Computing observed p-value" << std::endl;
	}
	{
		Progress prog(B3, verbose);
		#ifdef _OPENMP
		#pragma omp parallel for schedule(static) num_threads(ncores)
		#endif
		for (int b3 = 0; b3 < B3; ++b3) {
			if (!Progress::check_abort()) {
				mat Sig = arma::iwishrnd(S, df, Schol);
				mat Sigchol = arma::chol(Sig);
				mat beta(J, K, fill::randn);
				beta = (Sigchol.t() * beta) / std::sqrt(static_cast<double>(T)) + betahat;
				rowvec num = arma::sum(beta % arma::solve(arma::trimatu(Sigchol), arma::solve(arma::trimatl(Sigchol.t()), y_Kp1 - beta)));
				rowvec denom = arma::sqrt(arma::sum(beta % arma::solve(arma::trimatu(Sigchol), arma::solve(arma::trimatl(Sigchol.t()), beta))));
				rowvec Ws = num / denom;
				uword i = Ws.index_max();
				idx(b3) = i;
				delta += (Ws(i) > z_alpha_K) ? 1.0 : 0.0;

				rowvec Us = 1.0 - arma::normcdf(Ws);

				int counter = 0;
				if (do_bonferroni) {
					rowvec p_adj = bonferroni(Us);
					urowvec test_multiple_testing = p_adj <= sig_level;
					delta_(counter) += arma::any(test_multiple_testing) ? 1.0 : 0.0;
					delta_multiple_testing.row(counter) += arma::conv_to<rowvec>::from(test_multiple_testing);
					++counter;
				}
				if (do_holm) {
					rowvec p_adj = holm(Us);
					urowvec test_multiple_testing = p_adj <= sig_level;
					delta_(counter) += arma::any(test_multiple_testing) ? 1.0 : 0.0;
					delta_multiple_testing.row(counter) += arma::conv_to<rowvec>::from(test_multiple_testing);
					++counter;
				}
				if (do_hochberg) {
					rowvec p_adj = hochberg(Us);
					urowvec test_multiple_testing = p_adj <= sig_level;
					delta_(counter) += arma::any(test_multiple_testing) ? 1.0 : 0.0;
					delta_multiple_testing.row(counter) += arma::conv_to<rowvec>::from(test_multiple_testing);
					++counter;
				}
				if (do_hommel) {
					rowvec p_adj = hommel(Us);
					urowvec test_multiple_testing = p_adj <= sig_level;
					delta_(counter) += arma::any(test_multiple_testing) ? 1.0 : 0.0;
					delta_multiple_testing.row(counter) += arma::conv_to<rowvec>::from(test_multiple_testing);
					++counter;
				}
				if (do_BH) {
					rowvec p_adj = BH(Us);
					urowvec test_multiple_testing = p_adj <= sig_level;
					delta_(counter) += arma::any(test_multiple_testing) ? 1.0 : 0.0;
					delta_multiple_testing.row(counter) += arma::conv_to<rowvec>::from(test_multiple_testing);
					++counter;
				}
				if (do_BY) {
					rowvec p_adj = BY(Us);
					urowvec test_multiple_testing = p_adj <= sig_level;
					delta_(counter) += arma::any(test_multiple_testing) ? 1.0 : 0.0;
					delta_multiple_testing.row(counter) += arma::conv_to<rowvec>::from(test_multiple_testing);
					++counter;
				}
				Wss(b3) = Ws(i);
			}
		}
		prog.increment();
	}
	double phat = delta / static_cast<double>(B3);
	vec phat_ = delta_ / static_cast<double>(B3);
	mat phat_multiple_testing = delta_multiple_testing / static_cast<double>(B3);
	
	return Rcpp::List::create(Rcpp::Named("test") = (phat > gamma_alpha) ? 1 : 0, // test function (= decision whether to reject H0)
							  Rcpp::Named("test_grand") = phat_ > gamma_alpha_,
							  Rcpp::Named("idx") = idx,
							  Rcpp::Named("Wmax") = Wss,
							  Rcpp::Named("emp_dist") = P, // empirical distribution of p-values
							  Rcpp::Named("emp_dist_grand") = P_, // empirical distribution of p-values
							  Rcpp::Named("emp_dist_indiv") = P_multiple_testing, // empirical distribution of p-values
							  Rcpp::Named("cr_val") = gamma_alpha, // critical value from empirical distribution
							  Rcpp::Named("cr_val_grand") = gamma_alpha_, // critical value from empirical distribution
							  Rcpp::Named("cr_val_indiv") = gamma_alpha_multiple_testing, // critical value from empirical distribution
							  Rcpp::Named("obs_phat") = phat, // observed p-value
							  Rcpp::Named("obs_phat_grand") = phat_, // observed p-value
							  Rcpp::Named("obs_phat_indiv") = phat_multiple_testing);
}