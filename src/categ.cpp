#include "categ.h"
#include "util.h"

// Try to make this version as fast as possible. E.g. avoid
// calling the one above.
//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector rcateg_mat(const Rcpp::NumericMatrix& P)
{
	unsigned int n = P.rows();
	unsigned int k = P.cols();

	const Rcpp::NumericVector& u = Rcpp::runif(n);
	Rcpp::IntegerVector z(n);

	for (unsigned int i = 0; i < n; i++) {
		const Rcpp::NumericVector& cp = Rcpp::cumsum(P.row(i));
		if (cp(k-1) <= 0) {
			// Watch out for the case where all probs are zero
			z(i) = NA_INTEGER;
		} else {
			z(i) = (u(i) < 1)*Rcpp::sum(u(i) > cp) + (u(i) >= 1)*(k-1);
		}
	}

	return z;
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector rcateg_logp(unsigned int n, const Rcpp::NumericVector& logprob)
{
	const Rcpp::IntegerVector& idx_sort = order(logprob, true);
	const Rcpp::NumericVector& lp = logprob[idx_sort];

	const Rcpp::NumericVector& logu = Rcpp::log(Rcpp::runif(n));
	Rcpp::IntegerVector idx_draw(n);

	unsigned int k = lp.size();
	Rcpp::NumericVector log_cp(k);
	log_cp(0) = lp(0);
	for (unsigned int j = 1; j < k; j++) {
		Rcpp::IntegerVector idx = Rcpp::seq(0, j);
		log_cp(j) = logsumprobs(lp[idx]);
	}

	for (unsigned int i = 0; i < n; i++) {
		idx_draw(i) = (logu(i) < 0)*Rcpp::sum(logu(i) > log_cp) + (logu(i) >= 0)*(k-1);
	}
	return idx_sort[idx_draw];
}

Rcpp::IntegerVector rcateg_logp_mat(const Rcpp::NumericMatrix& logprob)
{
	unsigned int n = logprob.rows();
	Rcpp::IntegerVector z(n);

	for (unsigned int i = 0; i < n; i++) {
		Rcpp::IntegerVector res = rcateg_logp(1, logprob.row(i));
		z(i) = res(0);
	}

	return z;
}

unsigned int rcateg(const Rcpp::NumericVector& p)
{
	unsigned int k = p.size();
	double u = R::runif(0, 1);
	const Rcpp::NumericVector& cp = Rcpp::cumsum(p);

	if (cp(k-1) <= 0) {
		// Watch out for the case where all probs are zero
		return NA_INTEGER;
	}

	// This is a sneaky way of finding the max index such that u(i)> cp.
	// It is possible for the RNG to draw exactly 1; we have to handle this
	// case specially.
	return (u < 1)*Rcpp::sum(u > cp) + (u >= 1)*(k-1);
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector rcateg(unsigned int n, const Rcpp::NumericVector& p)
{
	unsigned int k = p.size();
	const Rcpp::NumericVector& u = Rcpp::runif(n);
	const Rcpp::NumericVector& cp = Rcpp::cumsum(p);

	if (cp(k-1) <= 0) {
		// Watch out for the case where all probs are zero
		return Rcpp::rep(NA_INTEGER, n);
	}
	
	Rcpp::IntegerVector idx(n);
	for (unsigned int i = 0; i < n; i++) {
		// This is a sneaky way of finding the max index such that u(i)> cp.
		// It is possible for the RNG to draw exactly 1; we have to handle this
		// case specially.
		idx(i) = (u(i) < 1)*Rcpp::sum(u(i) > cp) + (u(i) >= 1)*(k-1);
	}
	return idx;
}
