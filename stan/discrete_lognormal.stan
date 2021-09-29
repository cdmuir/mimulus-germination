/////////////////////////////////////
// Discrete LogNormal distribution //
/////////////////////////////////////

// mu is the mean of log(y)
// sigma is the standard deviation of log(y)
// ntd is shift

// Probability mass function
  real discrete_lognormal_lpmf(int n, real mu, real sigma, real ntd) {
    return log_diff_exp(
          lognormal_lcdf(n - ntd | mu, sigma), 
          lognormal_lcdf(n - 1 - ntd | mu, sigma)
    );
  }

// Log of the complementary cumulative distribution function
  real discrete_lognormal_lccdf(int n, real mu, real sigma, real ntd) {
    return lognormal_lccdf(n - ntd | mu, sigma);
  }

// Log of the cumulative distribution function
  real discrete_lognormal_lcdf(int n, real mu, real sigma, real ntd) {
    return lognormal_lcdf(n - ntd | mu, sigma);
  }

// Random number generator
  real discrete_lognormal_rng(real mu, real sigma, real ntd) {
    return ceil(ntd + lognormal_rng(mu, sigma));
  }
