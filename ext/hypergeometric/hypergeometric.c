#include <ruby.h>

#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_gamma.h>

VALUE cHypergeometric;


inline double exp(const double x) {
  return gsl_sf_exp(x);
}

inline double lnfac(const unsigned int n) {
  return gsl_sf_lnfact(n);
}


static inline double cyn2(unsigned int k, unsigned int m, unsigned int n, unsigned int total) {
  return lnfac(n) + lnfac(total-n) + lnfac(m) + lnfac(total-m)
         - lnfac(n-k) - lnfac(k) - lnfac(m-k) - lnfac(total-n-m+k) - lnfac(total);
}


static VALUE hypg(VALUE self, VALUE k_, VALUE m_, VALUE n_, VALUE total_) {
  unsigned int k = FIX2INT(k_),
               m = FIX2INT(m_),
               n = FIX2INT(n_),
               total = FIX2INT(total_);

  unsigned int min = n;
  if (m < n)   min = m;

  if (k > m) rb_raise(rb_eArgError, "k > m");
  if (k > n) rb_raise(rb_eArgError, "k > n");

  double sum_p = 0.0;
  for (unsigned int i = k; i <= min; ++i) {
    sum_p += exp(cyn2(i, m, n, total));
  }

  if (sum_p < 0)      return rb_float_new(0.0);
  else if (sum_p > 1) return rb_float_new(1.0);
  else                return rb_float_new(sum_p);
}


void Init_hypergeometric() {
  cHypergeometric = rb_define_module("Hypergeometric");

  rb_define_singleton_method(cHypergeometric, "cdf", hypg, 4);
}

