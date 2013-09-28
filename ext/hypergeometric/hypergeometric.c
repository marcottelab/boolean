#include <ruby.h>

#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_gamma.h>

VALUE cHypergeometric;

/*
 * Warning: e**-708 = 3.308E-308. e**-709 = underflow error!
 */
inline double exp(const double x) {
  return gsl_sf_exp(x);
}

static VALUE rb_exp(VALUE self, VALUE x) {
  return rb_float_new(exp(NUM2DBL(x)));
}

inline double lnfac(const unsigned int n) {
  return gsl_sf_lnfact(n);
}

static VALUE rb_lnfac(VALUE self, VALUE n) {
  return rb_float_new(lnfac(FIX2INT(n)));
}


static inline double cyn2(unsigned int k, unsigned int m, unsigned int n, unsigned int total) {
  return lnfac(n) + lnfac(total-n) + lnfac(m) + lnfac(total-m)
         - lnfac(n-k) - lnfac(k) - lnfac(m-k) - lnfac(total-n-m+k) - lnfac(total);
}


static VALUE rb_cyn2(VALUE self, VALUE k, VALUE m, VALUE n, VALUE total) {
  return rb_float_new(cyn2(FIX2INT(k), FIX2INT(m), FIX2INT(n), FIX2INT(total)));
}


/* See the note above on underflow error! If that happens, we'll just skip the rest of
 * the hypergeometric calculation and print a warning.
 *
 * However, let this be a warning to ye, all who enter here. This lib does a better job
 * than GSL of calculating hypergeometric CDF, but
 */
static VALUE rb_hypg(VALUE self, VALUE k_, VALUE m_, VALUE n_, VALUE total_) {
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
    //double c = cyn2(i, m, n, total);

    // Check for underflow before calculating the exponent. We'll let the user know if it's irrecoverable.
    //if (c < -708.3964) {
    //  if (i == k) rb_raise(rb_eNotImpError, "GSL underflowed on p-value calculation, can't recover");
    //  rb_warn("GSL underflow with %u-th component of p-value calculation, returning %f (you may wish to use ruby_cdf instead)", i-k, sum_p);
    //  break;
    //}
    // More efficient to just handle this by calling ruby_cdf regardless.

    sum_p += exp(cyn2(i, m, n, total));
  }

  if (sum_p < 0)      return rb_float_new(0.0);
  else if (sum_p > 1) return rb_float_new(1.0);
  else                return rb_float_new(sum_p);
}


void Init_hypergeometric() {
  cHypergeometric = rb_define_module("Hypergeometric");

  rb_define_singleton_method(cHypergeometric, "cdf", rb_hypg, 4);
  rb_define_singleton_method(cHypergeometric, "cyn2", rb_cyn2, 4);
  rb_define_singleton_method(cHypergeometric, "exp", rb_exp, 1);
  rb_define_singleton_method(cHypergeometric, "lnfac", rb_lnfac, 1);
}

