require 'distribution'
require_relative "hypergeometric.so"

module Math
  class << self
    def ln_factorial(z)
      #Distribution::MathExtension::ApproxFactorial.stieltjes_ln_factorial(z)
      GSL::Sf::lnfact(z)     # about 100x faster.
    end
  end
end

module Hypergeometric
  class << self
    # No idea what 'cyn2' means. The function's original name was calc_rs_prob_cyn2. It calculates
    # the probability of a "Rosetta Stone fusion" (when two genes are glued together in the genome)
    # occurring by random chance. According to the comments on that function, the
    # arguments are:
    #
    # * n ::  unique hits by protein 1 in database of +total+ sequences
    # * m ::  unique hits by protein 2 in database of +total+ sequences
    # * k ::  potential rosetta stone fusion proteins linking protein 1 and protein 2
    # * total :: total number of sequences
    #
    # Returns the probability of protein 1 and protein 2 getting +k+ fusions by random chance.
    #
    # Source: Edward Marcotte, the University of Texas at Austin, 5/14/2001
    def cyn2(k, m, n, total)
      p = [n, total-n, m, total-m].inject(0.0) { |res,t| res + Math.ln_factorial(t) }
      [n-k, k, m-k, total-n-m+k, total].inject(p) { |res,t| res - Math.ln_factorial(t) }
    end

    # Slow version -- use Hypergeometric.cdf for the C version. Note that this
    # isn't technically the CDF. It's close, though. One is p <= and the other is p <.
    def ruby_cdf(k, m, n, total)
      raise(ArgumentError, "k>m") if k>m
      raise(ArgumentError, "k>n") if k>n

      min = n < m ? n : m

      sum_p = (k..min).inject(0.0) do |sum,i|
        sum + Math.exp(cyn2(i, m, n, total))
      end

      return 0.0 if sum_p < 0
      return 1.0 if sum_p > 1
      sum_p
    end

  end
end