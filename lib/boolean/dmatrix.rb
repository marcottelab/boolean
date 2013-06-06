require "nmatrix"

module Boolean
  # A distance matrix
  class DMatrix < NMatrix

    # Create a new distance matrix, for predicting to from from
    def initialize to, from
      super(:dense, [to.shape[0], from.shape[0]], :float64)

      (0...to.shape[0]).each do |i|
        m_set = to.orthogroups_for_phenotype(i)

        (0...from.shape[0]).each do |j|
          n_set = from.orthogroups_for_phenotype(j)
          k_set = m_set & n_set
          self[i,j] = 1.0 - Distribution::Hypergeometric.cdf(k_set.size, m_set.size, n_set.size, to.shape[1])
        end
      end
    end
  end
end