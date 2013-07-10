require "nmatrix"

module Boolean
  # A distance matrix
  class DMatrix < NMatrix

    # Create a new distance matrix, for predicting to from from.
    #
    # Also accounts for skippable rows in real matrices. Won't matter for random matrices since
    # these don't actually have any skippable rows.
    def initialize to, from
      raise(ArgumentError, "invalid shape! [#{to.shape[0]}, #{from.shape[0]}]") if to.shape[0] == 0 || from.shape[0] == 0
      super(:dense, [to.shape[0], from.shape[0]], :float64)

      skippable_rows = to.skippable_rows
      skippable_cols = from.skippable_rows
      rows = (0...to.shape[0]).to_a - skippable_rows
      cols = (0...from.shape[0]).to_a - skippable_cols

      rows.each do |i|
        m_set = to.orthogroups_for_phenotype(i)

        cols.each do |j|
          # The actual calculation
          n_set = from.orthogroups_for_phenotype(j)
          k_set = m_set & n_set
          #self[i,j] = 1.0 - Distribution::Hypergeometric.cdf(k_set.size-1, m_set.size, n_set.size, to.shape[1])
          self[i,j] = Hypergeometric.cdf(k_set.size, m_set.size, n_set.size, to.shape[1])
        end
      end

      # Fill in skipped entries with infinity.
      skippable_rows.each do |i|
        (0...shape[1]).each do |j|
          self[i,j] = Float::INFINITY
        end
      end

      skippable_cols.each do |j|
        (0...shape[0]).each do |i|
          self[i,j] = Float::INFINITY
        end
      end
    end


    # Write matrix as a binary file; does not compress by default.
    def write filename, compress=nil
      STDERR.puts("Writing #{self.inspect}")
      x = super(filename)
      `gzip #{filename}` if compress
      x
    end

    class << self
      # Preferentially read a non-compressed file by this name, then a gzipped one.
      def read filename
        if File.exists?(filename)
          super(filename)
        else
          `gunzip -c #{filename}.gz > #{filename}`
          x = super(filename)
          `rm #{filename}`
          return x
        end
      end
    end
  end
end