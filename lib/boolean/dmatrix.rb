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
      super([to.shape[0], from.shape[0]], dtype: :float64, stype: :dense)

      skippable_rows = to.skippable_rows
      skippable_cols = from.skippable_rows
      rows = (0...to.shape[0]).to_a - skippable_rows
      cols = (0...from.shape[0]).to_a - skippable_cols

      STDERR.puts "Filling in non-skippable entries"

      rows.each do |i|
        m_set = to.orthogroups_for_phenotype(i)

        cols.each do |j|
          # The actual calculation
          n_set = from.orthogroups_for_phenotype(j)
          k_set = m_set & n_set
          # k_set = to.yale_row_keys_intersection(i, from, j) # Can also do this instead

          #self[i,j] = 1.0 - Distribution::Hypergeometric.cdf(k_set.size-1, m_set.size, n_set.size, to.shape[1])
          begin
            self[i,j] = Hypergeometric.cdf(k_set.size, m_set.size, n_set.size, to.shape[1])
          rescue => e
            self[i,j] = Hypergeometric.ruby_cdf(k_set.size, m_set.size, n_set.size, to.shape[1])
          end
        end
      end

      STDERR.puts "Filling in skippable rows"

      # Fill in skipped entries with infinity.
      skippable_rows.each do |i|
        self[i,0...shape[1]] = Float::INFINITY # Use a slice-set to set all at once.
      end

      STDERR.puts "Filling in skippable columns"

      skippable_cols.each do |j|
        self[0...shape[0],j] = Float::INFINITY # Set the whole column at once.
      end

      STDERR.puts "Done creating DMatrix"
    end


    # Write matrix as a binary file; does not compress by default.
    def write filename, compress=nil
      STDERR.puts("Writing #{self.inspect} to #{filename}")
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
