require 'nmatrix'
require 'set'

class Array
  # This will only work efficiently for mathematically sparse rows.
  def draw_without_replacement(n)
    res = SortedSet.new()

    while res.size < n
      res << self[rand(self.size)]
    end

    res
  end
end


module Boolean
  # An Orthogroup-Phenotype matrix.
  #
  # Like GPMatrix, the rows and columns are swapped relative to the 2010 and 2013 papers.
  class OPMatrix < NMatrix
    # Need yale_row_as_sorted_set function
    include NMatrix::YaleFunctions

    # Create a new OPMatrix. The three arguments are the number of consecutive phenotype IDs, the number of consecutive
    # orthogroup IDs, and the initial storage capacity.
    #
    # Remember, after entering values in a non-random OPMatrix, you need to call #setup_skip_table!
    def initialize phenotype_count, orthogroup_count, storage_count
      # Make sure to provide storage_count + phenotype_count as the initial capacity, since otherwise a resize will be needed.
      super(:yale, [phenotype_count, orthogroup_count], storage_count+phenotype_count+1, :byte)
      @skip_table = {}
    end

    # Quickly tells whether a row has another just like it.
    def has_duplicate_row? row_id
      @skip_table.has_key?(row_id)
    end

    # Quickly tells whether a row has a duplicate with a lesser index.
    def is_skippable? row_id
     has_duplicate_row?(row_id) && @skip_table[row_id].is_a?(Fixnum)
    end

    # Get an array of skippable row IDs for this matrix.
    def skippable_rows
      @skip_table.values.select { |a| a.is_a?(Array) }.flatten
    end

    # Determine the number of skippable rows in the matrix. This will be used for construction of the
    # shuffled matrix.
    def count_skippables
      skippable_rows.size
    end

    # Basically the same as doing self[p,o] == 1
    def associated? orthogroup_id, phenotype_id
      self[phenotype_id, orthogroup_id] == 1
    end

    # Returns a hash. Keys are row IDs, values are either (a) arrays of other row IDs, or (b) row IDs
    # which point to where in the hash we should look for the lowest numbered duplicate row.
    #
    # Example of what gets returned:
    #   { 0 => [3,5,6],
    #     3 => 0,
    #     5 => 0,
    #     6 => 0,
    #     7 => [10,11] }
    def find_duplicate_rows_table
      h = {}
      (0...self.shape[0]).each do |i|
        i_column_indices = self.yale_row_as_array(i)
        ((i+1)...self.shape[0]).each do |j|  # j always > i
          j_column_indices = self.yale_row_as_array(j)

          next unless i_column_indices == j_column_indices

          # naive approach -- new duplicate key
          if !h.has_key?(i)
            h[i] = [j]
            h[j] = i
            next
          end

          # h already has i as a key
          if h[i].is_a?(Array)
            h[i] << j
            h[j] = i
            next
          end

          # h[i] is a single index (a Fixnum)
          h[h[i]] = j # Update the original
          h[j] = h[i] # Add a reference to the original (a Fixnum)

        end
      end
      h
    end

    def setup_skip_table!
      @skip_table = find_duplicate_rows_table
    end


    # Make a copy of the matrix with the contents of each row non-independently shuffled. The copy
    # will exclude duplicate rows from the original, so indices will no longer line up.
    def shuffle_rows
      t = clone_structure_without_skippables

      skip_offset = 0

      ary = (0...self.shape[1]).to_a.shuffle

      (0...self.shape[0]).each do |i|

        if is_skippable?(i)
          skip_offset += 1
          next
        end


        column_indices = self.yale_row_as_array(i)
        new_i = i - skip_offset

        diag = 0

        shuffled_indices = SortedSet.new()
        column_indices.each do |j|
          if ary[j] == new_i
            diag = 1
          else
            shuffled_indices << ary[j]
          end
        end

        t[new_i,new_i] = diag
        t.yale_vector_insert(new_i, shuffled_indices.to_a, [1]*shuffled_indices.size)
      end

      #t.setup_skip_table!

      return t
    end

    # Make a copy of the matrix with the contents of each row independently shuffled. Indices will
    # no longer line up due to exclusion of skipped rows
    def shuffle_each_row
      t = clone_structure_without_skippables   # Start with the constructor for OPMatrix

      skip_offset = 0

      # Create an array we can use for shuffling
      ary = (0...self.shape[1]).to_a

      (0...self.shape[0]).each do |i|

        if is_skippable?(i)
          skip_offset += 1
          next
        end

        row_size       = self.yale_d(i) + self.yale_nd_row_size(i) # diagonal plus non-diagonal elements
        new_indices    = ary.draw_without_replacement(row_size) # draw without replacement r elements from the index array

        new_i = i - skip_offset

        # Get non-diagonal indices and mark the diagonal if appropriate
        if new_indices.include?(new_i)
          t[new_i,new_i] = 1
          new_indices.delete(new_i)
        end

        # Fast-insert the indices in the Yale storage
        t.yale_vector_insert(new_i, new_indices.to_a, [1]*new_indices.size)
      end

      return t # return the copy with the shuffled rows! Hooray!
    end

    # Gives the set of orthogroup indices for some phenotype index (index means internal/renumbered, not the ID from
    # data files).
    alias_method :orthogroups_for_phenotype, :yale_row_as_array #:yale_row_as_sorted_set

    class << self
      # Make sure we re-calculate duplicate rows upon matrix loading
      def read *args
        x = super(*args)
        binding.pry
        x.send :eval, "@skip_table = {}"
        x
      end
    end
  protected

    def clone_structure
      self.class.new(self.shape[0], self.shape[1], self.capacity-self.shape[0]-1)
    end

    def clone_structure_without_skippables
      new_shape0 = self.shape[0] - count_skippables
      # This is going to allocate more space than is strictly necessary, since we're removing many
      # duplicates.
      self.class.new(new_shape0, self.shape[1], self.capacity-new_shape0-1)
    end
  end
end