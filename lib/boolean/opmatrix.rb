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
    def initialize phenotype_count, orthogroup_count, storage_count
      # Make sure to provide storage_count + phenotype_count as the initial capacity, since otherwise a resize will be needed.
      super(:yale, [phenotype_count, orthogroup_count], storage_count+phenotype_count+1, :byte)
    end

    def clone_structure
      self.class.new(self.shape[0], self.shape[1], self.capacity-self.shape[0]-1)
    end

    # Basically the same as doing self[p,o] == 1
    def associated? orthogroup_id, phenotype_id
      self[phenotype_id, orthogroup_id] == 1
    end

    # Make a copy of the matrix with the row contents shuffled.
    def shuffle_rows
      t = clone_structure   # Start with the constructor for OPMatrix

      # Get the ija pointer and the diagonal contents. We'll use these to figure out the size of each row.
      STDERR.puts "shuffle_rows"

      # Create an array we can use for shuffling
      ary = (0...self.shape[1]).to_a

      (0...self.shape[0]).each do |i|
        row_size       = self.yale_d(i) + self.yale_nd_row_size(i) # diagonal plus non-diagonal elements
        new_indices    = ary.draw_without_replacement(row_size) # draw without replacement r elements from the index array

        # Get non-diagonal indices and mark the diagonal if appropriate
        if new_indices.include?(i)
          t[i,i] = 1
          new_indices.delete(i)
        end

        # Fast-insert the indices in the Yale storage
        t.yale_vector_insert(i, new_indices.to_a, [1]*new_indices.size)
      end

      return t # return the copy with the shuffled rows! Hooray!
    end

    # Gives the set of orthogroup indices for some phenotype index (index means internal/renumbered, not the ID from
    # data files).
    alias_method :orthogroups_for_phenotype, :yale_row_as_array #:yale_row_as_sorted_set
  end
end