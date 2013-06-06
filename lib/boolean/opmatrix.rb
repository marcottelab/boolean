require 'nmatrix'

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

    # Basically the same as doing self[p,o] == 1
    def associated? orthogroup_id, phenotype_id
      self[phenotype_id, orthogroup_id] == 1
    end

    # Gives the set of orthogroup indices for some phenotype index (index means internal/renumbered, not the ID from
    # data files).
    alias_method :orthogroups_for_phenotype, :yale_row_as_array #:yale_row_as_sorted_set
  end
end