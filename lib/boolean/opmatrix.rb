require 'nmatrix'

# An Orthogroup-Phenotype matrix.
#
# Like GPMatrix, the rows and columns are swapped relative to the 2010 and 2013 papers.
class OPMatrix < NMatrix
  # Create a new OPMatrix. The three arguments are the number of consecutive phenotype IDs, the number of consecutive
  # orthogroup IDs, and the initial storage capacity.
  def initialize phenotype_count, orthogroup_count, storage_count
    super(:yale, [phenotype_count, orthogroup_count], storage_count, :byte)
  end

  # Basically the same as doing self[p,o] == 1
  def associated? orthogroup_id, phenotype_id
    self[phenotype_id, orthogroup_id] == 1
  end

  # Gives the set of orthogroup indices for some phenotype index (index means internal/renumbered, not the ID from
  # data files).
  def orthogroups_for_phenotype i
    self.row(i, :copy).to_h[0].keys
  end
end