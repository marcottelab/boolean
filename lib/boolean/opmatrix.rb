require 'nmatrix'

class OPMatrix < NMatrix
  def initialize phenotype_count, orthogroup_count, storage_count
    super(:yale, [phenotype_count, orthogroup_count], storage_count, :byte)
  end

  def associated? orthogroup_id, phenotype_id
    self[phenotype_id] == orthogroup_id
  end

  def orthogroups_for_phenotype i
    self.row(i, :copy).cast(:list, :byte).to_h[0].keys
  end
end