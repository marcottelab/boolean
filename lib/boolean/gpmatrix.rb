require "nmatrix"

class GPMatrix < NMatrix
  def initialize filename, species_id
    @genes      = {}
    @phenotypes = {}
    @species_id = species_id
    @opmatrices = {} # Orthogroup-Phenotype matrices by species

    associations = {}

    gene_count = 0
    phenotype_count = 0

    f = File.new(filename, 'r')
    while line = f.gets
      line.chomp!
      phenotype, gene = line.split

      gene = gene.split(':')[1].to_i if gene.include?(':')

      # Count genes from 0 and supply each with an actual-id to count mapping
      @genes[gene] ||= begin
        x = gene_count
        gene_count += 1
        x
      end

      # Number the phenotypes.
      @phenotypes[phenotype] ||= begin
        x = phenotype_count
        phenotype_count += 1
        x
      end

      # Map the assigned gene ID to the assigned phenotype ID
      associations[@genes[gene]] = @phenotypes[phenotype]
    end

    # Call the actual matrix creation constructor with the appropriate size
    super(:yale, [phenotype_count, gene_count], associations.size, :byte)

    # Add the associations
    associations.each_pair do |j, i|
      self[i,j] = 1
    end
  end


  def opmatrix reader, species_id
    @opmatrices[species_id] ||= begin

      opm = OPMatrix.new(@phenotypes.keys.size, reader.r_to_g.keys.size,    self.capacity)

      # genes maps actual gene IDs to assigned gene IDs. Need to be able to get back the actual gene ID.
      inverted_genes = @genes.invert

      # Walk through the non-zero entries (and technically also the diagonals)
      self.each_stored_with_indices do |val,pid,gid|
        orthogroup_id = reader.orthogroup_id( inverted_genes[gid] )  # returns the renumbered orthogroup ID
        next if orthogroup_id.nil? # many genes won't exist in other species; skip these.
        opm[pid, orthogroup_id] = val
      end

      opm
    end
  end

  def associated? gene, phenotype
    self[@phenotypes[phenotype], @genes[gene]] == 1
  end

  attr_reader :genes, :phenotypes, :species_id
end
