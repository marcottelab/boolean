require "nmatrix"

module Boolean
  # A Gene-Phenotype matrix.
  #
  # Note that this version differs from that described in McGary et al., 2010, and Woods et al., 2013. Namely, the rows
  # and columns are swapped, as NMatrix's yale storage type is more efficient at extracting rows than columns. Since we're
  # comparing phenotypes between species, we want to be able to quickly extract those -- which means we need efficient
  # access to rows.
  class GPMatrix < NMatrix
    def initialize filename, species_id
      @genes      = {}
      @phenotypes = {}
      @species_id = species_id

      associations = Hash.new { |h,k| h[k] = [] }
      association_count = 0

      gene_count = 0
      @phenotype_count = 0

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
          x = @phenotype_count
          @phenotype_count += 1
          x
        end

        # Map the assigned gene ID to the assigned phenotype ID
        associations[@genes[gene]] << @phenotypes[phenotype]
        association_count += 1
      end

      # Call the actual matrix creation constructor with the appropriate size
      super(:yale, [@phenotype_count, gene_count], association_count + @phenotype_count, :byte)

      # Add the associations
      associations.each_pair do |j, i_array|
        i_array.each do |i|
          self[i,j] = 1
        end
      end

      require 'yaml'
      File.open( 'associations.yaml', 'w' ) do |out|
        YAML.dump(associations, out)
      end
    end


    def opmatrix reader
      opm = OPMatrix.new(@phenotype_count, reader.orthogroup_count, self.capacity)

      # genes maps actual gene IDs to assigned gene IDs. Need to be able to get back the actual gene ID.
      inverted_genes = @genes.invert

      # Walk through the non-zero entries (and technically also the diagonals)
      self.each_stored_with_indices do |val,pid,gid|
        next if val == 0
        orthogroup_id = reader.orthogroup_id( inverted_genes[gid] )  # returns the renumbered orthogroup ID
        next if orthogroup_id.nil? # many genes won't exist in other species; skip these.
        opm[pid, orthogroup_id] = val
      end

      opm.setup_skip_table!

      opm
    end

    def associated? gene, phenotype
      self[@phenotypes[phenotype], @genes[gene]] == 1
    end

    attr_reader :genes, :phenotypes, :species_id
  end
end