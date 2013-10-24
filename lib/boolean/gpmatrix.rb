require "nmatrix"

module Boolean
  # A Gene-Phenotype matrix.
  #
  # Note that this version differs from that described in McGary et al., 2010, and Woods et al., 2013. Namely, the rows
  # and columns are swapped, as NMatrix's yale storage type is more efficient at extracting rows than columns. Since we're
  # comparing phenotypes between species, we want to be able to quickly extract those -- which means we need efficient
  # access to rows.
  class GPMatrix < NMatrix
    include NMatrix::YaleFunctions

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
        #associations[@genes[gene]] << @phenotypes[phenotype]
        associations[@phenotypes[phenotype]] << @genes[gene]
        association_count += 1
      end

      # Call the actual matrix creation constructor with the appropriate size
      super([@phenotype_count, gene_count], association_count + @phenotype_count, dtype: :byte, stype: :yale)

      # Add the associations
      #associations.each_pair do |i, j_array|
      #  j_array.each do |j|
      #    self[i,j] = 1
      #  end
      #end

      associations.each_pair do |i, j_array|
        self[i,i] = 1 if j_array.delete(i) # Set diagonal separately
        self.__yale_vector_set__(i, j_array, [1]*j_array.size)
      end
    end

    alias_method :genes_for_phenotype, :yale_ja_d_keys_at


    # Do not call this function until construction is finished.
    def phenotype_ids
      @phenotype_ids ||= @phenotypes.invert
    end

    def gene_ids
      @gene_ids ||= @genes.invert
    end


    def opmatrix reader
      opm = OPMatrix.new(@phenotype_count, reader.orthogroup_count, self.capacity)

      # genes maps actual gene IDs to assigned gene IDs. Need to be able to get back the actual gene ID.
      inverted_genes = gene_ids

      # Walk through the non-zero entries (and technically also the diagonals).
      # Note that this could be done more efficiently with GPMatrix#genes_for_phenotype and #NMatrix#__yale_vector__set.
      # But it's going quickly enough relative to other operations that we'll just leave it here for now.
      self.each_stored_with_indices do |val,pid,gid|
        next if val == self.default_value
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