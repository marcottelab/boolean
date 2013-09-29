module Boolean
  # Orthogroup Reader
  #
  # This class parses the sqltable file (INPARANOID output) for a given species pair. It is needed for converting a
  # gene-phenotype matrix (GPMatrix) to an orthogroup-phenotype matrix (OPMatrix).
  class OrthoReader

    def inspect #:nodoc:
      ["#<Boolean::OrthoReader:0x00#{(self.object_id << 1).to_s(16)}",
       "g_to_o.size=#{@g_to_o.size}",
       "r_to_g.size=#{@r_to_g.size}",
       "@renumber.size=#{@renumber.size}",
      "@orthogroup_count=#{@orthogroup_count}"].join(' ') + ">"
    end


    attr_reader :g_to_o, :r_to_g, :renumber, :orthogroup_count, :genes_by_species

    # Convert a (renumbered) gene ID (a GPMatrix index) to our renumbered orthogroup ID (an OPMatrix index)
    def orthogroup_id gid
      @renumber[@g_to_o[gid]]
    end

    # Convert a renumbered orthogroup ID (OPMatrix index) to a renumbered gene ID (GPMatrix index).
    def gene_id rid
      @r_to_g[rid]
    end

    def o_to_g
      @o_to_g ||= begin
        h = Hash.new { |h,k| h[k] = [] }
        @g_to_o.each_pair do |gene,orthogroup|
          h[orthogroup] << gene
        end
        h
      end
    end

    # Map internal orthogroup IDs to INPARANOID orthogroup IDs.
    def unrenumber
      @unrenumber ||= begin
        h = Hash.new { |h,k| h[k] = [] }
        @renumber.each_pair do |inparanoid_id, oid|
          h[oid] << inparanoid_id
        end
        h
      end
    end


    # Create an OrthoReader for a pair of species (given by +species_ids+). Optionally, can supply the data directory as
    # second argument, which defaults to "data"
    def initialize species_ids, dir="data"
      @renumber = {} # map INPARANOID orthogroup IDs to renumbered oids

      # Keep track of genes for each species in case we want to do predictions.
      @genes_by_species = {
          species_ids[0] => SortedSet.new,
          species_ids[1] => SortedSet.new
      }

      Dir.chdir(dir) do
        possible_filenames = ["sqltable.#{species_ids[0]}-#{species_ids[1]}", "sqltable.#{species_ids[1]}-#{species_ids[0]}"]
        filename = File.exists?(possible_filenames.first) ? possible_filenames.first : possible_filenames.last
        File.open(filename, "r") do |f|

          # genes to orthogroups and orthogroups to genes
          @g_to_o = {}
          @r_to_g = Hash.new { |h,k| h[k] = [] }

          count = 0

          while line = f.gets
            line.chomp!
            orthogroup_id, score, sp, conf, isoform = line.split("\t")
            gene, num = isoform.split('-').map { |x| x.to_i }

            # Note which gene comes from which species
            @genes_by_species[sp] << gene

            # Don't need to assign new IDs to orthogroups since INPARANOID basically guarantees they'll be 0-counted without
            # skips.
            orthogroup_id = orthogroup_id.to_i

            # If a gene has isoforms that are split between two orthogroups, mark that the orthogroups are merged.
            # Use the first-encountered orthogroup ID for both.
            if @g_to_o.has_key?(gene) && @g_to_o[gene] != orthogroup_id
              @renumber[orthogroup_id] = @renumber[@g_to_o[gene]]
            else
              @renumber[orthogroup_id] = count
              count += 1
            end
            @g_to_o[gene] = orthogroup_id  # Need to keep track of the orthogroup ID so we can continue to renumber.

            # Map renumbered (merged) orthogroup IDs to lists of genes that belong in them.
            @r_to_g[@renumber[orthogroup_id]] << gene unless @r_to_g[@renumber[orthogroup_id]].include?(gene)

            # Save the total number of merged orthogroups.
            @orthogroup_count = count
          end

        end
      end
    end
  end
end