# Orthogroup Reader
#
# This class parses the sqltable file (INPARANOID output) for a given species pair. It is needed for converting a
# gene-phenotype matrix (GPMatrix) to an orthogroup-phenotype matrix (OPMatrix).
class OrthoReader

  attr_reader :g_to_o, :r_to_g, :renumber, :orthogroup_count

  # Convert a (renumbered) gene ID (a GPMatrix index) to our renumbered orthogroup ID (an OPMatrix index)
  def orthogroup_id gid
    @renumber[@g_to_o[gid]]
  end

  # Convert a renumbered orthogroup ID (OPMatrix index) to a renumbered gene ID (GPMatrix index).
  def gene_id rid
    @r_to_g[rid]
  end

  # Create an OrthoReader for a pair of species (given by +species_ids+). Optionally, can supply the data directory as
  # second argument, which defaults to "data"
  def initialize species_ids, dir="data"
    @renumber = {} # map INPARANOID orthogroup IDs to renumbered oids

    Dir.chdir(dir) do
      possible_filenames = ["sqltable.#{species_ids[0]}-#{species_ids[1]}", "sqltable.#{species_ids[1]}-#{species_ids[0]}"]
      filename = File.exists?(possible_filenames.first) ? possible_filenames.first : possible_filenames.last
      f = File.new(filename, "r")

      # genes to orthogroups and orthogroups to genes
      @g_to_o = {}
      o_to_g = Hash.new { |h,k| h[k] = [] }
      @r_to_g = Hash.new { |h,k| h[k] = [] }

      count = 0

      while line = f.gets
        line.chomp!
        orthogroup_id, score, sp, conf, isoform = line.split("\t")
        gene, num = isoform.split('-').map { |x| x.to_i }

        # Don't need to assign new IDs to orthogroups since INPARANOID basically guarantees they'll be 0-counted without
        # skips.
        orthogroup_id = orthogroup_id.to_i

        if @g_to_o.has_key?(gene) && @g_to_o[gene] != orthogroup_id
          @renumber[orthogroup_id] = @renumber[@g_to_o[gene]]
        else
          @renumber[orthogroup_id] = count
          count += 1
        end

        @g_to_o[gene] = orthogroup_id
        o_to_g[orthogroup_id] << gene unless o_to_g[orthogroup_id].include?(gene)
        @r_to_g[@renumber[orthogroup_id]] << gene unless @r_to_g[@renumber[orthogroup_id]].include?(gene)
        @orthogroup_count = count
      end

      f.close
    end
  end
end