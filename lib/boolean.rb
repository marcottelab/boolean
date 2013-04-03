require_relative "boolean/gpmatrix.rb"
require_relative "boolean/opmatrix.rb"
require_relative "boolean/ortho_reader.rb"

module Boolean
  class << self
    # Loads the sqltable file for a given species pair and assigns numbers which will be used as indices
    # for the matrices.
    def reader species=%w{Hs Mm}
      raise(ArgumentError, "Expected pair of species in an array (e.g., %w{Hs Mm})") unless species.is_a?(Array) && species.size == 2
      OrthoReader.new species
    end

    def gp_matrix basename="phenotypes.2.mcgary", species
      GPMatrix.new "data/#{basename}.#{species}", species
    end

    #opm = gpm.opmatrix reader, "Mm"
  end
end

# Each row of opm represents a mouse phenotype