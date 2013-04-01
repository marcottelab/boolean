require "./gpmatrix.rb"
require "./opmatrix.rb"

# Loads the sqltable file for a given species pair and assigns numbers which will be used as indices
# for the matrices.
reader = OrthoReader.new(["Hs", "Mm"])

gpm = GPMatrix.new "data/phenotypes.2.mcgary.Hs", "Hs"
opm = gpm.opmatrix reader, "Mm"



# Each row of opm represents a mouse phenotype