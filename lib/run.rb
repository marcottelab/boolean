require "./lib/boolean.rb"
# require "gsl"
require "distribution"

reader = Boolean.reader(%w{Hs Dr}); puts "Done"

hs_gpm = Boolean.gp_matrix("phenotypes.2.mcgary", "Hs")
hs_opm = hs_gpm.opmatrix(reader) # human orthogroup-phenotype matrix

dr_gpm = Boolean.gp_matrix("phenotypes.2.woods", "Dr")
dr_opm = dr_gpm.opmatrix(reader) # zebrafish orthogroup-phenotype matrix

dr_bopm = Boolean::BOPMatrix.new(dr_opm, :|)

# Create dense matrix, don't initialize values by offering a default -- we're going to do so momentarily.
distances = Boolean::DMatrix.new(hs_opm, dr_bopm)

