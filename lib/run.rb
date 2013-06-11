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
distances.write("distances")

# Find the best distance for each
best_d    = Array.new(hs_opm.shape[0], nil)
distances.each_row.with_index do |target,i|
  best_d[i] = target.min
end

# Keep track of the number of times a better score is reported
counter  = Array.new(hs_opm.shape[0], 0)

# Permute 1,000 times
(0...1000).each do |n|
  permute_bopm = dr_bopm.shuffle_rows
  d            = Boolean::DMatrix.new(hs_opm, permute_bopm)

  (0...hs_opm.shape[0]).each do |phi| # for each human phenotype, see how

  end
end
