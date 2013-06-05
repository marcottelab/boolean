require "./lib/boolean.rb"
# require "gsl"
require "distribution"

reader = Boolean.reader(%w{Hs Dr}); puts "Done"

hs_gpm = Boolean.gp_matrix("phenotypes.2.mcgary", "Hs")
hs_opm = hs_gpm.opmatrix(reader) # human orthogroup-phenotype matrix

dr_gpm = Boolean.gp_matrix("phenotypes.2.woods", "Dr")
dr_opm = dr_gpm.opmatrix(reader) # zebrafish orthogroup-phenotype matrix

dr_bopm = Boolean::BOPMatrix.new(dr_opm, :|)

# Initialize dense matrix to 0.0, float64
distances = NMatrix.new([hs_opm.shape[0],dr_bopm.shape[0]], 0.0, :float64)

(0...hs_opm.shape[0]).each do |i|
  m_set = hs_opm.orthogroups_for_phenotype(i)
  STDERR.puts "i=#{i}"
  (0...dr_bopm.shape[0]).each do |j|
    n_set = mm_opm.orthogroups_for_phenotype(j)
    k_set = m_set & n_set
    distances[i,j] = Distribution::Hypergeometric.cdf(k_set.size, m_set.size, n_set.size, hs_opm.shape[1])
  end
end
