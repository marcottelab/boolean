require "./lib/boolean.rb"
# require "gsl"
require "distribution"

reader = Boolean.reader(%w{Hs Mm})

hs_gpm = Boolean.gp_matrix("phenotypes.2.mcgary", "Hs")
hs_opm = hs_gpm.opmatrix(reader) # human orthogroup-phenotype matrix

mm_gpm = Boolean.gp_matrix("phenotypes.2.mcgary", "Mm")
mm_opm = mm_gpm.opmatrix(reader) # mouse orthogroup-phenotype matrix

distances = NMatrix.new(:dense, [hs_opm.shape[0],mm_opm.shape[0]], 0, :float64)

(0...hs_opm.shape[0]).each do |i|
  m_set = hs_opm.orthogroups_for_phenotype(i)
  (0...mm_opm.shape[0]).each do |j|
    n_set = mm_opm.orthogroups_for_phenotype(j)
    k_set = m_set & n_set
    distances[i,j] = Distribution::Hypergeometric.cdf(k_set.size, m_set.size, n_set.size, hs_opm.shape[1])
    STDERR.puts "#{i}, #{j}: #{distances[i,j]}"
  end
end
