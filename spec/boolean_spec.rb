require File.join(File.dirname(__FILE__), "spec_helper.rb")
require 'pry'

describe :boolean do
  it "should generate an OrthoReader" do
    Boolean.reader(%w{Hs Mm})
  end

  it "should generate a gene-phenotype matrix (GPMatrix)" do
    Boolean.gp_matrix("phenotypes.2.mcgary", "Hs")
  end

  it "should generate an orthogroup-phenotype matrix (OPMatrix) from a gene-phenotype matrix" do
    reader = Boolean.reader(%w{Hs Mm})
    gpm    = Boolean.gp_matrix("phenotypes.2.mcgary", "Hs")
    opm    = gpm.opmatrix(reader, "Mm")
  end

end