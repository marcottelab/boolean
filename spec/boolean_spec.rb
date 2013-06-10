require File.join(File.dirname(__FILE__), "spec_helper.rb")
require 'pry'

describe :boolean do
  #it "should generate an OrthoReader" do
  #  Boolean.reader(%w{Hs Mm})
  #end

  #it "should generate a gene-phenotype matrix (GPMatrix)" do
  #  Boolean.gp_matrix("phenotypes.2.mcgary", "Hs")
  #end

  #it "should generate two comparable orthogroup-phenotype matrices (OPMatrix) from two gene-phenotype matrices" do
  #  reader = Boolean.reader(%w{Hs Mm})
  #
  #  hs_gpm = Boolean.gp_matrix("phenotypes.2.mcgary", "Hs")
  #  hs_opm = hs_gpm.opmatrix(reader) # human orthogroup-phenotype matrix
  #
  #  mm_gpm = Boolean.gp_matrix("phenotypes.2.mcgary", "Mm")
  #  mm_opm = mm_gpm.opmatrix(reader) # mouse orthogroup-phenotype matrix
  #
  #  hs_opm.shape[1].should == mm_opm.shape[1]
  #end


  it "should generate a boolean orthogroup-phenotype matrix and shuffle-copy it properly" do
    reader = Boolean.reader(%w{Hs Dr})

    dr_gpm = Boolean.gp_matrix("phenotypes.2.woods", "Dr")
    dr_opm = dr_gpm.opmatrix(reader) # zebrafish orthogroup-phenotype matrix

    dr_bopm = Boolean::BOPMatrix.new(dr_opm, :|)
    shuffled= dr_bopm.shuffle_rows

    (0...dr_bopm.shape[0]).each do |i|
      a = dr_bopm.yale_row_as_array(i)
      b = shuffled.yale_row_as_array(i)
      a.size.should equal(b.size)
    end
  end

end