require File.join(File.dirname(__FILE__), "spec_helper.rb")
#require 'pry'

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

  it "should quickly run a Boolean::Analysis" do
    `rm real`
    a = Boolean::Analysis.new(:op => :|)
  end

  it "should quickly read a DMatrix file 'real'" do
    d = Boolean::DMatrix.read('real')
  end
=begin
  it "should generate a boolean orthogroup-phenotype matrix and shuffle-copy it properly" do
    reader = Boolean.reader(%w{Hs Dr})

    dr_gpm = Boolean.gp_matrix("phenotypes.2.woods", "Dr")
    dr_opm = dr_gpm.opmatrix(reader) # zebrafish orthogroup-phenotype matrix

    #dr_bopm = Boolean::BOPMatrix.new(dr_opm, :|)
    # Test non-independent shuffling
    shuffled= dr_opm.shuffle_rows

    i_offset = 0
    (0...dr_opm.shape[0]).each do |i|
      if dr_opm.is_skippable?(i)
        i_offset += 1
        next
      end
      a = dr_opm.yale_row_as_array(i)
      b = shuffled.yale_row_as_array(i-i_offset)
      require 'pry'; binding.pry if a.size != b.size
      a.size.should equal(b.size)
    end

    # Test independent shuffling
    shuffled= dr_opm.shuffle_each_row

    i_offset = 0
    (0...dr_opm.shape[0]).each do |i|
      if dr_opm.is_skippable?(i)
        i_offset += 1
        next
      end
      a = dr_opm.yale_row_as_array(i)
      b = shuffled.yale_row_as_array(i-i_offset)
      a.size.should equal(b.size)
    end
  end
=end

end
