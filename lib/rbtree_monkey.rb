# These are monkey patches for RBTree, mostly for YAML saving and loading (see also psych_monkey.rb).
# There are also some functions for calculating FDRs, which is relatively easy to do with an RBTree.
#
# Modified from: http://markmail.org/message/ntlhxtacsqjh7rd4 (thanks to Why the lucky stiff!)
# Ended up just basing it mostly on Hash though.
class RBTree

  # This part doesn't work. See psych.rb.
  #YAML.add_ruby_type /RBTree/ do |type,val|
  #  r = RBTree.new
  #  val.each { |k,v| r[k] = v }
  #  r
  #end

  def to_yaml( opts = {} )
    YAML::quick_emit( self, opts ) do |out|
      out.map( taguri, to_yaml_style ) do |map|
        each do |k,v|
          map.add(k.to_s, v)
        end
      end
    end
  end


  # Iterate across every p-value in the tree (which isn't 1 or infinity) as if it
  # weren't binned. That is, if 0.05 is in a bin of size 5, iterate across 0.05
  # five times. Also yields a true-false value indicating if this is the last entry
  # in the bin.
  def each_relevant_pvalue
    return enum_for(:each_relevant_pvalue) unless block_given?

    self.each_pair do |p, bin_size|
      next if p >= 1.0
      bin_size.times do |bin_idx|
        yield(p, bin_idx+1 == bin_size)
      end
    end

    self
  end


  def fdr pvalue
    self.lower_bound(pvalue)[1]
  end


  # Normalize to produce FDRs for each bin.
  def normalize
    total = 0

    result = RBTree.new

    self.each_pair do |k,v|
      next if k == Float::INFINITY || k == Float::NAN
      total += v
    end
    total = total.to_f

    sum = 0
    self.each_pair do |k,v|
      next if k == Float::INFINITY || k == Float::NAN
      sum += v
      result[k] = sum / total
    end

    result
  end
end