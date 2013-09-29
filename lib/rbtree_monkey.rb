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