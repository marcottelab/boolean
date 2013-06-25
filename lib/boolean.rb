require "benchmark"
require "distribution"
require "rbtree" # Use as an ordered hash

require_relative "boolean/gpmatrix.rb"
require_relative "boolean/opmatrix.rb"
require_relative "boolean/bopmatrix.rb"
require_relative "boolean/ortho_reader.rb"
require_relative "boolean/dmatrix.rb"
require_relative "hypergeometric.rb"

class Hash
  # File activesupport/lib/active_support/core_ext/hash/reverse_merge.rb, line 17
  def reverse_merge!(other_hash)
    # right wins if there is no left
    merge!( other_hash ){|key,left,right| left }
  end
end


module Boolean
  class << self

    # Loads the sqltable file for a given species pair and assigns numbers which will be used as indices
    # for the matrices.
    def reader species=%w{Hs Mm}
      raise(ArgumentError, "Expected pair of species in an array (e.g., %w{Hs Mm})") unless species.is_a?(Array) && species.size == 2
      OrthoReader.new species
    end

    def gp_matrix basename="phenotypes.2.mcgary", species
      GPMatrix.new "data/#{basename}.#{species}", species
    end

    def permutation_test(opts = {})
      opts.reverse_merge!({
        :start => 0,
        :end => 1000,
        :with => :shuffle_each_row, # or :shuffle_rows (doesn't seem to work right)
        :from => ["phenotypes.2.woods", "Dr"],
        :to => ["phenotypes.2.mcgary", "Hs"],
        :op => nil
      })

      STDERR.puts "Initial setup..."
      reader   = reader([opts[:to][1], opts[:from][1]])

      to_gpm   = gp_matrix(*opts[:to])
      to       = to_gpm.opmatrix(reader)

      from_gpm = gp_matrix(*opts[:from])
      from_opm = from_gpm.opmatrix(reader)

      from     = opts[:op].nil? ? from_opm : BOPMatrix.new(from_opm, opts[:op])

      real     = DMatrix.new(to, from)
      real.write("real")

      start_i  = opts[:start]
      end_i    = opts[:end]

      say_with_time "Permuting #{opts[:end]} times" do
        (start_i...end_i).each do |i|
          say_with_time "(#{i}/#{end_i-start_i})" do
            random_from = from.send(opts[:with])
            random      = DMatrix.new(to, random_from)
            random.write("random.#{i}")
          end
        end
      end

      true
    end

    # Analyze the results of a permutation test. Only argument is +n+, the number of randomizations.
    # Writes to a matrix file called "counts". Returns the real matrix and the counts together.
    def analyze_permutation_test(n)
      real_matrix = say_with_time "Reading 'real' matrix" do
        DMatrix.read("real") # rows:to; columns:from
      end

      # Track the distributions by p-value
      real_dist   = RBTree.new   { |h,k| h[k] = 0 }
      random_dist = RBTree.new   { |h,k| h[k] = 0 }

      real_matrix.each do |v|
        real_dist[v]  += 1
      end

      #counts = NMatrix.new(:dense, [real.shape[0], real.shape[1]], 0, :int16)

      (0...n).each do |t|
        say_with_time "Analyzing (#{t}/#{n})" do
          STDERR.puts "\tReading 'random' matrix"
          random_matrix = DMatrix.read("random.#{t}")
          STDERR.puts "\tUpdating counts..."
          (0...real_matrix.shape[0]).each do |to_phi|
            (0...real_matrix.shape[1]).each do |from_phi|

              random_dist[ random_matrix[to_phi,from_phi] ] += 1

            end
          end
        end
      end

      [real_dist, random_dist]
    end

    def say_with_time msg
      puts msg
      result = nil
      time = Benchmark.measure { result = yield }
      puts "%.4fs" % time.real
      result
    end

    #opm = gpm.opmatrix reader, "Mm"
  end
end

# Each row of opm represents a mouse phenotype