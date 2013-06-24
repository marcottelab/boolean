require "benchmark"
require "distribution"

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
        :n => 1000,
        :with => :shuffle_rows,
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

      from = opts[:op].nil? ? from_opm : BOPMatrix.new(from_opm, opts[:op])

      real = DMatrix.new(to, from)
      real.write("real")

      say_with_time "Permuting #{opts[:n]} times" do
        (0...opts[:n]).each do |n|
          say_with_time "(#{n}/#{opts[:n]})" do
            random_from = from.send(opts[:with])
            random      = DMatrix.new(to, random_from)
            random.write("random.#{n}")
          end
        end
      end

      true
    end

    # Analyze the results of a permutation test. Only argument is +n+, the number of randomizations.
    # Writes to a matrix file called "counts". Returns the real matrix and the counts together.
    def analyze_permutation_test(n)
      real = say_with_time "Reading 'real' matrix" do
        DMatrix.read("real") # rows:to; columns:from
      end

      counts = NMatrix.new(:dense, [real.shape[0], real.shape[1]], 0, :int16)

      (0...n).each do |t|
        say_with_time "Analyzing (#{t}/#{n})" do
          STDERR.puts "\tReading 'random' matrix"
          random = DMatrix.read("random.#{t}")
          STDERR.puts "\tUpdating counts..."
          (0...counts.shape[0]).each do |to_phi|
            (0...counts.shape[1]).each do |from_phi|
              counts[to_phi,from_phi] += 1 if random[to_phi,from_phi] <= real[to_phi,from_phi]
            end
          end
        end
      end

      say_with_time "Writing 'counts' matrix" do
        counts.write("counts")
      end

      [real, counts]
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