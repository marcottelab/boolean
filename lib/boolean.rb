require "benchmark"
require "distribution"
require "rbtree" # Use as an ordered hash
require 'yaml'

require_relative "boolean/analysis.rb"
require_relative "boolean/gpmatrix.rb"
require_relative "boolean/opmatrix.rb"
require_relative "boolean/bopmatrix.rb"
require_relative "boolean/ortho_reader.rb"
require_relative "boolean/dmatrix.rb"
require_relative "boolean/plot.rb"
require_relative "hypergeometric.rb"

class Hash
  # File activesupport/lib/active_support/core_ext/hash/reverse_merge.rb, line 17
  def reverse_merge!(other_hash)
    # right wins if there is no left
    merge!( other_hash ){|key,left,right| left }
  end
end

unless File.respond_to?(:write)
  class File
    def self.write filename, obj
      File.open(filename, 'w') do |f|
        f.write(obj.to_yaml)
      end
    end
  end
end

# Modified from: http://markmail.org/message/ntlhxtacsqjh7rd4 (thanks to Why the lucky stiff!)
# Ended up just basing it mostly on Hash though.
class RBTree
  YAML.add_ruby_type /RBTree/ do |type,val|
    r = RBTree.new
    val.each { |k,v| r[k] = v }
    r
  end

  def to_yaml( opts = {} )
    YAML::quick_emit( self, opts ) do |out|
      out.map( taguri, to_yaml_style ) do |map|
        each do |k,v|
          map.add(k.to_s, v)
        end
      end
    end
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

    def parallel_permutation_test(to, from, real, opts)
      start_i  = opts[:start]
      end_i    = opts[:end]

      random_dist = RBTree.new   { |h,k| h[k] = 0 }
      say_with_time "Permuting #{end_i-start_i} times" do
        (start_i...end_i).each do |i|
          if File.exists?("random.#{i}.gz")
            puts "Iteration #{i} already appears to exist; skipping."
          else
            say_with_time "(#{i}/#{end_i})" do
              random_from = from.send(opts[:with])
              random      = DMatrix.new(to, random_from)

              # Write the file
              random_filename = "random.#{i}"
              random.write(random_filename, :compress)

              random.each do |v|
                random_dist[v] += 1
              end
            end
          end
        end
      end

      say_with_time "Writing random dist to file" do
        File.write("random.dist.yml.#{end_i}", random_dist.to_yaml)
      end
    end

    # Deprecated. Use parallel_permutation_test
    def permutation_test(opts = {})
      opts.reverse_merge!({
        :start => 0,
        :end => 1000,
        :with => :shuffle_each_row, # or :shuffle_rows
      })


      analysis = say_with_time "Creating Boolean::Analysis object" do
        Boolean::Analysis.new(opts)
      end

      to       = analysis.to
      from     = analysis.from
      real     = analysis.distances

      real_dist   = analysis.pvalue_distribution
      random_dist = RBTree.new   { |h,k| h[k] = 0 }


      start_i  = opts[:start]
      end_i    = opts[:end]

      say_with_time "Permuting #{opts[:end]} times" do
        (start_i...end_i).each do |i|
          say_with_time "(#{i}/#{end_i})" do
            random_from = from.send(opts[:with])
            random      = DMatrix.new(to, random_from)

            random.each do |v|
              random_dist[v] += 1
            end

            # Write the file
            random_filename = "random.#{i}"
            random.write(random_filename, :compress)
          end
        end
      end

      say_with_time "Writing distributions to files" do
        File.write("real.dist.yml", real_dist.to_yaml)
        File.write("random.dist.yml", random_dist.to_yaml)
      end

      [real_dist, random_dist]
    end

    # Analyze the results of a permutation test. Only argument is +n+, the number of randomizations.
    # Writes to a matrix file called "counts". Returns the real matrix and the counts together.
    def analyze_permutation_test(n)
      real_matrix = say_with_time "Reading 'real' matrix" do
        #`gunzip -c real.gz > real`
        x = DMatrix.read("real") # rows:to; columns:from
        x
        #`rm real`
      end

      # Track the distributions by p-value
      real_dist   = RBTree.new   { |h,k| h[k] = 0 }
      random_dist = RBTree.new   { |h,k| h[k] = 0 }

      real_matrix.each do |v|
        real_dist[v]  += 1
      end

      # Remove skipped values from the distribution.
      real_dist.delete(Float::INFINITY)

      #counts = NMatrix.new(:dense, [real.shape[0], real.shape[1]], 0, :int16)

      (0...n).each do |t|
        say_with_time "Analyzing (#{t}/#{n})" do
          STDERR.puts "\tReading 'random' matrix"
          random_matrix = DMatrix.read("random.#{t}")
          STDERR.puts "\tUpdating counts..."
          random_matrix.each do |v|
            random_dist[v] += 1
          end
        end
      end

      say_with_time "Writing distributions to files" do
        File.write('real.dist.yml', real_dist.to_yaml)
        File.write('random.dist.yml', random_dist.to_yaml)
      end

      [real_dist, random_dist]
    end

    def load_permutation_test
      [YAML::load(File.read('real.dist.yml')),
       YAML::load(File.read('random.dist.yml'))]
    end

    def merge_permutation_tests
      merged = RBTree.new { |h,k| h[k] = 0 }
      Dir::glob("random.dist.yml.*").each do |filename|
        current = YAML::load(File.read(filename))
        current.each_pair do |pvalue, count|
          merged[pvalue] += count
        end
      end
      File.write("random.dist.yml", merged.to_yaml)
      merged
    end

    def plot_permutation_test(*args)
      real = nil
      ran = nil

      if args.size > 1
        real = args.shift
        ran = args.shift
      else
        real, ran = analyze_permutation_test(args.shift)
      end

      return Boolean::Plot.fig_2b(real, ran)
    end

    def analyze_best_hits(opts = {})
      opts.reverse_merge!({
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

      real     = File.exist?("real") ? DMatrix.read("real") : DMatrix.new(to, from).tap { |r| r.write("real", false) }

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
