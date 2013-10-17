require "benchmark"
require "distribution"
require "rbtree" # Use as an ordered hash
require "psych"

require_relative "boolean/analysis.rb"
require_relative "boolean/gpmatrix.rb"
require_relative "boolean/opmatrix.rb"
require_relative "boolean/bopmatrix.rb"
require_relative "boolean/ortho_reader.rb"
require_relative "boolean/dmatrix.rb"
require_relative "boolean/plot.rb"
require_relative "boolean/ppv.rb"
require_relative "hypergeometric.rb"
require_relative "rbtree_monkey.rb"
require_relative "psych_monkey.rb"

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

      op_piece = Boolean::Analysis.filename_friendly_op(opts[:op]) if opts[:op]

      random_dist = RBTree.new   { |h,k| h[k] = 0 }
      say_with_time "Permuting #{end_i-start_i} times" do
        (start_i...end_i).each do |i|
          filename = opts[:op] ? "random.#{op_piece}.#{i}" : "random.#{i}"
          if File.exists?("#{filename}.gz")
            puts "Iteration #{i} already appears to exist (#{filename}); skipping."
          else
            say_with_time "(#{i}/#{end_i})" do
              random_from = from.send(opts[:with])
              random      = DMatrix.new(to, random_from)

              # Write the file
              random.write(filename, :compress)

              random.each do |v|
                random_dist[v] += 1
              end
            end
          end
        end
      end

      filename = opts[:op] ? "random.#{op_piece}.dist.yml.#{end_i}" : "random.dist.yml.#{end_i}"
      say_with_time "Writing random dist to file '#{filename}'" do
        File.write(filename, random_dist.to_yaml)
      end
    end


    # Merge multiple permutation tests together.
    def merge_permutation_tests op=nil
      file_prefix = op.nil? ? "random.dist.yml" : "random.#{Boolean::Analysis::filename_friendly_op(op)}.dist.yml"

      merged = RBTree.new { |h,k| h[k] = 0 }
      Dir::glob("#{file_prefix}.*").each do |filename|
        say_with_time "Merging file '#{filename}'" do
          current = Psych::load(File.read(filename))
          current.each_pair do |pvalue, count|
            merged[pvalue] += count
          end
        end
      end

      say_with_time "Writing final output, '#{file_prefix}'" do
        File.write(file_prefix, merged.to_yaml)
      end

      merged
    end


    # Analyze the results of a permutation test. Only argument is +n+, the number of randomizations.
    # Writes to a matrix file called "counts". Returns the real matrix and the counts together.
    def analyze_permutation_test(opts)
      real_dist   = load_real_as_dist(opts)
      random_dist = load_random_permutation_test(opts[:op])

      [real_dist, random_dist]
    end

    def load_real_as_dist(opts)
      force = opts[:force]
      suffix = opts[:op] ? ".#{Boolean::Analysis::filename_friendly_op(opts[:op])}" : ""
      d = DMatrix.read("real#{suffix}", force)
      return d.counts if d.respond_to?(:counts)
      c = RBTree.new { |h,k| h[k] = 0 }
      d.each do |v|
        next if v == Float::INFINITY || v == Float::NAN
        c[v] += 1
      end
      c
    end


    def load_permutation_test(opts)
      suffix = opts[:op] ? "#{Boolean::Analysis::filename_friendly_op(opts[:op])}.dist.yml" : "dist.yml"
      [Psych::load(File.read("real.#{suffix}")),
       Psych::load(File.read("random.#{suffix}"))]
    end


    def load_random_permutation_test(op=nil)
      filename = op ? "random.#{Boolean::Analysis::filename_friendly_op(op)}.dist.yml" : "random.dist.yml"
      r = Psych::load(File.read(filename))
      r.delete(Float::INFINITY)
      r
    end


    def plot_permutation_test(opts)
      real, ran = analyze_permutation_test(opts)

      return Boolean::Plot.fig_2b(real, ran)
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
