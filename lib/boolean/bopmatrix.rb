require "benchmark"
require "ruby-progressbar"

def say_with_time msg
  puts msg
  result = nil
  time = Benchmark.measure { result = yield }
  puts "%.4fs" % time.real
  result
end


module Boolean
  class BOPMatrix < OPMatrix
    # Create a new BOPMatrix, which stands for boolean orthogroup-phenotype matrix.
    #
    # It's created from another OPMatrix, but generates phenotypes by producing boolean combinations of the phenotypes
    # from the input matrix.
    #
    # operation is the type of boolean operation to use.
    def initialize opmatrix, operation
      associations = []
      @decipher = {} # convert from count to the boolean combination
      count = 0
      matrix_size = 0

      comb = say_with_time("generating numeric combinations") do
        (0...opmatrix.shape[0]).to_a.combination(2)
      end

      say_with_time("generating gene set combinations") do
        # Generate ID combinations of all phenotypes
        comb.each do |pair|
          left        = opmatrix.orthogroups_for_phenotype(pair[0])
          right       = opmatrix.orthogroups_for_phenotype(pair[1])
          next if left.size == 0 || right.size == 0

          left_right  = left.send(operation, right)
          if left_right != left && left_right != right && left_right.size > 2 # no point in duplicating existing phenotypes
            associations[count] = left_right
            @decipher[count]     = "#{pair[0]}#{operation}#{pair[1]}"
            count += 1
            matrix_size += left_right.size
          end

          next unless operation == :-

          right_left = right.send(operation, left)
          if right_left != right && right_left != right && right_left.size > 2
            associations[count] = right_left
            @decipher[count]     = "#{pair[1]}#{operation}#{pair[0]}"
            count += 1
            matrix_size += left_right.size
          end
        end
      end

      # Create the data structure.
      say_with_time("generating matrix data structure of capacity #{matrix_size}") do
        super(count, opmatrix.shape[1], matrix_size)
        self.extend(NMatrix::YaleFunctions)  # enable yale_vector_insert
        STDERR.puts "ija size=#{self.yale_ija.size}"
      end

      say_with_time("filling matrix data structure") do
        i = 0
        bar = ProgressBar.create(:title => "Fill", :total => matrix_size)
        last_capacity = self.capacity

        begin
          j_array = associations.shift

          # Insert the whole array at once using a specially exposed helper function.
          self.yale_vector_insert(i, j_array.to_a, [1]*j_array.size)
          #j_array.each do |j|
          #  self[i,j] = 1
          #end

          i            += 1
          bar.progress += j_array.size
        end while i < count

        bar.finish
      end
    end

    attr_reader :decipher

  end
end