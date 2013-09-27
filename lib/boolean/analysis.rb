module Boolean
  class Analysis

    def inspect #:nodoc:
      ["#<Boolean::OrthoReader:0x00#{(self.object_id << 1).to_s(16)}", "@to @to_gpm @from @from_gpm", "distances=#{@distances.inspect}", "op=:#{@op.to_s}", "reader=#{@reader.inspect}"].join(' ') + " >"
    end


    def initialize opts = {}
      opts.reverse_merge!({
         :from => ["phenotypes.2.woods", "Dr"],
         :to => ["phenotypes.2.mcgary", "Hs"],
         :op => nil
       })

      STDERR.puts "Initial setup..."
      @reader    = Boolean.reader([opts[:to][1], opts[:from][1]])

      @to_gpm    = Boolean.gp_matrix(*opts[:to])
      @to        = @to_gpm.opmatrix(@reader)

      @from_gpm   = Boolean.gp_matrix(*opts[:from])
      from_opm   = from_gpm.opmatrix(@reader)

      @from      = opts[:op].nil? ? from_opm : BOPMatrix.new(from_opm, opts[:op])

      @distances = if File.exist?("real")
                     Boolean.say_with_time "Reading existing 'real' matrix" do
                       Boolean::DMatrix.read("real")
                     end
                   else
                     d = Boolean.say_with_time "Creating new DMatrix" do
                       Boolean::DMatrix.new(@to, @from)
                     end
                     Boolean.say_with_time "Writing 'real' matrix" do
                       d.write("real", false)
                     end
                   end

      @op = opts[:op]
    end

    attr_reader :to, :to_gpm, :from, :from_gpm, :distances, :reader, :op

    # Get p-value distribution for the real matrix
    def pvalue_distribution
      distrib = RBTree.new { |h,k| h[k] = 0 }
      distances.each do |v|
        distrib[v] += 1
      end
      distrib
    end

    # For i in @to, get the nearest k (default 1) bins in the @from matrix based on @distances
    def binned_nearest i, k: 1, cutoff: 0.0001
      bins = distances.row(i).binned_sorted_indices
      bins = bins[0...[k,bins.size].min]
      within_cutoff = []
      bins.each do |bin|
        break if distances[i,bin[0]] > cutoff
        within_cutoff << bin
      end
      within_cutoff
    end

    def unitary_phenotype_description id, which=:from
      which = which == :from ? :from_gpm : :to_gpm
      phenotype_identifier = self.send(which).phenotypes.select { |k,v| v == id }.keys[0]
      `grep '#{phenotype_identifier}\t' data/PhenotypeDescriptions.*`.rstrip.split("\t").last
    end

    def binary_phenotype_descriptions id
      from.decipher[id].map do |combo|
        left = unitary_phenotype_description(combo.left, :from)
        right = unitary_phenotype_description(combo.right, :from)
        [left,right].join(friendly_op)
      end
    end

    def from_phenotype_description id
      if op.nil?
        from_gpm.phenotypes.select { |k,v| v == id }
      else
      end
    end

    def friendly_op
      if op == :|
        return " **OR** "
      elsif op == :&
        return " **AND** "
      elsif op == :-
        return " **LESS** "
      else
        return " **???** "
      end
    end

    # Display the phenotypes we find that might be good.
    def display_binned_nearest i, k: 1, cutoff: 0.0001
      bins = binned_nearest(i,k,cutoff)
      if bins.empty?
        puts "Nothing found for #{i} with k=#{k} and a cutoff of #{cutoff}"
        return nil
      end

      puts "i=#{i}: #{unitary_phenotype_description(i, :to)}"

      to_set   = to.orthogroups_for_phenotype(i)

      bins.each.with_index do |bin,rank|
        puts "  in #{rank}-ranked bin with distance: #{distances[i,bin[0]]}"
        bin.each do |j|
          from_set = from.orthogroups_for_phenotype(j)
          op_set = to_set & from_set
          puts "  j=#{j}: ( #{to_set.size} | #{op_set.size} | #{from_set.size} )"
          if self.op.nil?
            puts "   - #{unitary_phenotype_description(j)}"
          else
            binary_phenotype_descriptions(j).each do |desc|
              puts "   - #{desc}"
            end
          end
        end
      end

    end


    # Helper function for filtering and displaying all of the likely phenologs.
    def filter_and_display_all_binned_nearest k: 1, cutoff: 0.0001
      (0...to.shape[0]).each do |i|
        to_set   = to.orthogroups_for_phenotype(i)
        next if to_set.size < 3 # Skip for phenotypes that don't have at least 3 items

        bins = binned_nearest(i, k: k, cutoff: cutoff)
        next if bins.empty? # Nothing found for this one.
        bins.each.with_index do |bin,rank|
          bin.each do |j|
            from_set = from.orthogroups_for_phenotype(j)
            op_set   = to_set & from_set
            binding.pry if op_set.size == 0
            next if op_set.size < 2 || from_set.size == op_set.size

            # Display the set information.
            puts "  j=#{j}: ( #{to_set.size} | #{op_set.size} | #{from_set.size} )"
            # Display them.
            if self.op.nil?
              puts "   - #{unitary_phenotype_description(j, :from)}"
            else
              binary_phenotype_descriptions(j, :from).each do |desc|
                puts "   - #{desc}"
              end
            end

          end
        end
      end
    end

  end
end