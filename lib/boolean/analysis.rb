module Boolean
  class Analysis

    def inspect #:nodoc:
      ["#<Boolean::Analysis:0x00#{(self.object_id << 1).to_s(16)}", "@to @to_gpm @from @from_gpm", "distances=#{@distances.inspect}", "op=:#{@op.to_s}", "reader=#{@reader.inspect}"].join(' ') + " >"
    end


    def initialize opts = {}
      opts.reverse_merge!({
         :from => ["phenotypes.2.woods", "Dr"],
         :to => ["phenotypes.2.mcgary", "Hs"],
         :op => nil,
         :components => true
       })

      STDERR.puts "Initial setup..."
      @reader    = Boolean.reader([opts[:to][1], opts[:from][1]])

      @to_gpm    = Boolean.gp_matrix(*opts[:to])
      @to        = @to_gpm.opmatrix(@reader)
      @to_species= opts[:to][1]

      @from_gpm   = Boolean.gp_matrix(*opts[:from])
      @from_opm   = from_gpm.opmatrix(@reader)

      @op = opts[:op]

      # Get the basic distances (no operation).
      # We only create a Proc here to keep from having to type this twice.
      read_or_create_components = Proc.new { Boolean::Analysis.read_or_create_file(@to, @from_opm, nil) }

      if @op.nil?
        @from      = @from_opm
        @distances = read_or_create_components.call
      else # Get the distances to the boolean combinations in addition.
        @from       = BOPMatrix.new(@from_opm, @op)
        @components = read_or_create_components.call if opts[:components]
        @distances  = Boolean::Analysis.read_or_create_file(@to, @from, @op)
      end

      stat_suffix = @op ? "#{Boolean::Analysis::filename_friendly_op(@op)}.stat.txt" : "stat.txt"
      File.open("#{@to_species}-#{opts[:from][1]}.#{stat_suffix}", "w") do |f|
        f.puts    "OrthoReader:"
        f.puts    "  Entrez ID to inparanoid orthogroup mappings: #{@reader.g_to_o.size}"
        f.puts    "    Unique inparanoid orthogroups: #{@reader.g_to_o.values.sort.uniq.size}"
        f.puts    "  Unique renumbered orthogroups: #{@reader.orthogroup_count}"
        f.puts    "  Inparanoid-observed genes by species:"

        @reader.genes_by_species.each_pair do |sp,set|
          f.puts  "    #{sp}: #{set.size}"
        end

        f.puts    "'To' gene-phenotype matrix:"
        f.puts    "  phenotypes (rows): #{@to_gpm.shape[0]}"
        f.puts    "  genes (columns): #{@to_gpm.shape[1]}"

        f.puts    "'From' gene-phenotype matrix:"
        f.puts    "  phenotypes (rows): #{@from_gpm.shape[0]}"
        f.puts    "  genes (columns): #{@from_gpm.shape[1]}"

        f.puts    "'To' orthogroup-phenotype matrix:"
        f.puts    "  phenotypes (rows): #{@to.shape[0]}"
        f.puts    "  orthogroups (columns): #{@to.shape[1]}"
        f.puts    "  unique phenotypes (rows - skippables): #{@to.shape[0] - @to.count_skippables}"

        f.puts    "'From' orthogroup-phenotype matrix:"
        f.puts    "  phenotypes (rows): #{@from_opm.shape[0]}"
        f.puts    "  orthogroups (columns): #{@from_opm.shape[1]}"
        f.puts    "  unique phenotypes (rows - skippables): #{@from_opm.shape[0] - @from_opm.count_skippables}"

        if @from.is_a?(BOPMatrix)
          f.puts  "'From' boolean orthogroup-phenotype matrix:"
          f.puts  "  unique phenotypes (rows): #{@from.shape[0]}"
          f.puts  "  total combinations: #{@from.decipher.values.flatten.size}"
          f.puts  "  orthogroups (columns): #{@from.shape[1]}"
        end

      end
    end

    def self.read_or_create_file(to, from, op)
      filename = op.nil? ? "real" : "real.#{filename_friendly_op(op)}"
      if File.exist?(filename)
        Boolean.say_with_time "Reading existing '#{filename}' matrix" do
          Boolean::DMatrix.read(filename)
        end
      else
        d = Boolean.say_with_time "Creating new DMatrix" do
          Boolean::DMatrix.new(to, from)
        end
        Boolean.say_with_time "Writing '#{filename}' matrix" do
          d.write(filename, false)
        end
        d
      end
    end

    attr_reader :to, :to_gpm, :from, :from_gpm, :from_opm, :distances, :components, :reader, :op

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
      which_gpm = which == :from ? :from_gpm : :to_gpm
      which_opm = which == :from ? :from_opm : :to
      phenotype_identifiers = self.send(which_opm).equivalent_phenotype_ids(id).map do |eq_id|
        self.send(which_gpm).phenotype_ids[eq_id]
      end
      phenotype_identifiers.map do |external_identifier|
        `grep '#{external_identifier}\t' data/PhenotypeDescriptions.*`.rstrip.split("\t").last
      end.join('[=]')
    end

    def phenotype_description_grep str, which=:to
      which = which == :from ? :from_gpm : :to_gpm
      `grep '#{str}' data/PhenotypeDescriptions.*`.split("\n").map do |g|
        phenotype_id = g.split("\t").first.split(':').last
        self.send(which).phenotypes[phenotype_id]
      end
    end

    def binary_phenotype_descriptions id
      from.decipher[id].map do |combo|
        left = unitary_phenotype_description(combo.left, :from)
        right = unitary_phenotype_description(combo.right, :from)
        [left,right].join(friendly_op)
      end
    end

    def binary_phenotype_distances i, id
      from.decipher[id].map do |combo|
        [components[i,combo.left], components[i,combo.right]]
      end
    end

    #def from_phenotype_description id
    #  if op.nil?
    #    from_gpm.phenotypes.select { |k,v| v == id }
    #  else
    #  end
    #end


    def self.filename_friendly_op op
      if op == :|
        return "OR"
      elsif op == :&
        return "AND"
      elsif op == :-
        return "LESS"
      else
        return "UNKNOWN"
      end
    end

    def friendly_op
      "**#{Boolean::Analysis.filename_friendly_op(op)}**"
    end

    # Display the phenotypes we find that might be good.
    def display_binned_nearest i, k: 1, cutoff: 0.0001
      bins = binned_nearest(i, k: k, cutoff: cutoff)
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


    # For phenotype combination i (to) and j (from), display all gene sets by orthogroup -- candidates, known, and overlap.
    def display_gene_sets i, j
      to_set   = to.orthogroups_for_phenotype(i)
      from_set = from.orthogroups_for_phenotype(j)
      op_set   = to_set & from_set

      hash_orthogroups_to_genes_for_set = Proc.new do |set|
        # There are multiple inparanoid orthogroups for each internal oid, sometimes.
        # This happens because of orthogroups getting merged. So we need to get all
        # of them and squish them together.
        inparanoid_set = set.map { |oid| @reader.unrenumber[oid] }.flatten.sort.uniq

        h = {}
        inparanoid_set.each do |inp_oid|
          entrez_ids = @reader.o_to_g[inp_oid]
          # Only include the Entrez IDs from the species to be predicted.
          h[inp_oid] = entrez_ids.select { |g| @reader.genes_by_species[@to_species].include?(g) }
        end
        h
      end

      to_h      = hash_orthogroups_to_genes_for_set.call(to_set)
      overlap_h = hash_orthogroups_to_genes_for_set.call(op_set)
      from_h    = hash_orthogroups_to_genes_for_set.call(from_set)

      puts "* #{unitary_phenotype_description(i, :to)}"
      puts "  #{to_h}"

      if from.is_a?(BOPMatrix)
        binary_phenotype_descriptions(j).each do |desc|
          puts "* #{desc}"
        end
      else
        puts "* #{unitary_phenotype_description(j, :from)}"
      end

      puts "  #{from_h}"
      puts "* overlap"
      puts "  #{overlap_h}"

      [to_h, overlap_h, from_h]
    end


    # Helper function for filtering and displaying all of the likely phenologs.
    def filter_and_display_all_binned_nearest k: 1, cutoff: 0.0001
      #random_dist = ::Boolean.load_random_permutation_test(op).normalize

      (0...to.shape[0]).each do |i|
        to_set   = to.orthogroups_for_phenotype(i)
        next if to_set.size < 3 # Skip for phenotypes that don't have at least 3 items

        bins = binned_nearest(i, k: k, cutoff: cutoff)
        next if bins.empty? # Nothing found for this one.

        str = "* #{unitary_phenotype_description(i, :to)}"

        bins.each.with_index do |bin,rank|
          bin.each do |j|
            from_set = from.orthogroups_for_phenotype(j)
            op_set   = to_set & from_set
            next if op_set.size < 2 || from_set.size == op_set.size

            #binary_phenotype_distances i, id
            unless str.nil? # Only print the phenotype once, and only do it if stuff has been found.
              puts str
              str = nil
            end

            # Calculate the candidate groups to display. Want to hash from INPARANOID ID to Entrez ID,
            # which may be tricky since we've used so many mappings to coerce the data into the right
            # format. There's a complete explanation in my lab notebook (9/29/13), which I may digitize.
            candidate_groups = begin
              candidate_set  = from_set - op_set

              # There are multiple inparanoid orthogroups for each internal oid, sometimes.
              # This happens because of orthogroups getting merged. So we need to get all
              # of them and squish them together.
              inparanoid_set = candidate_set.map do |oid|
                @reader.unrenumber[oid]
              end.flatten.sort.uniq

              # Now hash from INPARANOID orthogroup ID to Entrez gene ID
              h = {}
              inparanoid_set.each do |inp_oid|
                entrez_ids = @reader.o_to_g[inp_oid]
                # Only include the Entrez IDs from the species to be predicted.
                h[inp_oid] = entrez_ids.select { |g| @reader.genes_by_species[@to_species].include?(g) }
              end
              h
            end

            # Display the set information.
            puts "  j=#{j}: ( #{to_set.size} | #{op_set.size} | #{from_set.size} )\tp=#{'%1.5E' % distances[i,j]}"
            puts "  candidate orthogroups: #{candidate_groups.inspect}"
            # Display them.
            if self.op.nil?
              puts "   - #{unitary_phenotype_description(j, :from)}"
            elsif defined?(@components)
              # Get the individual phenotype distances for the components. We want to make sure
              # the distance for the boolean combination isn't actually worse than the individual
              # components. We only want where the boolean combination is less probable.
              unitary_distances = binary_phenotype_distances(i, j)

              # Write a + or - for each whole number difference in log-pvalues. The first is free. If they're equal, = is written.
              mark = Proc.new do |d|
                log_difference = Math.log(d / distances[i,j])
                if log_difference < 0
                  '-' * (1-log_difference)
                elsif log_difference > 0
                  '+' * (1+log_difference)
                else
                  '='
                end
              end

              binary_phenotype_descriptions(j).each.with_index do |desc,idx|
                d1, d2 = unitary_distances[idx]
                puts "   - #{desc} : #{'%1.2E' % d1} : #{'%1.2E' % d2}\t#{mark.call(d1)} : #{mark.call(d2)} "
              end
            else # No filtering -- just print them all.
              binary_phenotype_descriptions(j).each.with_index do |desc,idx|
                puts "   - #{desc}"
              end
            end

          end
        end
      end
    end

  end
end
