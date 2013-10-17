# Extension for dense NMatrix for calculating PPV.

module Boolean
  module PPV
    # Generate a p-value count RBTree (excluding Infinity and NaN).
    # The result can easily be turned into a distribution using
    #     x.normalize
    def counts
      c = RBTree.new { |h,k| h[k] = 0 }
      self.each do |v|
        next if v == Float::INFINITY || v == Float::NAN
        c[v] += 1
      end
      c
    end

    def dist
      counts.normalize
    end

    # Generate an RBTree of PPVs from this distance matrix, using another RBTree as
    # background.
    def print_ppv_table background, threshold = 0.2, permutations = 1000
      #real_val = self.sort.select { |v| v < 1.0 } # sorted p-values below 1.0
      #real_val << 1.0

      true_count  = 0
      false_count = 0
      total_count = 0

      current_real_size = 0

      ppv_count = RBTree.new { |h,k| h[k] = 0 }

      real = self.counts

      background_iter = background.each_relevant_pvalue
      current_background = background_iter.next

      real.each_relevant_pvalue.with_index do |pval_and_end_of_bin, i|
        pval, end_of_bin = pval_and_end_of_bin
        true_count += 1
        break if pval == 1.0
        current_real_size += 1

        if end_of_bin #if pval != real_val[i+1]

          while pval >= current_background[0] # 0 because 1 is end_of_bin
            false_count += 1
            current_background = background_iter.next
          end

          total_count = true_count + false_count / permutations.to_f
          cumulative_fdr = (false_count / permutations.to_f) / total_count
          break if cumulative_fdr > threshold
          ppv = 1.0 - cumulative_fdr

          # Print the details
          puts("#{pval}\t#{true_count}\t#{false_count}\t#{ppv}")

          ppv_count[ppv] += current_real_size
          current_real_size = 0
        end
      end

      num_above_ppv_threshold = 0
      ppv_count.reverse_each do |ppv, count|
        num_above_ppv_threshold += count
        puts "#{ppv}\t#{num_above_ppv_threshold}" # non-monotonic (note: I think this actually means NOT non-monotonic, just copied from Kris' code)
      end
    end
  end
end