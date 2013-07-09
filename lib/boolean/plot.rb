require 'pry'
require 'rubyvis'


class Float
  def exponent
    self.to_s.split("e")[1].to_i
  end
end


module Boolean
  module Plot
    class << self
      # Bin by order of magnitude.
      # Pairs should consist of a p-value and a number of items with that p-value.
      #
      # Returns fractions
      def rebin_by_exponent n, pairs
        bins = Array.new(n+1, 0) # 0.1-1, 0.01-0.1, 0.001-0.01
        denom = 0
        pairs.each do |pair|
          denom += pair[1]
          (-n...0).each do |i|
            if pair[0] < 10**i
              bins[-i] += pair[1]
              break
            elsif i == 0 && pair[0] >= 1
              STDERR.puts "Warning: item fell out of bin range (>1): #{pair[0]}"
            end
          end
        end

        bins.map { |bin| bin / denom.to_f }
      end

      def fig_2b(real, ran)
        min_exp = [real.keys.min.exponent, ran.keys.min.exponent].min
        real_bins = rebin_by_exponent(min_exp.abs, real)
        ran_bins  = rebin_by_exponent(min_exp.abs, ran)
        #binding.pry
        #binding.pry
        #real_denom = 0
        #real_ary = real.keys
        #real.each_pair { |k,v| real_denom += v }

        #ran_denom = 0
        #ran_ary = ran.keys
        #ran.each_pair { |k,v| ran_denom += v}

        #real_a = real.map { |a| [a[0], a[1].quo(real_denom)] }
        #ran_a  = ran.map  { |a| [a[0], a[1].quo(ran_denom)]  }

        log_min = 10**([(real_bins - [0.0]).min, (ran_bins - [0.0]).min].min.exponent-1)

        w = 500
        h = 500
        x = pv.Scale.linear(0, 30).range(0, w)
        y = pv.Scale.log(log_min, 1.0).range(0, h)

        vis = Rubyvis::Panel.new do
          width w
          height h
          bottom 20
          left 60
          right 50
          top 5

          line do
            data real_bins[2...-1]
            left(lambda { x.scale(self.index) })
            bottom(lambda { |d| d == 0 ? 0 : y.scale(d) })
            stroke_style 'blue'
          end

          line do
            data ran_bins[2...-1]
            left(lambda { x.scale(self.index) })
            bottom(lambda { |d| d == 0 ? 0 : y.scale(d) })
            stroke_style 'red'
            stroke_dasharray(3) if self.respond_to?(:stroke_dasharray)
          end

          rule do
            data y.ticks
            bottom y
            visible(lambda { |d| d.to_s =~ /^1/ } )
          end.anchor('left').add(pv.Label).visible(lambda { |d| d<1 && d.to_s =~ /^1/ }).
            text(y.tick_format)

          rule do
            data x.ticks
            left x
            stroke_style(lambda { |d| d != 0 ? "#eee" : "#000" })
          end.anchor('bottom').add(pv.Label).visible(lambda { |d| d.to_i == d && d.to_i % 5 == 0 }).
              text(lambda { |d| d == 0 ? "[1,0.1)" : "[10^-#{d-1},10^-#{d})" })

          #line do
          #  data ran_a
          #   left(lambda   { |d| x.scale(d[0]) })
          #  bottom(lambda { |d| y.scale(d[1]) })
          #  stroke_style 'red'
          #end
        end

        vis.render
        File.open("2b.svg", "w") do |f|
          f.write(vis.to_svg)
        end

        real_bins
      end
    end
  end
end