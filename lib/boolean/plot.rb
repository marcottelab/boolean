require 'pry'
require 'rubyvis'


class Float
  def exponent
    return nil if self == Float::INFINITY || self == -Float::INFINITY
    return 0 if self == 0
    Math.log10(self).floor
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
          ((-n-1)..0).each do |i|
            if pair[0] <= 10**i
              bins[n+i-1] += pair[1]
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

        log_min = 10**([(real_bins - [0.0]).min, (ran_bins - [0.0]).min].min.exponent-1)

        w = 500
        h = 400
        x = pv.Scale.linear(0, real_bins.size).range(0, w)
        y = pv.Scale.log(log_min, 1.0).range(0, h)

        vis = Rubyvis::Panel.new do
          width w
          height h
          bottom 80
          left 60
          right 20
          top 5

          rule do
            line_width 2
            data y.ticks
            bottom y
            visible(lambda { |d| d.to_s =~ /^1/ } )
          end.anchor('left').add(pv.Label).visible(lambda { |d| d<1 && d.to_s =~ /^1/ }).
            text(y.tick_format)

          rule do
            line_width 2
            data( (0..(real_bins.size)).to_a)
            left x
            stroke_style(lambda { |d| d != 0 ? "#eee" : "#000" })
            visible(lambda { |d| d.to_i == d && d.to_i % 5 == 1 })
          end.anchor('bottom').add(pv.Label).visible(lambda { |d| d.to_i == d && d.to_i % 5 == 1 }).
              text(lambda { |d| d == 1 ? "[1,0.1)" : "[10^-#{d-1},10^-#{d})" }).text_angle(-Math::PI/2).text_align('right')

          line do
            data real_bins.reverse
            line_width 3
            left(lambda { x.scale(self.index) })
            bottom(lambda { |d| d == 0 ? 0 : y.scale(d) })
            stroke_style 'blue'
          end

          line do
            data ran_bins.reverse
            line_width 3
            left(lambda { x.scale(self.index) })
            bottom(lambda { |d| d == 0 ? 0 : y.scale(d) })
            stroke_style 'red'
            stroke_dasharray(3) #if self.respond_to?(:stroke_dasharray)
          end
        end

        vis.render
        File.open("2b.svg", "w") do |f|
          f.write(vis.to_svg)
        end

        [real_bins, ran_bins]
      end
    end
  end
end
