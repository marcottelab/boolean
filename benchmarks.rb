require 'set'
require 'benchmark'

f = 10_000
ar1 = (1..(10*f)).to_a # 100_000 elements
ar2 = ((5*f)..(15*f)).to_a # also 100_000 elements
set1 = ar1.to_set
set2 = ar2.to_set
sset1 = SortedSet.new(ar1)
sset2 = SortedSet.new(ar2)
n = 10 #20000

Benchmark.bm(10) do |testcase|
  testcase.report('Array'){ n.times{ ar1 & ar2 } }
  testcase.report('Set'){ n.times{ set1 & set2 } }
  testcase.report('SortedSet') { n.times{ sset1 & sset2 } }
  testcase.report('Set2'){ n.times{ ar1.select{ |element| set2.include? element } } }
  testcase.report('Set2present'){ n.times{ ar1.any?{ |element| set2.include? element } } }

  testcase.report('SortedSet2'){ n.times{ ar1.select{ |element| sset2.include? element } } }
  testcase.report('SortedSet2present'){ n.times{ ar1.any?{ |element| sset2.include? element } } }
end


require './lib/hypergeometric'

Benchmark.bm(10) do |testcase|
  testcase.report('Native') { n.times { Hypergeometric.cdf(2,3,5,3000) }}
  testcase.report('Ruby') { n.times { Hypergeometric.ruby_cdf(2,3,5,3000) }}
  testcase.report('GSL') { n.times { 1.0 - Distribution::Hypergeometric.cdf(1,3,5,3000)}}
end

n = 10000
Benchmark.bm(10) do |testcase|
  testcase.report('GSL') { n.times { GSL::Sf.lnfact(rand(100000)) }}
  testcase.report('Native') { n.times { Distribution::MathExtension::ApproxFactorial.stieltjes_ln_factorial(rand(100000)) }}
end


