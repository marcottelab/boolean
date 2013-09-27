lib = File.expand_path('../lib/', __FILE__)
$:.unshift lib unless $:.include?(lib)


Gem::Specification.new do |gem|
  gem.name = "boolean"
  gem.version = "0.0.1"
  gem.summary = "Not really a gem."
  gem.description = "Just playing around with gemspec."
  gem.homepage = 'http://marcottelab.org'
  gem.authors = ['John Woods']
  gem.email =  ['john.woods@marcottelab.org']

  gem.files         = `git ls-files`.split("\n")
  gem.test_files    = `git ls-files -- {test,spec,features}/*`.split("\n")
  gem.executables   = `git ls-files -- bin/*`.split("\n").map{ |f| File.basename(f) }
  gem.extensions = ['ext/nmatrix/extconf.rb']
  gem.require_paths = ["lib"]

  gem.required_ruby_version = '>= 1.9.2'

  gem.add_dependency 'rdoc', '>=4.0.1'

  gem.add_development_dependency 'rake', '~>0.9'
  gem.add_development_dependency 'bundler'
  gem.add_development_dependency 'rspec'
  gem.add_development_dependency 'rspec-longrun'
  gem.add_development_dependency 'pry', '~>0.9.9'
  gem.add_development_dependency 'guard-rspec', '~>0.7.0'
  gem.add_development_dependency 'rake-compiler', '~>0.8.1'
end

