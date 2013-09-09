require 'rubygems'
require 'bundler'


#begin
#  Bundler.setup(:default, :development)
#rescue Bundler::BundlerError => e
#  $stderr.puts e.message
#  $stderr.puts "Run `bundle install` to install missing gems"
#  exit e.status_code
#end
require 'rake'

#require 'jeweler'
#Jeweler::Tasks.new do |gem|
  # gem is a Gem::Specification... see http://docs.rubygems.org/read/chapter/20 for more options
#  gem.name = "bio-polyploid-tools"
#  gem.homepage = "http://github.com/homonecloco/bioruby-polyploid-tools"
#  gem.license = "MIT"
#  gem.summary = %Q{Tool to work with polyploids, NGS and molecular biology}
#  gem.description = %Q{Repository of tools developed in TGAC and Crop Genetics in JIC to work with polyploid wheat}
#  gem.email = "ricardo.ramirez-gonzalez@tgac.ac.uk"
#  gem.authors = ["Ricardo Ramirez-Gonzalez"]
  # Include your dependencies below. Runtime dependencies are required when using your gem,
  # and development dependencies are only needed for development (ie running rake tasks, tests, etc)
  #  gem.add_runtime_dependency 'jabber4r', '> 0.1'
  #  gem.add_development_dependency 'rspec', '> 1.2.3'
 # gem.extensions = "ext/mkrf_conf.rb"
#end
#Jeweler::RubygemsDotOrgTasks.new

require 'rake/testtask'
Rake::TestTask.new(:test) do |test|
  test.libs << 'lib' << 'test'
  test.pattern = 'test/**/test_*.rb'
  test.verbose = true
end


if RUBY_VERSION.start_with?("1.8")
  require 'rcov/rcovtask'
  Rcov::RcovTask.new do |test|
    test.libs << 'test'
    test.pattern = 'test/**/test_*.rb'
    test.verbose = true
  end
end

task :default => :test

#require 'rdoc/task'
##RDoc::Task.new do |rdoc|
#  version = File.exist?('VERSION') ? File.read('VERSION') : ""

#  rdoc.rdoc_dir = 'rdoc'
#  rdoc.title = "bio-samtools #{version}"
#  rdoc.rdoc_files.include('README*')
#  rdoc.rdoc_files.include('lib/**/*.rb')
#end
