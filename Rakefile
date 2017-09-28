require 'rubygems'
require 'bundler'
#
#require 'bundler/version'

begin
  Bundler.setup(:default, :development)
  rescue Bundler::BundlerError => e
  $stderr.puts e.message
  $stderr.puts "Run `bundle install` to install missing gems"
  exit e.status_code
end
require 'rake'


if RUBY_VERSION.start_with?("2.1") or RUBY_VERSION.start_with?("2.2") or RUBY_VERSION.start_with?("2.0")
  require 'jeweler'
  @taskClass = Jeweler
else
  require 'juwelier'
  @taskClass = Juwelier
end



@taskClass::Tasks.new do |gem|
  # gem is a Gem::Specification... see http://docs.rubygems.org/read/chapter/20 for more options
   gem.name = "bio-polyploid-tools"
  gem.homepage = "http://github.com/tgac/bioruby-polyploid-tools"
  gem.license = "MIT"
  gem.summary = %Q{Tool to work with polyploids, NGS and molecular biology}
  gem.description = %Q{Repository of tools developed in TGAC and Crop Genetics in JIC to work with polyploid wheat}
   gem.email = "ricardo.ramirez-gonzalez@tgac.ac.uk"
  gem.authors = ["Ricardo H.  Ramirez-Gonzalez"]
  # Include your dependencies below. Runtime dependencies are required when using your gem,
  # and development dependencies are only needed for development (ie running rake tasks, tests, etc)
  #gem.add_runtime_dependency 'bio-samtools', '= 0.6.2'
  #  gem.add_development_dependency 'rspec', '> 1.2.3'
#  gem.extensions = "ext/mkrf_conf.rb"
end
@taskClass::RubygemsDotOrgTasks.new

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

