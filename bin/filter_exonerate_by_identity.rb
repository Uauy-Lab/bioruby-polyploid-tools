#!/usr/bin/env ruby
require 'bio'
require 'optparse'
$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
$: << File.expand_path('.')
path= File.expand_path(File.dirname(__FILE__) + '/../lib/bioruby-polyploid-tools.rb')
require path

options = {}
options[:identity] = 95
options[:covered] = 90
OptionParser.new do |opts|
  
  opts.banner = "Usage: filter_exonerate_by_identity.rb [options]"

  opts.on("-e", "--exo FILE", "Exonerate alignment produced by polymarker or with the following ryo: 'RESULT:\\t%S\\t%pi\\t%ql\\t%tl\\t%g\\t%V\\n'") do |o|
    options[:exo_file] = o.upcase
  end
  opts.on("-i", "--identity FLOAT", "Minimum percentage identity") do |o|
    options[:identity] = o.to_f
  end
  opts.on("-c", "--covered FLOAT", "Minimum percentage coverage") do |o|
    options[:covered] = o.to_f
  end
  
end.parse!


exo_file = options[:exo_file]
min_identity = options[:identity];
min_coverage = options[:covered]
File.foreach(exo_file) do |line|  
  aln = Bio::DB::Exonerate::Alignment.parse_custom(line)
  if aln.identity > min_identity and aln.query_coverage > min_coverage
	puts aln.line
  end
end

