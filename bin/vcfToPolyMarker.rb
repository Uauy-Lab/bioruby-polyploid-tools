#/usr/bin/env ruby
require 'optparse'

require 'csv'
require 'bio'
require 'bio-samtools'


OptionParser.new do |opts|
  opts.banner = "Usage: polymarker.rb [options]"

  opts.on("-c", "--contigs FILE", "File with contigs to use as database") do |o|
    options[:path_to_contigs] = o
  end
  
  opts.on("-m", "--marker_list FILE", "File with the list of markers to search from") do |o|
    options[:marker_list] = o
  end

  
end