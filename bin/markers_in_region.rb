#!/usr/bin/env ruby

#This uses the map output from map_markers_to_contigs.rb 
#You need a reference with the name of the contigs, containing the chromosome 
#arm and a list of sequences to map. The algorithm creates a smaller reference 
#file, so the search only spans across the contigs in the region. This should
#allow to use a refined mapping algorithm. 
require 'bio'
require 'optparse'

$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
$: << File.expand_path('.')
path= File.expand_path(File.dirname(__FILE__) + '/../lib/bioruby-polyploid-tools.rb')
require path


def log(msg)
  time=Time.now.strftime("%Y-%m-%d %H:%M:%S.%L")
  puts "#{time}: #{msg}"
end

markers = nil

options = {}
OptionParser.new do |opts|
  
  opts.banner = "Usage: markers_in_region.rb [options]"

  opts.on("-c", "--chromosome CHR", "chromosome (1A, 3B, etc)") do |o|
    options[:chromosome] = o.upcase
  end
  opts.on("-r", "--reference FASTA", "reference with the contigs") do |o|
    options[:reference] = o
  end
  opts.on("-m", "--map CSV", "File with the map and sequence \n Header: INDEX_90K,SNP_ID,SNP_NAME,CHR,COORDINATES_CHR,MAP_ORDER,CHR_ARM,DISTANCE_CM,SEQUENCE") do |o|
    options[:map] = o
  end
  
end.parse!
#reference="/Users/ramirezr/Documents/TGAC/references/Triticum_aestivum.IWGSP1.21.dna_rm.genome.fa"
reference = options[:reference] if options[:reference]
throw raise Exception.new(), "Reference has to be provided" unless reference