#!/usr/bin/env ruby
require 'bio'
$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
$: << File.expand_path('.')
path= File.expand_path(File.dirname(__FILE__) + '/../lib/bioruby-polyploid-tools.rb')
require path

blat_file=ARGV[0]

blat_aln = Bio::Blat::Report.new(Bio::FlatFile.open(blat_file).to_io)
blat_aln.each_hit() do |hit|
  if hit.percentage_covered >= 50
    puts hit.data.join("\t")
  end
end