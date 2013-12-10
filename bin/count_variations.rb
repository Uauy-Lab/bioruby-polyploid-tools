#!/usr/bin/env ruby

require 'bio'
require 'rubygems'
require 'pathname'
require 'bio-samtools'

require 'set'

$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
$: << File.expand_path('.')
path= File.expand_path(File.dirname(__FILE__) + '/../lib/bioruby-polyploid-tools.rb')
require path

puts  ARGV[0]

fasta_db = Bio::DB::Fasta::FastaFile.new( ARGV[0])
fasta_db.load_fai_entries
bam1 =  Bio::DB::Sam.new({:fasta=>ARGV[0], :bam=>ARGV[1]})

fasta_db.index.entries.each do | r |
  #Np r.get_full_region
  #container.process_region( { :region => r.get_full_region.to_s, :output_file => output_file } )
  region=r.get_full_region
  
  

  cons_1 = bam1.consensus_with_ambiguities({:region=>region, :case=>true})

  snps = cons_1.count_ambiguities
 
  snps_per_1k = (1000 * snps.to_f ) / region.size
  
  puts "#{r.id}\t#{region.size}\t#{snps}\t#{snps_per_1k}\n#{cons_1}"
  
end