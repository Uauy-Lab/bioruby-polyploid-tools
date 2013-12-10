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
bam2 =  Bio::DB::Sam.new({:fasta=>ARGV[0], :bam=>ARGV[2]})

fasta_db.index.entries.each do | r |
  #Np r.get_full_region
  #container.process_region( { :region => r.get_full_region.to_s, :output_file => output_file } )
  region=r.get_full_region
  
  

  cons_1 = bam1.consensus_with_ambiguities({:region=>region, :case=>true})
  cons_2 = bam2.consensus_with_ambiguities({:region=>region, :case=>true})
  if cons_1 != cons_2
    
    snps_1 = cons_1.count_ambiguities
    snps_2 = cons_2.count_ambiguities
    
    snps_tot = Bio::Sequence.snps_between(cons_1, cons_2)
  
    snps_per_1k_1   = (1000 * snps_1.to_f   ) / region.size
    snps_per_1k_2   = (1000 * snps_2.to_f   ) / region.size
    snps_per_1k_tot = (1000 * snps_tot.to_f ) / region.size
  
    puts "#{r.id}\t#{region.size}\t#{snps_1}\t#{snps_per_1k_1}\t#{snps_2}\t#{snps_per_1k_2}\t#{snps_tot}\t#{snps_per_1k_tot}"
  end
end