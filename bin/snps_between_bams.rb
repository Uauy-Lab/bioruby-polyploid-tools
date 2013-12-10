#!/usr/bin/env ruby

require 'bio'
require 'rubygems'
require 'pathname'
require 'bio-samtools'

require 'set'

$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
$: << File.expand_path('.')
path=File.expand_path(File.dirname(__FILE__) + '/../lib/bioruby-polyploid-tools.rb')
$stderr.puts "Loading: #{path}"
require path

puts  ARGV[0]

fasta_db = Bio::DB::Fasta::FastaFile.new( ARGV[0])
fasta_db.load_fai_entries
bam1 =  Bio::DB::Sam.new({:fasta=>ARGV[0], :bam=>ARGV[1]})
bam2 =  Bio::DB::Sam.new({:fasta=>ARGV[0], :bam=>ARGV[2]})

output_prefix = ARGV[3]

block_size=300

main_table="#{output_prefix}_#{block_size}_table.csv"

table_file = File.open(main_table, "w")
table_file.puts "gene\tlength\tsnps_1\tcalled_1\tsnps_per_#{block_size}_1\tsnps_2\tcalled_2\tsnps_per_#{block_size}_2\tsnps_tot\tsnps_per_1k_tot"

hist_1= Hash.new(0)
hist_2= Hash.new(0)

fasta_file = File.open("#{output_prefix}.fa", "w")
fasta_db.index.entries.each do | r |
  #Np r.get_full_region
  #container.process_region( { :region => r.get_full_region.to_s, :output_file => output_file } )
  region=r.get_full_region


  begin

    cons_1 = bam1.consensus_with_ambiguities({:region=>region, :case=>true})
    cons_2 = bam2.consensus_with_ambiguities({:region=>region, :case=>true})
    if cons_1 != cons_2

      snps_1 = cons_1.count_ambiguities
      snps_2 = cons_2.count_ambiguities

      called_1 = cons_1.upper_case_count
      called_2 = cons_2.upper_case_count

      snps_tot = Bio::Sequence.snps_between(cons_1, cons_2)

      snps_per_1k_1   = (block_size * snps_1.to_f   ) / called_1
      snps_per_1k_2   = (block_size * snps_2.to_f   ) / called_2
      snps_per_1k_tot = (block_size * snps_tot.to_f ) / region.size

      hist_1[snps_per_1k_1.to_i] += 1
      hist_2[snps_per_1k_2.to_i] += 1

      table_file.puts "#{r.id}\t#{region.size}\t#{snps_1}\t#{called_1}#{snps_per_1k_1}\t#{snps_2}\t#{called_2}\t#{snps_per_1k_2}\t#{snps_tot}\t#{snps_per_1k_tot}"
      fasta_file.puts ">#{r.id}_1"
      fasta_file.puts "#{cons_1}"
      fasta_file.puts ">#{r.id}_2"
      fasta_file.puts "#{cons_2}"
    end
  rescue Exception => e
    $stderr.puts "Unable to process #{region}: #{e.to_s}"
  end
end
fasta_file.close
table_file.close

hist_table="#{output_prefix}_#{block_size}_hist.csv"
hist_file = File.open(hist_table, "w")

all_keys = SortedSet.new(hist_1.keys)
all_keys.merge(hist_2.keys)
hist_file.puts "SNPs/#{block_size}\thist_1\thist_2\n"
all_keys.each do |k|
  hist_file.puts "#{k}\t#{hist_1[k]}\t#{hist_2[k]}"
end

hist_file.close



