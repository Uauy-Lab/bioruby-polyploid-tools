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



fasta_db = Bio::DB::Fasta::FastaFile.new( ARGV[0])
fasta_db.load_fai_entries
bam1 =  Bio::DB::Sam.new({:fasta=>ARGV[0], :bam=>ARGV[1]})
bam2 =  Bio::DB::Sam.new({:fasta=>ARGV[0], :bam=>ARGV[2]})


output_prefix = ARGV[3]

block_size=1000

min_cov = ARGV[4].to_i ? ARGV[4].to_i : 10
chunk = ARGV[5].to_i
chunk_size = ARGV[6].to_i




main_table="#{output_prefix}_#{block_size}_#{min_cov}_table.#{chunk}.csv"

table_file = File.open(main_table, "w")
table_file.puts "gene\tlength\tsnps_1\tcalled_1\tsnps_per_#{block_size}_1\tsnps_2\tcalled_2\tsnps_per_#{block_size}_2\tsnps_tot\tsnps_per_1k_tot"

hist_1= Hash.new(0)
hist_2= Hash.new(0)

fasta_file = File.open("#{output_prefix}_#{min_cov}.#{chunk}.fa", "w")
i = -1
min = chunk * chunk_size
max = min + chunk_size

fasta_db.index.entries.each do | r |
  i = i  + 1
  next if i < min or i >= max
  #Np r.get_full_region
  #container.process_region( { :region => r.get_full_region.to_s, :output_file => output_file } )
  region=r.get_full_region


  begin
    reg_a = bam1.fetch_region({:region=>region,  :min_cov=>min_cov, :A=>1})
    reg_b = bam2.fetch_region({:region=>region,  :min_cov=>min_cov, :A=>1})
    cons_1 = reg_a.consensus
    cons_2 = reg_b.consensus
   

    snps_1 = cons_1.count_ambiguities
    snps_2 = cons_2.count_ambiguities
    
    called_1 = reg_a.called
    called_2 = reg_b.called

    snps_tot = Bio::Sequence.snps_between(cons_1, cons_2)

    snps_per_1k_1   = (block_size * snps_1.to_f   ) / region.size
    snps_per_1k_2   = (block_size * snps_2.to_f   ) / region.size
    snps_per_1k_tot = (block_size * snps_tot.to_f ) / region.size

    hist_1[snps_per_1k_1.to_i] += 1
    hist_2[snps_per_1k_2.to_i] += 1

    table_file.print "#{r.id}\t#{region.size}\t"
    table_file.print "#{snps_1}\t#{called_1}\t#{snps_per_1k_1}\t"
    table_file.print "#{snps_2}\t#{called_2}\t#{snps_per_1k_2}\t"
    table_file.print "#{snps_tot}\t#{snps_per_1k_tot}\n"
    fasta_file.puts ">#{r.id}_1"
    fasta_file.puts "#{cons_1}"
    fasta_file.puts ">#{r.id}_2"
    fasta_file.puts "#{cons_2}"
    
  rescue Exception => e
    $stderr.puts "Unable to process #{region}: #{e.to_s}"
  end
end
fasta_file.close
table_file.close

hist_table="#{output_prefix}_#{block_size}_#{min_cov}_hist.#{chunk}.csv"
hist_file = File.open(hist_table, "w")

all_keys = SortedSet.new(hist_1.keys)
all_keys.merge(hist_2.keys)
hist_file.puts "SNPs/#{block_size}\thist_1\thist_2\n"
all_keys.each do |k|
  hist_file.puts "#{k}\t#{hist_1[k]}\t#{hist_2[k]}"
end

hist_file.close



