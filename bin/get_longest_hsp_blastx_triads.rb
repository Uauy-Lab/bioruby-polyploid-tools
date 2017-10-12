#!/usr/bin/env ruby
require 'optparse'
require 'bio'
require 'csv'
#$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
#$: << File.expand_path('.')
#path= File.expand_path(File.dirname(__FILE__) + '/../lib/bio-polyploid-tools.rb')
#require path

options = {}
options[:identity] = 50
options[:min_bases] = 200
options[:blastx] = "-"

OptionParser.new do |opts|
  
  opts.banner = "Usage: filter_blat.rb [options]"

  opts.on("-p", "--blastx FILE", "BLAST XML  file") do |o|
    options[:blastx] = o
  end
  opts.on("-i", "--identity FLOAT", "Minimum percentage identity") do |o|
    options[:identity] = o.to_f
  end
  opts.on("-c", "--min_bases int", "Minimum alignment length (default 200)") do |o|
    options[:min_bases] = o.to_i
  end

  opts.on("-t", "--triads FILE", "CSV file with the gene triad names in the named columns 'A','B' and 'D' ") do |o|
    options[:triads] = o
  end
  
end.parse!

valid_pairs_A_B = Hash.new
valid_pairs_A_D = Hash.new
valid_pairs_B_D = Hash.new

CSV.foreach(options[:triads], headers:true ) do |row|
  valid_pairs_A_B[row['A']] = row['B']
  valid_pairs_A_D[row['A']] = row['D']
  valid_pairs_B_D[row['B']] = row['D']
end

stream = ARGF
stream = IO.open(options[:blastx]) unless  options[:blastx] == "-"
puts "Loaded #{valid_pairs_B_D.length} triads"
$stdout.flush

blast_report = Bio::FlatFile.new(Bio::Blast::Report, stream)

blast_report.each_entry do |report| 
  puts "Hits for " + report.query_def + " against " + report.db
  $stdout.flush
  report.each do |hit|
    query  = hit.query_id.split("-")[0]
    target = hit.target_id.split("-")[0]
    if valid_pairs_A_B[query] == target or valid_pairs_A_D[query] == target or valid_pairs_B_D[query] == target  
      puts hit.target_id, "\t", hit.evalue, "\n" if hit.evalue < 0.001
      puts hit.inspect
    end
    
  end
end

stream.close unless  options[:blat_file] == "-"
