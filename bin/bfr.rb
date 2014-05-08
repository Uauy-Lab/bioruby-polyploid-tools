#!/usr/bin/env ruby
require 'rubygems'
#require 'extensions/all'
require 'bio-samtools'
require 'optparse'

$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
$: << File.expand_path('.')
path=File.expand_path(File.dirname(__FILE__) + '/../lib/bioruby-polyploid-tools.rb')
$stderr.puts "Loading: #{path}"
require path

options = {}

options[:chunk] = 0
options[:chunk_size] = 0
options[:bucket] = 1

OptionParser.new do |opts|
  opts.banner = "Usage: bfr.rb [options]"
  
  opts.on("-r", "--reference FILE", "Fasta file with the reference sequence. Make sure to run faidx before running bfr in parallel") do |o|
     options[:reference] = o
   end
   
  opts.on("-a", "--parent_1 FILE", "Sorted BAM file with the alginments from parental 1") do |o|
    options[:parent_1] = o
  end
  
  opts.on("-b", "--parent_2 FILE", "Sorted BAM file with the alginments from parental 2") do |o|
    options[:parent_2] = o
  end
  
  opts.on("-c", "--bulk_1 FILE", "Sorted BAM file with the alginments from bulk1 1 (corresponding to the phenotype of parental 1)") do |o|
    options[:bulk_1] = o
  end
  
  opts.on("-d", "--bulk_2 FILE", "Sorted BAM file with the alginments from bulk1 2 (corresponding to the phenotype of parental 2)") do |o|
    options[:bulk_2] = o
  end
  
  opts.on("-o", "--bfr FILE", "Output file with the BFRs in the chunck") do |o|
    options[:output_filename] = o
  end

  opts.on("-s", "--stats FILE", "Output with the summary of the run. Only writes at the end, so in principle, paralell process should be able to write on it to get a status of how much has been completed.") do |o|
    options[:stats_file] = o
  end
  opts.on("-d", "--bulk_2 FILE", "Sorted BAM file with the alginments from bulk1 2 (corresponding to the phenotype of parental 2)") do |o|
    options[:bulk_2] = o
  end
  
  opts.on("-m", "--chunk_size FILE", "Sorted BAM file with the alginments from bulk1 2 (corresponding to the phenotype of parental 2)") do |o|
    options[:chunk_size] = o.to_i
  end
  
  opts.on("-n", "--chunk FILE", "Sorted BAM file with the alginments from bulk1 2 (corresponding to the phenotype of parental 2)") do |o|
    options[:chunk] = o.to_i
  end
  
    
end.parse!

p options
p ARGV


reference = options[:reference] 
chunk =  options[:chunk] 
chunk_size = options[:chunk_size]
output_filename =  options[:output_filename] 
stats_file = options[:stats_file]


min = chunk * chunk_size
max = min + chunk_size


parental_1=options[:parent_1]
parental_2=options[:parent_2]


bulk_1 = options[:bulk_1]
bulk_2 = options[:bulk_2]


fasta_db = Bio::DB::Fasta::FastaFile.new({:fasta=>reference})
fasta_db.load_fai_entries


if chunk_size == 0
  min = 0
  max = fasta_db.index.entries.size
end

container = Bio::BFRTools::BFRContainer.new

container.reference reference
container.parental_1  ( {:path => parental_1 } )
container.parental_2  ( {:path => parental_2 } )
container.bulk_1 ( {:path => bulk_1  })
container.bulk_2 ( {:path => bulk_2  })

i = -1

container.init_counters
output_file =  File.open(output_filename, "w")
puts "Range: #{min}:#{max}"
fasta_db.index.entries.each do | r |
  i = i  + 1
  #puts r
  #puts i
  next if i < min or i >= max 
  container.process_region({:region => r.get_full_region.to_s,:output_file => output_file } )
  #puts "Processed"
end
output_file.close

file_h = nil
if !File.exists? stats_file
  file_h = File.open(stats_file, "w")
  container.print_header({:output_file_stats => file_h})
else
  file_h =  File.open(stats_file, "a")
end
container.print_stats({:output_file_stats => file_h})

file_h.close
