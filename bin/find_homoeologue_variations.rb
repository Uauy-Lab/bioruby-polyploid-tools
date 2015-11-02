#!/usr/bin/env ruby
require 'bio'
require 'rubygems'
require 'pathname'
require 'bio-samtools'
require 'optparse'
require 'set'
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
options[:model] = "est2genome"
options[:min_identity] = 90
options[:extract_found_contigs] = false
OptionParser.new do |opts|
  
  opts.banner = "Usage: find_homoeologue_variations.rb [options]"

  opts.on("-c", "--sequences FASTA", "Sequence of the region to searc") do |o|
    options[:sequences] = o
  end
  opts.on("-r", "--reference FASTA", "reference with the contigs") do |o|
    options[:reference] = o
  end
  opts.on("-o", "--output DIR", "Directory to write the output") do |o|
    options[:output] = o
  end

  opts.on("-x", "--extract_found_contigs", "If present, save in a separate file the contigs with matches. Useful to debug.") do |o|
    options[:extract_found_contigs] = true
  end
  
end.parse!
#reference="/Users/ramirezr/Documents/TGAC/references/Triticum_aestivum.IWGSP1.21.dna_rm.genome.fa"
reference = options[:reference] if options[:reference]
throw raise Exception.new(), "Reference has to be provided" unless reference
sequences = options[:sequences] if options[:sequences]
throw raise Exception.new(), "Fasta file with sequences has to be provided" unless sequences
output_folder = options[:output] if options[:output]
throw raise Exception.new(), "An output directory has to be provided" unless output_folder
model=options[:model] 
Dir.mkdir(output_folder)
min_identity= options[:min_identity]

exonerate_file="#{output_folder}/exonerate_tmp.tab"
temp_contigs="#{output_folder}/contigs_tmp.fa"

fasta_file = Bio::DB::Fasta::FastaFile.new({:fasta=>reference})
fasta_file.load_fai_entries

log "Searching markers in genome"
found_contigs = Set.new
exo_f = File.open(exonerate_file, "w")
contigs_f = File.open(temp_contigs, "w") if options[:extract_found_contigs]
Bio::DB::Exonerate.align({:query=>sequences, :target=>reference, :model=>model}) do |aln|
	if aln.identity > min_identity
    exo_f.puts aln.line
    unless found_contigs.include?(aln.target_id) #We only add once each contig. Should reduce the size of the output file. 
      found_contigs.add(aln.target_id)
      entry = fasta_file.index.region_for_entry(aln.target_id)
      raise ExonerateException.new,  "Entry not found! #{aln.target_id}. Make sure that the #{target_id}.fai was generated properly." if entry == nil
      region = entry.get_full_region
      seq = fasta_file.fetch_sequence(region)
      contigs_f.puts(">#{aln.target_id}\n#{seq}") if options[:extract_found_contigs]
    end
  end  
end
exo_f.close() 
contigs_f.close() if options[:extract_found_contigs]


log "Reading best alignment on each chromosome"

container= Bio::PolyploidTools::ExonContainer.new
container.flanking_size=options[:flanking_size] 
container.gene_models(temp_fasta_query)
container.chromosomes(target)
container.add_parental({:name=>"A"})
container.add_parental({:name=>"B"})
snps.each do |snp|
  snp.container = container
  snp.flanking_size = container.flanking_size
  snp.variation_free_region = options[:variation_free_region]
  container.add_snp(snp)
end
container.add_alignments({:exonerate_file=>exonerate_file, :arm_selection=>options[:arm_selection] , :min_identity=>min_identity})


