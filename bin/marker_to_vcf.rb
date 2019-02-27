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

options = {}
options[:min_identity] = 90
options[:filter_best]  = false

OptionParser.new do |opts|
  opts.banner = "Usage: marler_to_vcf.rb [options]"

  opts.on("-c", "--contigs FILE", "File with contigs to use as database") do |o|
    options[:path_to_contigs] = o
  end
  
  opts.on("-m", "--marker_list FILE", "File with the list of markers to search from") do |o|
    options[:marker_list] = o
  end
  
  opts.on("-b", "--filter_best", "If set, only keep the best alignment for each chromosome") do 
    options[:filter_best]  = false
  end

  opts.on("-i", "--min_identity INT", "Minimum identity to consider a hit (default 90)") do |o|
    options[:min_identity] = o.to_i
  end
  
  opts.on("-o", "--output FOLDER", "Output folder") do |o|
    options[:output_folder] = o
  end

  opts.on("-a", "--arm_selection #{Bio::PolyploidTools::ChromosomeArm.getValidFunctions.join('|')}", "Function to decide the chromome arm") do |o|
    options[:arm_selection] = Bio::PolyploidTools::ChromosomeArm.getArmSelection(o)
   end
  
  opts.on("-A", "--aligner exonerate|blast", "Select the aligner to use. Default: blast") do |o|
    raise "Invalid aligner" unless o == "exonerate" or o == "blast" 
    options[:aligner] = o.to_sym
  end

  opts.on("-d", "--database PREFIX", "Path to the blast database. Only used if the aligner is blast. The default is the name of the contigs file without extension.") do |o|
    options[:database] = o
  end
end.parse!
options[:database] = options[:path_to_contigs] 
p options
p ARGV


path_to_contigs=options[:path_to_contigs]

original_name="A"
snp_in="B"

fasta_reference = nil
test_file=options[:marker_list]

output_folder="#{test_file}_primer_design_#{Time.now.strftime('%Y%m%d-%H%M%S')}" 
output_folder= options[:output_folder] if  options[:output_folder]
Dir.mkdir(output_folder)
#T
temp_fasta_query="#{output_folder}/to_align.fa"
temp_contigs="#{output_folder}/contigs_tmp.fa"
exonerate_file="#{output_folder}/exonerate_tmp.tab"

min_identity= options[:min_identity]

@status_file="#{output_folder}/status.txt"


def write_status(status)
  f=File.open(@status_file, "a")
  f.puts "#{Time.now.to_s},#{status}"
  f.close
end


snps = Hash.new


write_status "Loading Reference"
#0. Load the fasta index 
fasta_reference_db = nil
if fasta_reference
  fasta_reference_db = Bio::DB::Fasta::FastaFile.new({:fasta=>fasta_reference})
  fasta_reference_db.load_fai_entries
  write_status "Fasta reference: #{fasta_reference}"
end

#1. Read all the SNP files 
#chromosome = nil
write_status "Reading SNPs"

File.open(test_file) do | f |
  f.each_line do | line |
    snp = Bio::PolyploidTools::SNPSequence.parse(line)  
    snp.genomes_count = options[:genomes_count]
    snp.snp_in = snp_in
    snp.original_name = original_name
    if snp.position 
      snps[snp.gene] = snp
    else
      $stderr.puts "ERROR: #{snp.gene} doesn't contain a SNP"
    end
  end
end

#2. Generate all the fasta files
write_status "Writing sequences to align"
written_seqs = Set.new
file = File.open(temp_fasta_query, "w")
snps.each_pair do |k,snp|
  unless written_seqs.include?(snp.gene)
    written_seqs << snp.gene 
    file.puts snp.to_fasta
  end
end
file.close


#3. Run exonerate on each of the possible chromosomes for the SNP
#puts chromosome
#chr_group = chromosome[0]
write_status "Searching markers in genome"
exo_f = File.open(exonerate_file, "w")
contigs_f = File.open(temp_contigs, "w") if options[:extract_found_contigs]
filename=path_to_contigs 
#puts filename
target=filename

fasta_file = Bio::DB::Fasta::FastaFile.new({:fasta=>target})
fasta_file.load_fai_entries
found_contigs = Set.new

def do_align(aln, exo_f, found_contigs, min_identity,fasta_file,options)
  if aln.identity > min_identity
    exo_f.puts aln.line
  end  
end

Bio::DB::Blast.align({:query=>temp_fasta_query, :target=>options[:database]}) do |aln|
  do_align(aln, exo_f, found_contigs,min_identity, fasta_file,options)
end

exo_f.close() 

def print_positions(min_identity:90, filter_best:false, exonerate_filename:"test.exo", snps:{}) 
  File.open(exonerate_filename) do |f|
    f.each_line do | line |
      record = Bio::DB::Exonerate::Alignment.parse_custom(line)
      next unless  record and record.identity >= min_identity
      snp = snps[record.query_id]                           
      next unless snp != nil and snp.position.between?( (record.query_start + 1) , record.query_end)
      begin
        puts record.query_id
        position = record.query_position_on_target(snp.position)
        q_strand = record.query_strand
        t_strand = record.target_strand
        template = snp.template_sequence
        puts template
        puts position
        puts "\n"
        #puts exon.inspect
      rescue Bio::DB::Exonerate::ExonerateException
        $stderr.puts "Failed for the range #{record.query_start}-#{record.query_end} for position #{snp.position}"     
      end
    end
  end
end

puts "SNPs:"
#puts snps.inspect
print_positions(exonerate_filename:exonerate_file, min_identity:98, snps:snps)



