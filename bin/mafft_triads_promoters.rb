#!/usr/bin/env ruby
require 'optparse'
require 'bio'
require 'csv'
require 'bio-blastxmlparser'
require 'fileutils'
require 'tmpdir'


options = {}
options[:identity] = 50
options[:min_bases] = 200
options[:split_token] = "-"
options[:tmp_folder]  = Dir.mktmpdir
options[:program]  = "blastn"
options[:random_sample] = 0

OptionParser.new do |opts|
  
  opts.banner = "Usage: filter_blat.rb [options]"

  opts.on("-i", "--identity FLOAT", "Minimum percentage identity") do |o|
    options[:identity] = o.to_f
  end
  opts.on("-c", "--min_bases int", "Minimum alignment length (default 200)") do |o|
    options[:min_bases] = o.to_i
  end

  opts.on("-t", "--triads FILE", "CSV file with the gene triad names in the named columns 'A','B' and 'D' ") do |o|
    options[:triads] = o
  end

  opts.on("-f", "--sequences FILE" , "FASTA file containing all the possible sequences. ") do |o|
    options[:fasta] = o
  end

  opts.on("-s", "--split_token CHAR", "Character used to split the sequence name. The name will be evarything before this token on the name of the sequences") do |o|
    options[:split_token] = o
  end

  opts.on("-p", "--program blastn|blastp", "The program to use in the alignments. Currntly only supported blastn and blastp") do |o|
    options[:program] = o
  end

  opts.on("-r", "--random_sample INT", "Number of blast to run and keep. If set, only the number of subsets will be run") do |o|
    options[:random_sample] = o.to_i
  end


end.parse!

def promoter_alignment(sequences_to_align) 
  options = ['--maxiterate', '1000', '--ep', '0', '--genafpair', '--quiet']
  options = ['--maxiterate', '1000', '--localpair', '--quiet']
  @mafft = Bio::MAFFT.new( "mafft" , options) unless @mafft 
  report = @mafft.query_align(sequences_to_align)
  report.alignment
end

def write_fasta_from_hash(sequences, filename)
  out = File.new(filename, "w")
  sequences.each_pair do | chromosome, exon_seq | 
    out.puts ">#{chromosome}\n#{exon_seq}\n"
  end
  out.close
end

split_token = options[:split_token]

sequences = Hash.new
sequence_count=0
Bio::FlatFile.open(Bio::FastaFormat, options[:fasta]) do |fasta_file|
  fasta_file.each do |entry|
    gene_name = entry.entry_id.split(split_token)[0]  
    sequences[gene_name] = entry unless sequences[gene_name]
    sequences[gene_name] = entry if entry.length > sequences[gene_name].length
    sequence_count += 1
  end
end

$stderr.puts "#Loaded #{sequences.length} genes from #{sequence_count} sequences"
#FileUtils.mkdir_p(options[:tmp_folder])
#$stderr.puts "TMP dir: #{options[:tmp_folder]}"

i =0 
CSV.foreach(options[:triads], headers:true ) do |row|
   a = row['A']
   b = row['B']
   d = row['D']
   triad = row['group_id']

 
   to_align = Bio::Alignment::SequenceHash.new 
   to_align[a] = sequences[a]
   to_align[b] = sequences[b]
   to_align[d] = sequences[d]

   prom_aln = promoter_alignment to_align

   cent_triad = triad.to_i / 100
   folder = "prom_aln/#{cent_triad}/"
   FileUtils.mkdir_p folder

   save_prom = "#{folder}/#{triad}.rom.fa"
   write_fasta_from_hash(prom_aln, save_prom)
   i += 1
   break if i > 10

end
