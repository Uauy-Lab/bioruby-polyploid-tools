#!/usr/bin/env ruby
require 'optparse'

require 'csv'
require 'fileutils'
require 'tmpdir'
require 'bio-samtools'
require 'bio'

opts = {}
opts[:identity] = 50
opts[:min_bases] = 200
opts[:split_token] = "."
opts[:tmp_folder]  = Dir.mktmpdir
opts[:program]  = "blastn"
opts[:random_sample] = 0

OptionParser.new do |o|
  
  o.banner = "Usage: mask_triads.rb [options]"

  o.on("-t", "--triads FILE", "CSV file with the gene triad names in the named columns 'A','B' and 'D' ") do |o|
    opts[:triads] = o
  end

  o.on("-f", "--fasta FILE" , "FASTA file containing all the possible peptide sequences. ") do |o|
    opts[:fasta] = o
  end

  o.on("-s", "--split_token CHAR", "Character used to split the sequence name. The name will be evarything before this token on the name of the sequences") do |o|
    opts[:split_token] = o
  end

end.parse!


split_token = opts[:split_token]

fasta_reference_db = Bio::DB::Fasta::FastaFile.new(fasta: opts[:fasta])
fasta_reference_db.load_fai_entries
#puts fasta_reference_db.index.entries
cannonical = Hash.new 
fasta_reference_db.index.entries.each do |e|  
  gene = e.id.split(split_token)[0]
  cannonical[gene] = e unless cannonical[gene]
  cannonical[gene]  = e if   e.length > cannonical[gene].length
end

$stderr.puts "#Loaded #{cannonical.length} canonical sequences from #{fasta_reference_db.index.size} in reference"

$stderr.puts "TMP dir: #{opts[:tmp_folder]}"

def write_fasta_from_hash(sequences, filename)
  out = File.new(filename, "w")
  sequences.each_pair do | chromosome, exon_seq | 
    out.puts ">#{chromosome}\n#{exon_seq}\n"
  end
  out.close
end




CSV.foreach(opts[:triads], headers:true ) do |row|
  next unless row["cardinality_abs"] == "1:1:1"
  mafft_opts = ['--maxiterate', '1000', '--localpair', '--quiet']
  mafft = Bio::MAFFT.new( "mafft" , mafft_opts)
   a = row['A']
   b = row['B']
   d = row['D']
   triad = row['group_id']

   to_align = Bio::Alignment::SequenceHash.new 
   puts cannonical[a].get_full_region
   seq_a = fasta_reference_db.fetch_sequence(cannonical[a].get_full_region)
   seq_b = fasta_reference_db.fetch_sequence(cannonical[b].get_full_region)
   seq_d = fasta_reference_db.fetch_sequence(cannonical[d].get_full_region)
   to_align[a] = seq_a
   to_align[b] = seq_b
   to_align[d] = seq_d
  
   report = mafft.query_alignment(to_align)

   #cent_triad = triad.to_i / 100
   #folder = "alignments/#{cent_triad}/"
   #FileUtils.mkdir_p folder
   
   #report = mafft.query_align(to_align)
   #aln = report.alignment

   #save_cds = "#{folder}/#{triad}.cds.fa"
   #write_fasta_from_hash(cds_seqs, save_cds)
  #break
end
