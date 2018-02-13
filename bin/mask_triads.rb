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
  
  opts.banner = "Usage: mask_triads.rb [options]"

  opts.on("-t", "--triads FILE", "CSV file with the gene triad names in the named columns 'A','B' and 'D' ") do |o|
    options[:triads] = o
  end

  opts.on("-f", "--pep FILE" , "FASTA file containing all the possible peptide sequences. ") do |o|
    options[:pep] = o
  end

  opts.on("-s", "--cds FILE" , "FASTA file containing all the possible CDS sequences. ") do |o|
    options[:cds] = o
  end

  opts.on("-s", "--split_token CHAR", "Character used to split the sequence name. The name will be evarything before this token on the name of the sequences") do |o|
    options[:split_token] = o
  end

end.parse!


def peptide_alignment(sequences_to_align) 
  options = ['--maxiterate', '1000', '--localpair', '--quiet']
  mafft = Bio::MAFFT.new( "mafft" , options)
  report = mafft.query_align(sequences_to_align)
  report.alignment
end


split_token = options[:split_token]

pep_seq = Hash.new
pep_seq_count=0
Bio::FlatFile.open(Bio::FastaFormat, options[:pep]) do |fasta_file|
  fasta_file.each do |entry|
    gene_name = entry.entry_id.split(split_token)[0]  
    pep_seq[gene_name] = entry unless pep_seq[gene_name]
    pep_seq[gene_name] = entry if entry.length > pep_seq[gene_name].length
    pep_seq_count += 1
  end
end
$stderr.puts "#Loaded #{pep_seq.length} genes from #{pep_seq_count} pep_seq"

cds_seq = Hash.new
cds_seq_count=0
Bio::FlatFile.open(Bio::FastaFormat, options[:cds]) do |fasta_file|
  fasta_file.each do |entry|
    gene_name = entry.entry_id.split(split_token)[0]  
    cds_seq[gene_name] = entry unless cds_seq[gene_name]
    cds_seq[gene_name] = entry if entry.length > cds_seq[gene_name].length
    cds_seq_count += 1
  end
end
$stderr.puts "#Loaded #{cds_seq.length} genes from #{cds_seq_count} cds_seq"


$stderr.puts "TMP dir: #{options[:tmp_folder]}"

def write_fasta_from_hash(sequences, filename)
  out = File.new(filename, "w")
  #puts sequences.inspect
  sequences.each_pair do | chromosome, exon_seq | 
    out.puts ">#{chromosome}\n#{exon_seq}\n"
  end
  out.close
end


CSV.foreach(options[:triads], headers:true ) do |row|
   a = row['A']
   b = row['B']
   d = row['D']
   triad = row['group_id']

   to_align = Bio::Alignment::SequenceHash.new 
   to_align[a] = pep_seq[a]
   to_align[b] = pep_seq[b]
   to_align[d] = pep_seq[d]

   cds_seqs = Bio::Alignment::SequenceHash.new 
   cds_seqs[a] = cds_seq[a].to_biosequence
   cds_seqs[b] = cds_seq[b].to_biosequence
   cds_seqs[d] = cds_seq[d].to_biosequence

   cent_triad = triad.to_i / 100
   folder = "alignments/#{cent_triad}/"
   FileUtils.mkdir_p folder
   
   pep_aln = peptide_alignment(to_align)
   
   save_pep = "#{folder}/#{triad}.pep.fa"
   write_fasta_from_hash(pep_aln, save_pep)

   save_cds = "#{folder}/#{triad}.cds.fa"
   write_fasta_from_hash(cds_seqs, save_cds)
  #break
end
