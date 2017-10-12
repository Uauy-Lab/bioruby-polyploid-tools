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

OptionParser.new do |opts|
  
  opts.banner = "Usage: filter_blat.rb [options]"

#  opts.on("-p", "--psl FILE", "PSL file") do |o|
#    options[:blat_file] = o
#  end
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


end.parse!


def blast_pair_fast(path_a, path_b, out_path, program: "blastn")
  cmd = "#{program} -query #{path_a} -subject #{path_b} -task #{program} -out #{out_path} -outfmt '5' "
  #puts cmd
  executed = system cmd
  result = []
  blast_version = nil
  n = Bio::BlastXMLParser::XmlIterator.new(out_path).to_enum
  longest = nil
  max_length = 0
  max_pident = 0.0
  n.each do | iter |
    iter.each do | hit |
      hit.each do | hsp |
        if hsp.align_len > max_length
          max_length = hsp.align_len
          max_pident = 100 * hsp.identity.to_f / hsp.align_len.to_f
        end
      end
    end
  end
  [max_length, max_pident]
end

valid_pairs_A_B = Hash.new
valid_pairs_A_D = Hash.new
valid_pairs_B_D = Hash.new

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
$stderr.puts "TMP dir: #{options[:tmp_folder]}"

a_tmp   = options[:tmp_folder] + "/A.fa"
b_tmp   = options[:tmp_folder] + "/B.fa"
d_tmp   = options[:tmp_folder] + "/D.fa"
out_tmp = options[:tmp_folder] + "/out.blast"

puts [
  "group_id" , "query"      , "subject" , 
  "chr_query", "chr_subject", "aln_type",
  "length"   , "pident"    ].join("\t")

CSV.foreach(options[:triads], headers:true ) do |row|
   a = row['A']
   b = row['B']
   d = row['D']
   triad = row['group_id']
   seq_a = sequences[a]
   seq_b = sequences[b]
   seq_d = sequences[d]
   File.open(a_tmp, 'w') {|f| f.write(seq_a) } if seq_a
   File.open(b_tmp, 'w') {|f| f.write(seq_b) } if seq_b
   File.open(d_tmp, 'w') {|f| f.write(seq_d) } if seq_d
  if seq_a and seq_b
      to_print = [triad, a, b , "A","B","A->B"]
      to_print << blast_pair_fast(a_tmp, b_tmp, out_tmp, program:options[:program])
      puts to_print.join("\t")
   end
  if seq_a and seq_d
      to_print = [triad, a, b , "A","D","A->D"]
      to_print << blast_pair_fast(a_tmp, d_tmp, out_tmp, program:options[:program])
      puts to_print.join("\t")
  end
  if seq_b and seq_d
      to_print = [triad, a, b , "B","D","B->D"]
      to_print << blast_pair_fast(b_tmp, d_tmp, out_tmp, program:options[:program])
      puts to_print.join("\t")
  end
end
