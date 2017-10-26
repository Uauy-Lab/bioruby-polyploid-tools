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
options[:output_folder]  = "."
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

  opts.on("-o", "--output_folder DIR", "Folder to save the output") do |o|
    options[:output_folder] = o
  end


end.parse!

module Bio::Alignment::EnumerableExtension
  def each_base_alignment
    names = self.keys 

    i = 0
    len = 0 
    len = self[names[0]].length if names[0]
    total_alignments = names.size  
    while i < len  do
      yield names.map { | chr| self[chr][i]  }
      i += 1
    end
  end

  def cut_alignment(start, length)
    a = Bio::Alignment::SequenceHash.new
    a.set_all_property(get_all_property)
    each_pair do |key, str|
      seq = ""
      seq = str[start, length] if str != nil
      a.store(key, seq)
    end
    a
  end

  def best_block
    best_start = 0
    best_score = 0
    best_end = 0
    best_length = 0
    current_start = 0
    current_score = 0
    current_length = 0

    each_base_alignment_with_index do |bases, i|
      current_start = i if current_length == 0
      current_length += 1
      current_score += sum_of_pair bases
      if current_score > best_score
        best_score = current_score
        best_length = current_length
        best_end = i 
        best_start = current_start
      end

      if current_score < 0
        current_length = 0
        current_score = 0
      end

    end

    [best_start, best_length, len - best_start - best_length , len - best_start ]
  end

  def each_base_alignment_with_index
    names = self.keys 
    total_alignments = names.size  
    i = 0
    while i < len  do
      yield names.map { | chr| self[chr][i] } , i
      i += 1
    end
  end

  def each_base_alignment
    each_base_alignment_with_index do |chr, i|
      yield chr
    end
  end

  def sum_of_all_pairs
    return @sum_of_all_pairs if @sum_of_all_pairs
    @sum_of_all_pairs = 0
    self.each_base_alignment do |bases|
      @sum_of_all_pairs += sum_of_pair bases
    end
    @sum_of_all_pairs
  end

  def sum_of_identities
    return @sum_of_identities if @sum_of_identities
    @sum_of_identities = 0
    self.each_base_alignment do |bases|
      @sum_of_identities += s_o_i bases
    end
    @sum_of_identities
  end

  def len
    return @len if @len
    names = self.keys 
    @len = 0 
    @len = self[names[0]].length if names[0] and self[names[0]] != nil
    @len
  end

  def pairwise_comparaisons
    names = self.keys 
    n = names.size 
    c = n * (n-1)/2
    c
  end

  def identity
    max_score = len * pairwise_comparaisons
    sum_of_identities.to_f/max_score
  end

  def normalized_sum_of_all_pairs
    max_score = len * pairwise_comparaisons
    sum_of_all_pairs.to_f/max_score
  end

  def sum_of_pair(bases)
    x = bases.length - 1
    total  = 0
    for i in 0..x 
      y = i + 1
      for j in y..x
        case 
        when (bases[i] == "-" and bases[j] == "-")
          total += 0
        when (bases[i] == "N" and bases[j] == "N")
          total += 0
        when (bases[i] == "n" and bases[j] == "n")
          total += 0
        when (bases[i] == "-" or bases[j] == "-")
          total -= 2
        when bases[i] ==  bases[j]
          total += 1
        when  bases[i] !=  bases[j]
          total -= 1
        else
          $stderr.puts "Invalid comparaison! sum_of_all_pairs(#{bases})"
        end
      end
    end
    total 
  end

  def s_o_i(bases)
    x = bases.length - 1
    total  = 0
    for i in 0..x 
      y = i + 1
      for j in y..x
        total += 1 if bases[i] ==  bases[j]
      end
    end
    total 
  end

  def window_identities(window_size=100, offset=25)
    steps = (0..len).step(offset).to_a.map {|a| a + len%offset }.reverse
    ret = []
    steps.each_with_index do |e, i|
      start   = e - window_size
      tmp_aln = self.cut_alignment start, window_size
      tmp_arr = [
        i * offset, 
        i * offset + window_size,
        tmp_aln.sum_of_all_pairs, 
        tmp_aln.normalized_sum_of_all_pairs, 
        tmp_aln.sum_of_identities, 
        tmp_aln.identity]
      ret << tmp_arr
    end
    ret
  end
end

def promoter_alignment(sequences_to_align) 
  process = true
  sequences_to_align.each_value { |val| process &= val != nil }
  return sequences_to_align unless process
 #options = ['--maxiterate', '1000', '--ep', '0', '--genafpair', '--quiet']
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

def get_longest_aln(aln, max_gap: 10)
  names = aln.keys   
  i = 0
  len = 0  
  len = aln[names[0]].length if names[0] and aln[names[0]] != nil
  total_alignments = names.size
  masked_snps = "-" * len  
  longest_start = -1
  longest_length = 0
  current_start = -1
  current_length = 0
  current_gap = 0 
  longest_gaps = 0
  gaps = 0
  while i < len  do
    different = 0
    cov = 0
    names.each do | chr |
      if aln[chr][i]  != "-"
        cov += 1 
      end
    end
    if cov == total_alignments
      current_start = i if current_length == 0
      current_length += 1 
      current_gap = 0
    else
      gaps += 1
      current_gap += 1 
    end

    if current_length > longest_length
      longest_length = current_length
      longest_start  = current_start
      longest_gaps = gaps - current_gap
    end
    if current_gap > max_gap
      current_length = 0
      gaps = 0
    end
    i += 1
  end
  longest_length += longest_gaps
  [longest_start, longest_length, len - longest_start - longest_length, len - longest_start]
end

split_token = options[:split_token]

def read_alignments(fasta_path, split_token)
  sequences = Hash.new
  sequence_count=0
  Bio::FlatFile.open(Bio::FastaFormat, fasta_path) do |fasta_file|
    fasta_file.each do |entry|
      #puts entry
      gene_name = entry.entry_id.split(split_token)[0]  
      sequences[gene_name] = entry unless sequences[gene_name]
      sequences[gene_name] = entry if entry.length > sequences[gene_name].length
      sequence_count += 1
    end
  end
  [sequences,sequence_count]
end

sequences, sequence_count = read_alignments(options[:fasta], split_token)

$stderr.puts "#Loaded #{sequences.length} genes from #{sequence_count} sequences"
output_folder = options[:output_folder]

FileUtils.mkdir_p output_folder
summary_file    = "#{output_folder}/identities.txt" 
long_table_file = "#{output_folder}/sliding_window_identities.txt"

out = File.open(summary_file, "w")
long_table = File.open(long_table_file, "w")

i =0 

header =  ["triad", "total_aln_length"]
header << ["longest_start", "longest_length", "longest_start_from_CDS","longest_end_from_CDS", "longest_sum_of_all_pairs","longest_norm_sum_of_all_pairs","longest_sum_of_identities", "longest_identity"]
header << ["best_start",    "best_length"  ,  "best_start_from_CDS","best_end_from_CDS", "best_sum_of_all_pairs","best_norm_sum_of_all_pairs","best_sum_of_identities", "best_identity"]
out.puts header.join("\t")
long_table.puts ["triad", "type", "start_from_CDS", "end_from_cds" , "sum_of_all_pairs","norm_sum_of_all_pairs","sum_of_identities", "identity"].join("\t")
CSV.foreach( options[:triads], headers:true ) do |row|
 a = row['A']
 b = row['B']
 d = row['D']
 triad = row['group_id']

 cent_triad = triad.to_i / 100
 folder = "#{output_folder}/prom_aln/#{cent_triad}/"
 save_prom = "#{folder}/#{triad}.prom.fa"
 
 to_align = Bio::Alignment::SequenceHash.new 
 to_align[a] = sequences[a]
 to_align[b] = sequences[b]
 to_align[d] = sequences[d]

 prom_aln = nil
 unless File.file? save_prom
  prom_aln = promoter_alignment to_align 
 else
  ff, seqs_cnt = read_alignments save_prom, split_token
  seqs = Bio::Alignment::SequenceHash.new 
  prom_aln = Bio::Alignment.new(ff)
 end
 print_arr = [triad, prom_aln.len]
 aln_stats = get_longest_aln prom_aln
 print_arr << aln_stats
 cut_seqs = prom_aln.cut_alignment aln_stats[0], aln_stats[1]



 print_arr << cut_seqs.sum_of_all_pairs
 print_arr << cut_seqs.normalized_sum_of_all_pairs

 print_arr << cut_seqs.sum_of_identities
 print_arr << cut_seqs.identity

 best_aln_stats = prom_aln.best_block
 best_aln_cut = prom_aln.cut_alignment best_aln_stats[0], best_aln_stats[1]

 print_arr << best_aln_stats

 print_arr << best_aln_cut.sum_of_all_pairs
 print_arr << best_aln_cut.normalized_sum_of_all_pairs

 print_arr << best_aln_cut.sum_of_identities
 print_arr << best_aln_cut.identity

 base = [triad, "cut_longest_region"]
 cut_seqs.window_identities.each do |e|    
   long_table.puts [base, e].flatten.join("\t")
 end

  base = [triad, "cut_best_region"]
 best_aln_cut.window_identities.each do |e|    
   long_table.puts [base, e].flatten.join("\t")
 end

 base = [triad, "full_promoter"]
 prom_aln.window_identities.each do |e|    
  long_table.puts [base, e].flatten.join("\t")
 end

 out.puts print_arr.join("\t")

 FileUtils.mkdir_p folder

 write_fasta_from_hash(prom_aln, save_prom) unless File.file?(save_prom)

 save_prom_cut = "#{folder}/#{triad}.prom.cut.fa"
 write_fasta_from_hash(cut_seqs, save_prom_cut)  unless File.file?(save_prom)

 save_prom_cut_best = "#{folder}/#{triad}.prom.cut.best.fa"
 write_fasta_from_hash(best_aln_cut, save_prom_cut_best)  

 i += 1
 #break if i > 10
end
long_table.close
out.close
