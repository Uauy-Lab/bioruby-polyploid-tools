#!/usr/bin/env ruby
require 'optparse'

require 'csv'
require 'fileutils'
require 'tmpdir'
require 'bio-samtools'
require 'bio'

$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
$: << File.expand_path('.')
path= File.expand_path(File.dirname(__FILE__) + '/../lib/bioruby-polyploid-tools.rb')
require path
opts = {}
opts[:identity] = 50
opts[:min_bases] = 200
opts[:split_token] = "."
opts[:tmp_folder]  = Dir.mktmpdir
opts[:random_sample] = 0
opts[:output_folder] = "."

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

  o.on("-o", "--output_folder DIR", "Location to save the alignment masks. If the alignment exists, it is recycled to avoid calling MAFFT again") do |o|
    opts[:output_folder] = o
  end
end.parse!


split_token = opts[:split_token]
reference_name = File.basename opts[:fasta]
output_folder = opts[:output_folder]
@fasta_reference_db = Bio::DB::Fasta::FastaFile.new(fasta: opts[:fasta])
@fasta_reference_db.load_fai_entries
#puts @fasta_reference_db.index.entries
@cannonical = Hash.new
@fasta_reference_db.index.entries.each do |e|
  gene = e.id.split(split_token)[0]
  @cannonical[gene] = e unless @cannonical[gene]
  @cannonical[gene]  = e if   e.length > @cannonical[gene].length
end

$stderr.puts "#Loaded #{@cannonical.length} canonical sequences from #{@fasta_reference_db.index.size} in reference"

$stderr.puts "TMP dir: #{opts[:tmp_folder]}"

def write_fasta_from_hash(sequences, filename)
  out = File.new(filename, "w")
  sequences.each_pair do | chromosome, exon_seq |
    out.puts ">#{chromosome}\n#{exon_seq}\n"
  end
  out.close
end

def mafft_align(a, b, d)
  to_align = Bio::Alignment::SequenceHash.new
  seq_a = @fasta_reference_db.fetch_sequence(@cannonical[a].get_full_region)
  seq_b = @fasta_reference_db.fetch_sequence(@cannonical[b].get_full_region)
  seq_d = @fasta_reference_db.fetch_sequence(@cannonical[d].get_full_region)
  to_align[a] = seq_a
  to_align[b] = seq_b
  to_align[d] = seq_d
  report = mafft.query_alignment(to_align)
  aln = report.alignment
  aln
end

def read_alignment(path)
  aln = Bio::Alignment::SequenceHash.new
  i = 0
  Bio::FlatFile.open(Bio::FastaFormat, path) do |fasta_file|
    fasta_file.each do |entry|
      aln[entry.entry_id] = entry.seq if i < 3
      i += 1
    end
  end
  aln
end


mafft_opts = ['--maxiterate', '1000', '--localpair', '--quiet']
mafft = Bio::MAFFT.new( "mafft" , mafft_opts)
header_printed = false
stats = File.open("#{output_folder}/#{reference_name}.identity_stats.csv", "w")

CSV.foreach(opts[:triads], headers:true ) do |row|
  next unless row["cardinality_abs"] == "1:1:1" and row["HC.LC"] == "HC-only"
   a = row['A']
   b = row['B']
   d = row['D']
   triad = row['group_id']
   cent_triad = triad.to_i / 100
   folder = "#{output_folder}/alignments/#{reference_name}/#{cent_triad}/"
   save_cds = "#{folder}/#{triad}.fa"
   aligned = File.file?(save_cds)
   aln = aligned ? read_alignment(save_cds)  : mafft_align(a,b,d)
   folder = "#{output_folder}/alignments_new/#{reference_name}/#{cent_triad}/" if aligned
   FileUtils.mkdir_p folder
   save_cds = "#{folder}/#{triad}.fa"

   aln2 = Bio::Alignment.new aln
   aln2.add_seq(Bio::PolyploidTools::Mask.get(aln, target: a), "A")
   aln2.add_seq(Bio::PolyploidTools::Mask.get(aln, target: b), "B")
   aln2.add_seq(Bio::PolyploidTools::Mask.get(aln, target: d), "D")

   a_stats =  Bio::PolyploidTools::Mask.stats(aln2["A"], triad, a, "A", reference_name)
   b_stats =  Bio::PolyploidTools::Mask.stats(aln2["B"], triad, b, "B", reference_name)
   d_stats =  Bio::PolyploidTools::Mask.stats(aln2["D"], triad, d, "D", reference_name)
   stats.puts a_stats.keys.join(",") unless header_printed
   stats.puts a_stats.values.join(",")
   stats.puts b_stats.values.join(",")
   stats.puts d_stats.values.join(",")

   header_printed = true


   write_fasta_from_hash(aln2, save_cds)
end

stats.close
