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
options[:debug]  = false

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

    opts.on("-D", "--debug", "Validate that the flanking sequences are correct") do 
    options[:debug]  = true
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
vcf_file="#{output_folder}/snp_positions.vcf"

min_identity= options[:min_identity]

@status_file="#{output_folder}/status.txt"


def write_status(status)
  f=File.open(@status_file, "a")
  f.puts "#{Time.now.to_s},#{status}"
  f.close
end


snps = Hash.new

fasta_reference_db=nil

#if options[:debug]
write_status "Loading Reference"
fasta_reference_db = Bio::DB::Fasta::FastaFile.new({:fasta=>path_to_contigs})
fasta_reference_db.load_fai_entries
write_status "Fasta reference: #{fasta_reference}"
#end

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

def print_positions(min_identity:90, filter_best:false, exonerate_filename:"test.exo", snps:{}, reference:nil, out:$stdout) 
  marker_count=Hash.new { |h, k| h[k] = 1 }
  File.open(exonerate_filename) do |f|
    f.each_line do | line |
      record = Bio::DB::Exonerate::Alignment.parse_custom(line)
      next unless  record and record.identity >= min_identity
      snp = snps[record.query_id]                           
      next unless snp != nil and snp.position.between?( (record.query_start + 1) , record.query_end)
      begin

        position = record.query_position_on_target(snp.position)
        q_strand = record.query_strand
        t_strand = record.target_strand 
        template = snp.template_sequence

        vulgar = record.exon_on_gene_position(snp.position)
        tr = vulgar.target_region
        qr = vulgar.query_region
        template_pre = template[qr.start - 1 .. snp.position - 1 ]
        tr.orientation == :forward ? tr.end = position : tr.start = position
        region = tr
        target_seq = reference.fetch_sequence(region)
        target_seq[-1] = target_seq[-1].upcase
        ref_base = target_seq[-1]
        ma = ref_base
        alt_base = [snp.snp, snp.original].join(",")

        if snp.original == ref_base
          alt_base = snp.snp
        elsif snp.snp == ref_base
          alt_base = snp.original
        end

        if record.target_strand == :reverse
          alt_base = Bio::Sequence::NA.new(alt_base)
          ref_base = Bio::Sequence::NA.new(ref_base)
          alt_base.complement!.upcase!
          ref_base.complement!.upcase!
        end
        
        info =  ["OR=#{record.target_strand}"]
        info <<  "SC=#{record.score}"
        info <<  "PI=#{record.pi}"
        info <<  "MA=#{ma}"
        info <<  "TS=#{target_seq}"
        vcf_line="#{record.target_id}\t#{position}\t#{record.query_id}.path#{marker_count[record.query_id]}\t#{ref_base}\t#{alt_base}\t#{record.pi}\t.\t#{info.join(";")}"
        snp2 = Bio::PolyploidTools::SNP.parseVCF( vcf_line )
        snp2.setTemplateFromFastaFile(reference)
        seq2=snp2.to_polymarker_sequence(50)
        info << "PS=#{seq2}"
        vcf_line="#{record.target_id}\t#{position}\t#{record.query_id}.path#{marker_count[record.query_id]}\t#{ref_base}\t#{alt_base}\t#{record.pi}\t.\t#{info.join(";")}"
        out.puts(vcf_line)
        
        marker_count[record.query_id] += 1
      rescue Bio::DB::Exonerate::ExonerateException
        $stderr.puts "Failed for the range #{record.query_start}-#{record.query_end} for position #{snp.position}"     
      end
    end
  end
end


write_status "Printing VCF file"
#puts snps.inspect
out = File.open(vcf_file, "w")
out.puts "##fileformat=VCFv4.2"
out.puts "##fileDate=#{Time.now.strftime("%Y%m%d")}"
out.puts "##source=#{$0}"
out.puts "##reference=file://#{options[:path_to_contigs]}"
out.puts "##INFO=<ID=OR,Number=1,Type=String,Description=\"Orientation of the alignment of the marker\">"
out.puts "##INFO=<ID=SC,Number=1,Type=Float,Description=\"Alignment score of the marker\">"
out.puts "##INFO=<ID=PI,Number=1,Type=Float,Description=\"Percentage of identity of the alignment to the marker\">"
out.puts "##INFO=<ID=PS,Number=1,Type=String,Description=\"SNP sequence for PolyMarker\">"
out.puts "##INFO=<ID=MA,Number=1,Type=String,Description=\"Allele based on the original marker sequence\">"
out.puts "##INFO<ID=TS,Number=1,Type=String,Description=\"Target sequence before the SNP from the reference\""
out.puts "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
print_positions(exonerate_filename:exonerate_file, min_identity:95, snps:snps, reference: fasta_reference_db, out:out)
out.close
write_status "DONE"


