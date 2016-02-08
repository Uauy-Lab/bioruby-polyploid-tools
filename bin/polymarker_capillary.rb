#!/usr/bin/env ruby
require 'bio'
require 'bio-samtools'
require 'pathname'
require 'optparse'

$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
$: << File.expand_path('.')
path= File.expand_path(File.dirname(__FILE__) + '/../lib/bioruby-polyploid-tools.rb')
require path

def log(msg)
  time=Time.now.strftime("%Y-%m-%d %H:%M:%S.%L")
  puts "#{time}: #{msg}"
end



reference='wheat_6x_ty_mm_mutations_10mutants_for_validations/scaffolds_with_mm.fa'
markers='wheat_6x_ty_mm_mutations_10mutants_for_validations/CadMulitMap.fa'
output_folder='wheat_6x_ty_mm_mutations_10mutants_for_validations/PolyMarker'

options = Hash.new

options[:primer_3_preferences] = {
  :primer_product_size_range => "100-900" ,
  :primer_max_size => 25 , 
  :primer_lib_ambiguity_codes_consensus => 1,
  :primer_liberal_base => 1, 
  :primer_num_return =>1,
  :primer_explain_flag => 1,
  :primer_thermodynamic_parameters_path=>File.expand_path(File.dirname(__FILE__) + '../../conf/primer3_config/') + '/'
}
options[:genomes_count] = 3

OptionParser.new do |opts|
  opts.banner = "Usage: polymarker_capillary.rb [options]"

  opts.on("-r", "--reference FILE", "Fasta file with the assembly") do |o|
    options[:reference] = o
  end

  opts.on("-m", "--sequences FILE", "Fasta file with the sequences to amplify. the format must be Chromosome:start-end. Chromosome should match the names to the entries in the fasta files as it is used as main target") do |o|
    options[:markers] = o
  end

  opts.on("-o", "--output_folder FOLDER", "Path to a folder where the outputs are going to be stored") do |o|
    options[:output_folder] = o
  end
    opts.on("-g", "--genomes_count INT", "Number of genomes (default 3, for hexaploid)") do |o|
    options[:genomes_count] = o.to_i
  end
end.parse!


#puts options.inspect
reference     = options[:reference]
markers       = options[:markers]
output_folder = options[:output_folder]
log "Output folder: #{output_folder}"
exonerate_file="#{output_folder}/exonerate_tmp.tab"
Dir.mkdir(output_folder)

module Bio::PolyploidTools
  class SequenceToAmplify < SNP

    attr_accessor :sequence_original
    attr_accessor :rstart
    attr_accessor :rend
    #Format: 
    #A fasta entry with the id: contig:start-end
    #The sequence can be prodcued with samtools faidx
    def self.parse(fasta_entry)

      snp = SequenceToAmplify.new
      match_data = /(?<rname>\w*):(?<rstart>\w*)-(?<rend>\w*)/.match(fasta_entry.definition)
      
      rName = Regexp.last_match(:rname)
      rStart =  Regexp.last_match(:rstart).to_i
      rEnd =  Regexp.last_match(:rend).to_i
      snp.gene = fasta_entry.definition
      snp.chromosome=rName
      snp.sequence_original = fasta_entry.seq
      snp.template_sequence = fasta_entry.seq.upcase
      snp.snp_in = "B"
      snp.rstart = rStart
      snp.rend = rEnd

      snp.position   = 100
      snp.original   = snp.sequence_original[snp.position]
      
      tmp =  Bio::Sequence::NA.new(snp.original)
      rev = tmp.complement
      snp.snp = rev
      snp.exon_list = Hash.new()
      snp
    end

    def primer_3_all_strings(target_chromosome, parental) 
      #puts target_chromosome 
      #puts parental
      #puts aligned_sequences.to_fasta
      pr = primer_region(target_chromosome, parental )
      primer_3_propertes = Array.new

      seq_original = String.new(pr.sequence)
      #puts seq_original.size.to_s << "-" << primer_3_min_seq_length.to_s
      return primer_3_propertes if seq_original.size < primer_3_min_seq_length
      #puts "Sequence origina: #{ self.original}" 
      #puts pr.to_fasta
#      puts "Postion: #{pr.snp_pos}"
      seq_original[pr.snp_pos] = self.original
      seq_original_reverse = reverse_complement_string(seq_original)

      seq_snp =  String.new(pr.sequence)
      seq_snp[pr.snp_pos] =  self.snp
      seq_snp_reverse = reverse_complement_string(seq_snp)

      rev_pos = seq_snp.size - position

      if pr.homoeologous
        snp_type = "homoeologous"
      else
        snp_type = "non-homoeologous"
      end
      left_pos = Array.new
      right_pos = Array.new
      l_pos = pr.snp_pos
      pr.chromosome_specific.each {|pos| left_pos  << pos if pos < l_pos - 50 }
      pr.chromosome_specific.each {|pos| right_pos << pos if pos > l_pos + 50}
      
      pr.crhomosome_specific_intron.each {|pos| left_pos  << pos if pos < l_pos - 50}
      pr.crhomosome_specific_intron.each {|pos| right_pos << pos if pos > l_pos + 50}

      left_pos.each do |l|
        right_pos.each do |r|
          args = {:name =>"#{gene}:#{original}#{position}#{snp} #{original_name} chromosome_specific exon #{snp_type} #{chromosome}", 
          :left_pos => l, 
          :right_pos => r, 
          :sequence=>seq_original}
          primer_3_propertes << return_primer_3_string(args)
        end
      end
      primer_3_propertes
    end

    def parental_sequences
      return @parental_sequences if @parental_sequences
      gene_region = self.covered_region
      local_pos_in_gene = self.position
      
      @parental_sequences = Bio::Alignment::SequenceHash.new
      container.parents.each  do |name, bam|
        seq = self.sequence_original.clone.downcase
        
        if name == self.snp_in
          #puts self.snp
          seq[local_pos_in_gene] = self.snp
        else
          #puts self.original
          seq[local_pos_in_gene] = self.original
        end
        seq[local_pos_in_gene] = seq[local_pos_in_gene].upcase    
        @parental_sequences [name] = seq
        #puts name 
        #puts self.snp_in
        #puts seq
      end
      @parental_sequences
    end
  end
end


snps = Array.new
file = Bio::FastaFormat.open(markers)
file.each do |entry|

  begin
    tmp = Bio::PolyploidTools::SequenceToAmplify.parse(entry)
    snps << tmp if tmp
  rescue
    $stderr.puts "Unable to generate the marker for: #{entry.definition}"
  end

end
file.close



exo_f = File.open(exonerate_file, "w")
target=reference

fasta_file = Bio::DB::Fasta::FastaFile.new({:fasta=>target})
fasta_file.load_fai_entries
min_identity = 90
found_contigs = Set.new

Bio::DB::Exonerate.align({:query=>markers, :target=>reference, :model=>'ungapped'}) do |aln|
  if aln.identity > min_identity
    exo_f.puts aln.line
    #puts aln.line
    unless found_contigs.include?(aln.target_id) #We only add once each contig. Should reduce the size of the output file. 
      found_contigs.add(aln.target_id)
      entry = fasta_file.index.region_for_entry(aln.target_id)
      raise Exception.new,  "Entry not found! #{aln.target_id}. Make sure that the #{reference}.fai was generated properly." if entry == nil
    end
  end  
end
exo_f.close

arm_selection_functions = Hash.new

arm_selection_functions[:full_scaffold] = lambda do | contig_name |    
  return contig_name
end

container= Bio::PolyploidTools::ExonContainer.new
container.flanking_size=500 
container.gene_models(markers)
container.chromosomes(target)
container.add_parental({:name=>"A"})
container.add_parental({:name=>"B"})
snps.each do |snp|
  snp.snp_in = "B"
  snp.container = container
  snp.flanking_size = container.flanking_size
  snp.genomes_count = options[:genomes_count]
  container.add_snp(snp)
end
container.add_alignments({:exonerate_file=>exonerate_file, :arm_selection=>arm_selection_functions[:full_scaffold] , :min_identity=>min_identity})



exons_filename="#{output_folder}/localAlignment.fa"
file = File.open(exons_filename, "w")
container.print_fasta_snp_exones(file)
file.close



primer_3_input  ="#{output_folder}/primer3_input.txt"
primer_3_output ="#{output_folder}/primer3_output.txt"



file = File.open(primer_3_input, "w")
snp_in="B"
Bio::DB::Primer3.prepare_input_file(file, options[:primer_3_preferences])
added_exons = container.print_primer_3_exons(file, nil, snp_in)
file.close

Bio::DB::Primer3.run({:in=>primer_3_input, :out=>primer_3_output}) if added_exons > 0

masks_output = "#{output_folder}/masks_designed.fa"
output_file  = "#{output_folder}/primers.csv"
file = File.open(masks_output, "w")
out  = File.open(output_file,  "w")

class Bio::DB::Primer3::Primer3Record
  attr_accessor :primerPairs
end

Bio::DB::Primer3::Primer3Record.parse_file(primer_3_output) do | primer3record |
  #puts primer3record.inspect
  next if primer3record.primer_left_num_returned.to_i == 0
  
  seq_id = primer3record.sequence_id
  excluded = "-"
  exArr = excluded.split(",")
  st = exArr[0].to_i
  ed = exArr[1].to_i
  tot = ed + st

  excluded="#{st}-#{tot}"
  seq_len = primer3record.sequence_template.length
  printed = 0

  sequence_template = primer3record.sequence_template
  sequence_mask = "-" * st 
  sequence_mask << "*" * ed 
  sequence_mask << "-" * (seq_len - sequence_mask.length)

  file.puts ">#{seq_id}\n#{sequence_template}"
  file.puts ">#{seq_id}:mask\n#{sequence_mask}"
   #puts "FDFDS"

   #puts primer3record.primerPairs

   primer3record.primerPairs.each do |p| 
    #puts p.inspect
    printed += 1   
    lArr =  p.left.coordinates
    lArr[1] = lArr[0] + lArr[1]
    rArr =  p.right.coordinates
    rArr[1] = rArr[0] - rArr[1]
    toPrint = Array.new
    toPrint <<  seq_id.split(" ")
    #toPrint <<  seq_len
    toPrint <<  p.product_size
    toPrint <<  lArr.join("-")
    toPrint <<  p.left.tm
    toPrint <<  p.left.sequence
    toPrint <<  rArr.join("-")
    toPrint <<  p.right.tm
    toPrint <<  p.right.sequence

    middle = 501 
    toPrint << lArr[0]
    toPrint << rArr[0]
    toPrint << middle - lArr[0]
    toPrint << rArr[0] - middle
#Start End LeftDistance  RightDistance

    out.puts toPrint.join(",")

    sequence_primers = sequence_mask.clone
    a = lArr[0] 
    b = lArr[1] - 1
    #puts   sequence_template[a..b]
    sequence_primers[a..b] = sequence_template[a..b]
    b = rArr[0] 
    a = rArr[1] + 1

    sequence_primers[a..b] = sequence_template[a..b]

    file.puts ">#{seq_id}:primerPair:#{printed}\n#{sequence_primers}"
  end 

  if printed == 0
    toPrint = Array.new
    toPrint << seq_id.split(" ")
    toPrint << excluded
    toPrint << seq_len
    out.puts toPrint.join(",")
  end

end
out.close
file.close

