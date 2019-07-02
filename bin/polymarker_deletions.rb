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


class Bio::PolyploidTools::ExonContainer
   def add_alignments(opts=Hash.new) 
      opts = { :min_identity=>90 }.merge!(opts)
      exonerate_filename = opts[:exonerate_file]
      arm_selection = opts[:arm_selection]

      unless arm_selection
        arm_selection = lambda do | contig_name |
          ret = contig_name[0,3]       
          return ret
        end
      end

      File.open(exonerate_filename) do |f|
        f.each_line do | line |
          record = Bio::DB::Exonerate::Alignment.parse_custom(line)
          if  record and record.identity >= opts[:min_identity]
            snp_array = @snp_map[record.query_id]
            if snp_array != nil
              snp_array.each do |snp|                            
                if snp != nil and snp.position.between?( (record.query_start + 1) , record.query_end)
                  begin
                    exon = record.exon_on_gene_position(snp.position)
                    snp.add_exon(exon, arm_selection.call(record.target_id))
                  rescue Bio::DB::Exonerate::ExonerateException
                    $stderr.puts "Failed for the range #{record.query_start}-#{record.query_end} for position #{snp.position}"
                  end
                end
              end
            end
          end
        end
      end
    end
end

class Bio::DB::Primer3::SNP
  def to_s
     "#{gene}:#{snp_from.chromosome}"
  end
end

class Bio::DB::Primer3::Primer3Record

  def best_pair
    return @best_pair if @best_pair
    @best_pair = nil
    @total_caps = 100
    @primerPairs.each do | primer |
      capital_count = "#{primer.left.sequence}#{primer.right.sequence}".scan(/[A-Z]/).length
      if @best_pair.nil?
        @best_pair = primer 
        @total_caps = capital_count
        next
      end
      if capital_count < @total_caps
        @best_pair = primer 
        @total_caps = capital_count
      end
      if primer.size < @best_pair.size 
        @best_pair = primer 
        @total_caps = capital_count
      end
    end
    #@best_pair = @primerPairs.min
    @best_pair
  end

#CL3339Contig1:T509C AvocetS chromosome_specific exon 4D forward 
  def parse_header
    @snp, @line, @type, @in, @polymorphism, @chromosome, @orientation   = self.sequence_id.split(" ")  
    @type = @type.to_sym
    if @in
      @in = @in.to_sym == :exon 
    else
      @exon = false
    end

    if @polymorphism.to_sym == :homoeologous
      @homoeologous = true
    else
      @homoeologous = false
    end
    @parsed = true
    @orientation = @orientation.to_sym
  end

  def score
    best_pair
#    puts "score"
 #   puts self.inspect
    ret = 0
    ret += @scores[type]
    ret += @scores[:exon] if exon?
    ret -= @total_caps * 10  
    ret -= product_length
    ret
  end

  def to_s
      "#{gene}:#{snp_from.chromosome}"
  end

   def left_primer_snp(snp)
      tmp_primer = String.new(left_primer)
      return tmp_primer
    end

end

arm_selection_functions = Hash.new;


arm_selection_functions[:arm_selection_first_two] = lambda do | contig_name |
  ret = contig_name[0,2]       
  return ret
end

#Function to parse stuff like: "IWGSC_CSS_1AL_scaff_110"
#Or the first two characters in the contig name, to deal with 
#pseudomolecules that start with headers like: "1A"
#And with the cases when 3B is named with the prefix: v443
arm_selection_functions[:arm_selection_embl] = lambda do | contig_name|
  
  arr = contig_name.split('_')
  ret = "U"
  ret = arr[2][0,2] if arr.size >= 3
  ret = "3B" if arr.size == 2 and arr[0] == "v443"
  ret = arr[0][0,2] if arr.size == 1   
  return ret
end

arm_selection_functions[:arm_selection_morex] = lambda do | contig_name |
  ret = contig_name.split(':')[0].split("_")[1];       
  return ret
end

arm_selection_functions[:scaffold] = lambda do | contig_name |
  ret = contig_name;       
  return ret
end

markers = nil

options = {}
options[:model] = "est2genome"
options[:min_identity] = 90
options[:extract_found_contigs] = false
options[:arm_selection] = arm_selection_functions[:arm_selection_embl] ;
options[:genomes_count] = 3


options[:primer_3_preferences] = {
      :primer_product_size_range => "50-150" ,
      :primer_max_size => 25 , 
      :primer_lib_ambiguity_codes_consensus => 1,
      :primer_liberal_base => 1, 
      :primer_num_return=>5,
      :primer_explain_flag => 1,
      :primer_thermodynamic_parameters_path=>File.expand_path(File.dirname(__FILE__) + '../../conf/primer3_config/') + '/'
  }


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

  opts.on("-g", "--genomes_count INT", "Number of genomes (default 3, for hexaploid)") do |o|
    options[:genomes_count] = o.to_i
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
primer_3_input="#{output_folder}/primer_3_input_temp"
primer_3_output="#{output_folder}/primer_3_output_temp"
exons_filename="#{output_folder}/exons_genes_and_contigs.fa"
output_primers="#{output_folder}/primers.csv"
output_to_order="#{output_folder}/primers_to_order.csv"

fasta_file = Bio::DB::Fasta::FastaFile.new({:fasta=>reference})
fasta_file.load_fai_entries

original_name="A"
snp_in="B"

 arm_selection = options[:arm_selection]

unless arm_selection
   arm_selection = lambda do | contig_name |
      ret = contig_name[0,3]       
      return ret
    end
end
begin
log "Reading exons"
exons = Array.new
Bio::FlatFile.auto(sequences) do |ff|
  ff.each do |entry|
    fields = Array.new
    fields << entry.definition
    fields << arm_selection.call(entry.definition)
    fields << entry.seq
    
    line = fields.join(",")
    snp =  Bio::PolyploidTools::NoSNPSequence.parse(line)
    snp.genomes_count = options[:genomes_count]
    exons << snp
     
  end
end



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
container.gene_models(sequences)
container.chromosomes(temp_contigs)
container.add_parental({:name=>"A"})
container.add_parental({:name=>"B"})
exons.each do |exon|
  exon.container = container
  exon.flanking_size = 50
  exon.variation_free_region = options[:variation_free_region]
#  puts exon.inspect
  container.add_snp(exon)

end
container.add_alignments({:exonerate_file=>exonerate_file, :arm_selection=>options[:arm_selection] , :min_identity=>min_identity})

#4.1 generating primer3 file
log "Running primer3"
file = File.open(exons_filename, "w")
container.print_fasta_snp_exones(file)
file.close

file = File.open(primer_3_input, "w")

Bio::DB::Primer3.prepare_input_file(file, options[:primer_3_preferences])
added_exons = container.print_primer_3_exons(file, nil, snp_in)
file.close

Bio::DB::Primer3.run({:in=>primer_3_input, :out=>primer_3_output}) if added_exons > 0

#5. Pick the best primer and make the primer3 output
log "Selecting best primers"
kasp_container=Bio::DB::Primer3::KASPContainer.new
kasp_container.line_1= original_name
kasp_container.line_2= snp_in

if options[:scoring] == :het_dels
  kasp_container.scores = Hash.new
  kasp_container.scores[:chromosome_specific] = 0
  kasp_container.scores[:chromosome_semispecific] = 1000
  kasp_container.scores[:chromosome_nonspecific] = 100    
end

exons.each do |snp|
  snpk = kasp_container.add_snp(snp) 
end

kasp_container.add_primers_file(primer_3_output) if added_exons > 0
header = "Marker,SNP,RegionSize,chromosome,total_contigs,contig_regions,SNP_type,#{original_name},#{snp_in},common,primer_type,orientation,#{original_name}_TM,#{snp_in}_TM,common_TM,selected_from,product_size,errors"
File.open(output_primers, 'w') { |f| f.write("#{header}\n#{kasp_container.print_primers}") }

kasp_container.snp_hash.each_pair do |name, kaspSNP|  
  #puts kaspSNP.snp_from.surrounding_exon_sequences.inspect
  #puts kaspSNP.first_product
  #puts kaspSNP.realigned_primers

  out_fasta_products = "#{output_folder}/#{name}.fa"
  File.open(out_fasta_products, 'w') { |f| f.write(kaspSNP.realigned_primers_fasta) }


end

File.open(output_to_order, "w") { |io|  io.write(kasp_container.print_primers_with_tails()) }

log "DONE"
rescue StandardError => e
  log "ERROR\t#{e.message}"
  $stderr.puts e.backtrace
  raise e 
rescue Exception => e
  log "ERROR\t#{e.message}"
  $stderr.puts e.backtrace
  raise e  
end
#puts container.inspect 

#container.snp_map.each do | gene, snp_array|
#  snp_array.each do |e|
 #   puts e.inspect
#    puts e.aligned_sequences_fasta
#  end
#end

