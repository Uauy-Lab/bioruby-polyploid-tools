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

arm_selection_functions = Hash.new;


arm_selection_functions[:arm_selection_first_two] = lambda do | contig_name |
  ret = contig_name[0,2]       
  return ret
end
#Function to parse stuff like: IWGSC_CSS_1AL_scaff_110
arm_selection_functions[:arm_selection_embl] = lambda do | contig_name|
  ret = contig_name.split('_')[2][0,2]
  return ret
end

arm_selection_functions[:arm_selection_morex] = lambda do | contig_name |
  ret = contig_name.split(':')[0].split("_")[1];       
  return ret
end

def validate_files(o)

  [
    o[:path_to_contigs], 
    o[:marker_list], 
    o[:snp_list], 
    o[:mutant_list],
    o[:reference]
  ].flatten.compact.each do |f|  
        raise IOError "Unable to read #{f}" unless File.exists? f 
    end
end 

options = {}
options[:path_to_contigs] = "/tgac/references/external/projects/iwgsc/css/IWGSC_CSS_all_scaff_v1.fa"
options[:chunks] = 1
options[:bucket_size] = 0
options[:bucket] = 1
options[:model] = "est2genome"
options[:arm_selection] = arm_selection_functions[:arm_selection_embl] ;
options[:flanking_size] = 150;
options[:variation_free_region] = 0 
options[:extract_found_contigs] = false
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
  opts.banner = "Usage: polymarker.rb [options]"

  opts.on("-c", "--contigs FILE", "File with contigs to use as database") do |o|
    options[:path_to_contigs] = o
  end
  
  opts.on("-m", "--marker_list FILE", "File with the list of markers to search from") do |o|
    options[:marker_list] = o
  end
  
  opts.on("-s", "--snp_list FILE", "File with the list of snps to search from, requires --reference to get the sequence using a position") do |o|
    options[:snp_list] = o
  end

  opts.on("-t", "--mutant_list FILE", "File with the list of positions with mutation and the mutation line.\n\
    requires --reference to get the sequence using a position") do |o|
    options[:mutant_list] = o
  end
  
  opts.on("-r", "--reference FILE", "Fasta file with the sequence for the markers (to complement --snp_list)") do |o|
    options[:reference] = o
  end
  
  opts.on("-o", "--output FOLDER", "Output folder") do |o|
    options[:output_folder] = o
  end
  
  opts.on("-e", "--exonerate_model MODEL", "Model to be used in exonerate to search for the contigs") do |o|
     options[:model] = o
  end

  opts.on("-a", "--arm_selection arm_selection_embl|arm_selection_morex|arm_selection_first_two", "Function to decide the chromome arm") do |o|
    options[:arm_selection] = arm_selection_functions[o.to_sym];
   end
  
  opts.on("-p", "--primer_3_preferences FILE", "file with preferences to be sent to primer3") do |o|
    options[:primer_3_preferences] = Bio::DB::Primer3.read_primer_preferences(o, options[:primer_3_preferences] )
  end

  opts.on("-v", "--variation_free_region INT", "If present, avoid generating the common primer if there are homoeologous SNPs within the specified distance") do |o|
    options[:variation_free_region] = o.to_i
  end

  opts.on("-x", "--extract_found_contigs", "If present, save in a separate file the contigs with matches. Useful to debug.") do |o|
    options[:extract_found_contigs] = true
  end

  opts.on("-P", "--primers_to_order")do
    #TODO: have a string with the tails, optional. 
    options[:primers_to_order] = true
  end

  
end.parse!

validate_files(options)

if options[:primer_3_preferences][:primer_product_size_range]
  range = options[:primer_3_preferences][:primer_product_size_range]
  range_arr = range.split("-")
  min = range_arr[0].to_i
  max = range_arr[1].to_i
  raise  Bio::DB::Exonerate::ExonerateException.new "Range #{range} is invalid!" unless max > min
  options[:flanking_size] = max
end

p options
p ARGV


#TODO: Use temporary files somewhere in the file system and add traps to delete them/forward them as a result. 
#TODO: Make all this parameters

path_to_contigs=options[:path_to_contigs]

original_name="A"
snp_in="B"

fasta_reference = nil
#test_file="/Users/ramirezr/Dropbox/JIC/PrimersToTest/test_primers_nick_and_james_1.csv"
test_file=options[:marker_list]  if options[:marker_list]
test_file=options[:snp_list] if options[:snp_list]
test_file=options[:mutant_list] if options[:mutant_list]
fasta_reference = options[:reference]
output_folder="#{test_file}_primer_design_#{Time.now.strftime('%Y%m%d-%H%M%S')}" 
output_folder= options[:output_folder] if  options[:output_folder]
Dir.mkdir(output_folder)
#TODO Make this tmp files
temp_fasta_query="#{output_folder}/to_align.fa"
temp_contigs="#{output_folder}/contigs_tmp.fa"
exonerate_file="#{output_folder}/exonerate_tmp.tab"
primer_3_input="#{output_folder}/primer_3_input_temp"
primer_3_output="#{output_folder}/primer_3_output_temp"
exons_filename="#{output_folder}/exons_genes_and_contigs.fa"
output_primers="#{output_folder}/primers.csv"
output_to_order="#{output_folder}/primers_to_order.csv"
@status_file="#{output_folder}/status.txt"

primer_3_config=File.expand_path(File.dirname(__FILE__) + '/../conf/primer3_config')
model=options[:model] 

def write_status(status)
  f=File.open(@status_file, "a")
  f.puts "#{Time.now.to_s},#{status}"
  f.close
end

min_identity= 90
snps = Array.new

begin
  
write_status "Loading Reference"
#0. Load the fasta index 
fasta_reference_db = nil
if fasta_reference
  fasta_reference_db = Bio::DB::Fasta::FastaFile.new({:fasta=>fasta_reference})
  fasta_reference_db.load_fai_entries
  p "Fasta reference: #{fasta_reference}"
end

#1. Read all the SNP files 
#chromosome = nil
write_status "Reading SNPs"
File.open(test_file) do | f |
  f.each_line do | line |
    # p line.chomp!
    snp = nil
    if options[:marker_list] #List with Sequence
      snp = Bio::PolyploidTools::SNPSequence.parse(line)  
    elsif options[:snp_list] and options[:reference] #List and fasta file
      snp = Bio::PolyploidTools::SNP.parse(line)
      entry = fasta_reference_db.index.region_for_entry(snp.gene)
      if entry
       region = fasta_reference_db.index.region_for_entry(snp.gene).get_full_region
       snp.template_sequence = fasta_reference_db.fetch_sequence(region)
     else
        write_status "WARN: Unable to find entry for #{snp.gene}"
      end
    elsif options[:mutant_list] and options[:reference] #List and fasta file
      snp = Bio::PolyploidTools::SNPMutant.parse(line)
      entry = fasta_reference_db.index.region_for_entry(snp.contig)
      if entry
       region = fasta_reference_db.index.region_for_entry(snp.contig).get_full_region
       snp.full_sequence = fasta_reference_db.fetch_sequence(region)
     else
        write_status "WARN: Unable to find entry for #{snp.gene}"
      end
    else
      rise Bio::DB::Exonerate::ExonerateException.new "Wrong number of arguments. " 
    end
    rise Bio::DB::Exonerate::ExonerateException.new "No SNP for line '#{line}'" if snp == nil
  
    snp.snp_in = snp_in
    snp.original_name = original_name
    if snp.position 
      snps << snp
    else
      $stderr.puts "ERROR: #{snp.gene} doesn't contain a SNP"
    end

#    chromosome = snp.chromosome unless chromosome
  #  raise Bio::DB::Exonerate::ExonerateException.new "All the snps should come from the same chromosome" if chromosome != snp.chromosome
  end
end

#1.1 Close fasta file
#fasta_reference_db.close() if fasta_reference_db
#2. Generate all the fasta files
write_status "Writing sequences to align"
written_seqs = Set.new
file = File.open(temp_fasta_query, "w")
snps.each do |snp|
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
puts filename
target=filename

fasta_file = Bio::DB::Fasta::FastaFile.new({:fasta=>target})
fasta_file.load_fai_entries

found_contigs = Set.new
Bio::DB::Exonerate.align({:query=>temp_fasta_query, :target=>target, :model=>model}) do |aln|
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

#4. Load all the results from exonerate and get the input filename for primer3
#Custom arm selection function that only uses the first two characters. Maybe
#we want to make it a bit more cleaver
write_status "Reading best alignment on each chromosome"


container= Bio::PolyploidTools::ExonContainer.new
container.flanking_size=options[:flanking_size] 
container.gene_models(temp_fasta_query)
container.chromosomes(target)
container.add_parental({:name=>snp_in})
container.add_parental({:name=>original_name})
snps.each do |snp|
  snp.container = container
  snp.flanking_size = container.flanking_size
  snp.variation_free_region = options[:variation_free_region]
  container.add_snp(snp)
end
container.add_alignments({:exonerate_file=>exonerate_file, :arm_selection=>options[:arm_selection] , :min_identity=>min_identity})


#4.1 generating primer3 file
write_status "Running primer3"
file = File.open(exons_filename, "w")
container.print_fasta_snp_exones(file)
file.close

file = File.open(primer_3_input, "w")

Bio::DB::Primer3.prepare_input_file(file, options[:primer_3_preferences])
added_exons = container.print_primer_3_exons(file, nil, snp_in)
file.close

Bio::DB::Primer3.run({:in=>primer_3_input, :out=>primer_3_output}) if added_exons > 0


#5. Pick the best primer and make the primer3 output
write_status "Selecting best primers"
kasp_container=Bio::DB::Primer3::KASPContainer.new
kasp_container.line_1= original_name
kasp_container.line_2= snp_in

snps.each do |snp|
  kasp_container.add_snp(snp) 
end

kasp_container.add_primers_file(primer_3_output) if added_exons > 0
header = "Marker,SNP,RegionSize,chromosome,total_contigs,contig_regions,SNP_type,#{original_name},#{snp_in},common,primer_type,orientation,#{original_name}_TM,#{snp_in}_TM,common_TM,selected_from,product_size,errors"
File.open(output_primers, 'w') { |f| f.write("#{header}\n#{kasp_container.print_primers}") }

File.open(output_to_order, "w") { |io|  io.write(kasp_container.print_primers_with_tails()) }

write_status "DONE"
rescue StandardError => e
  write_status "ERROR\t#{e.message}"
  raise e 
rescue Exception => e
  write_status "ERROR\t#{e.message}"
  raise e  
end