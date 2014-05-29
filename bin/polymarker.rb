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

options = {}
options[:path_to_contigs] = "/tgac/references/external/projects/iwgsc/css/IWGSC_CSS_all_scaff_v1.fa"
options[:chunks] = 1
options[:bucket_size] = 0
options[:bucket] = 1
options[:model] = "est2genome"
options[:arm_selection] = arm_selection_functions[:arm_selection_embl] ;
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
  
    
end.parse!

p options
p ARGV


#TODO: Use temporary files somewhere in the file system and add traps to delete them/forward them as a result. 
#TODO: Make all this parameters

path_to_contigs=options[:path_to_contigs]

snp_in="A"
original_name="B"
fasta_reference = nil
#test_file="/Users/ramirezr/Dropbox/JIC/PrimersToTest/test_primers_nick_and_james_1.csv"
test_file=options[:marker_list]
test_file=options[:snp_list] if options[:snp_list]
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

write_status "Loading Reference"
#0. Load the fasta index 
fasta_reference_db = nil
if fasta_reference
  fasta_reference_db = Bio::DB::Fasta::FastaFile.new({:fasta=>fasta_reference})
  fasta_reference_db.load_fai_entries
  p "Fasta reference: #{fasta_reference}"
end




#1. Read all the SNP files 
#All the SNPs should be on the same chromosome as the first SNP. 
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
      region = fasta_reference_db.index.region_for_entry(snp.gene).get_full_region
      snp.template_sequence = fasta_reference_db.fetch_sequence(region)
    else
      rise Bio::DB::Exonerate::ExonerateException.new "Wrong number of arguments. " 
    end
    rise Bio::DB::Exonerate::ExonerateException.new "No SNP for line '#{line}'" if snp == nil
    snp.snp_in = snp_in
    snp.original_name = original_name
    snps << snp
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
contigs_f = File.open(temp_contigs, "w")
filename=path_to_contigs 
puts filename
target=filename

fasta_file = Bio::DB::Fasta::FastaFile.new({:fasta=>target})
fasta_file.load_fai_entries

found_cointigs = Set.new
Bio::DB::Exonerate.align({:query=>temp_fasta_query, :target=>target, :model=>model}) do |aln|
  if aln.identity > min_identity
    exo_f.puts aln.line
    unless found_cointigs.include?(aln.target_id) #We only add once each contig. Should reduce the size of the output file. 
      found_cointigs.add(aln.target_id)
      entry = fasta_file.index.region_for_entry(aln.target_id)
      raise ExonerateException.new,  "Entry not found! #{aln.target_id}. Make sure that the #{target_id}.fai was generated properly." if entry == nil
      region = entry.get_full_region
      seq = fasta_file.fetch_sequence(region)
      contigs_f.puts(">#{aln.target_id}\n#{seq}")
    end
  end  
end
 
exo_f.close()
contigs_f.close()

#4. Load all the results from exonerate and get the input filename for primer3
#Custom arm selection function that only uses the first two characters. Maybe
#we want to make it a bit more cleaver
write_status "Reading best alignment on each chromosome"


container= Bio::PolyploidTools::ExonContainer.new
container.flanking_size=100
container.gene_models(temp_fasta_query)
container.chromosomes(temp_contigs)
container.add_parental({:name=>snp_in})
container.add_parental({:name=>original_name})
snps.each do |snp|
  snp.container = container
  snp.flanking_size = container.flanking_size
  container.add_snp(snp)
end
container.add_alignments({:exonerate_file=>exonerate_file, :arm_selection=>options[:arm_selection] , :min_identity=>min_identity})


#4.1 generating primer3 file
write_status "Running primer3"
file = File.open(exons_filename, "w")
container.print_fasta_snp_exones(file)
file.close

file = File.open(primer_3_input, "w")
file.puts("PRIMER_PRODUCT_SIZE_RANGE=50-150")
file.puts("PRIMER_MAX_SIZE=25")
file.puts("PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS=1")
file.puts("PRIMER_LIBERAL_BASE=1")
file.puts("PRIMER_NUM_RETURN=5")
file.puts("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=#{primer_3_config}/")
container.print_primer_3_exons(file, nil, snp_in)
file.close

Bio::DB::Primer3.run({:in=>primer_3_input, :out=>primer_3_output})


#5. Pick the best primer and make the primer3 output
write_status "Selecting best primers"
kasp_container=Bio::DB::Primer3::KASPContainer.new
kasp_container.line_1=snp_in
kasp_container.line_2=original_name

snps.each do |snp|
  kasp_container.add_snp(snp) 
end

kasp_container.add_primers_file(primer_3_output)
header = "Marker,SNP,RegionSize,chromosome,total_contigs,contig_regions,SNP_type,#{snp_in},#{original_name},common,primer_type,orientation,#{snp_in}_TM,#{original_name}_TM,common_TM,selected_from,product_size"
File.open(output_primers, 'w') { |f| f.write("#{header}\n#{kasp_container.print_primers}") }

write_status "DONE"