#!/usr/bin/env ruby
require 'bio'
require 'rubygems'
require 'pathname'
require 'bio-samtools'

require 'set'

$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
$: << File.expand_path('.')
path= File.expand_path(File.dirname(__FILE__) + '/../lib/bioruby-polyploid-tools.rb')
require path


#TODO: Use temporary files somewhere in the file system and add traps to delete them/forward them as a result. 
#TODO: Make all this parameters
path_to_contigs="/Users/ramirezr/Documents/PHD/201305_Databases/iwgcs"
#path_to_contigs=path_to_chromosomes
snp_in="A"
original_name="B"
fasta_reference = nil
#test_file="/Users/ramirezr/Dropbox/JIC/PrimersToTest/test_primers_nick_and_james_1.csv"
test_file=ARGV[0]
fasta_reference = ARGV[1] if ARGV[1]
output_folder="#{test_file}_primer_design_#{Time.now.strftime('%Y%m%d-%H%M%S')}/"
Dir.mkdir(output_folder)
#TODO Make this tmp files
temp_fasta_query="#{output_folder}to_align.fa"
temp_contigs="#{output_folder}contigs_tmp.fa"
exonerate_file="#{output_folder}exonerate_tmp.tab"
primer_3_input="#{output_folder}primer_3_input_temp"
primer_3_output="#{output_folder}primer_3_output_temp"
exons_filename="#{output_folder}exons_genes_and_contigs.fa"
output_primers="#{output_folder}primers.csv"

primer_3_config=File.expand_path(File.dirname(__FILE__) + '/../conf/primer3_config')
model="est2genome"


min_identity= 92
snps = Array.new

#0. Load the fasta index 
fasta_reference_db = nil
if fasta_reference
  fasta_reference_db = Bio::DB::Fasta::FastaFile.new(fasta_reference)
  fasta_reference_db.load_fai_entries
  p "Fasta reference: #{fasta_reference}"
end


#1. Read all the SNP files 
#All the SNPs should be on the same chromosome as the first SNP. 
chromosome = nil
File.open(test_file) do | f |
  f.each_line do | line |
    # p line.chomp!
    snp = nil
    if ARGV.size == 1 #List with Sequence
      snp = Bio::PolyploidTools::SNPSequence.parse(line)  
    elsif ARGV.size == 2 #List and fasta file
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
    chromosome = snp.chromosome unless chromosome
    raise Bio::DB::Exonerate::ExonerateException.new "All the snps should come from the same chromosome" if chromosome != snp.chromosome
  end
end

#1.1 Close fasta file
#fasta_reference_db.close() if fasta_reference_db
#2. Generate all the fasta files

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
puts chromosome
chr_group = chromosome[0]
exo_f = File.open(exonerate_file, "w")
contigs_f = File.open(temp_contigs, "w")
Dir.foreach(path_to_contigs) do |filename |
  #puts filename
  if  File.fnmatch("#{chr_group}*.fa", filename)
    puts filename
    target="#{path_to_contigs}/#{filename}"
    
    fasta_file = Bio::DB::Fasta::FastaFile.new(target)
    fasta_file.load_fai_entries
    Bio::DB::Exonerate.align({:query=>temp_fasta_query, :target=>target, :model=>model}) do |aln|
      if aln.identity > min_identity
        exo_f.puts aln.line
        region = fasta_file.index.region_for_entry(aln.target_id).get_full_region
        seq = fasta_file.fetch_sequence(region)
        contigs_f.puts(">#{aln.target_id}\n#{seq}")
      end
      
    end
  end
end

exo_f.close()
contigs_f.close()

#4. Load all the results from exonerate and get the input filename for primer3
#Custom arm selection function that only uses the first two characters. Maybe
#we want to make it a bit more cleaver
arm_selection = lambda do | contig_name |
  ret = contig_name[0,2]       
  return ret
end

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
container.add_alignments({:exonerate_file=>exonerate_file, :arm_selection=>arm_selection})

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
container.print_primer_3_exons(file, chromosome,snp_in)
file.close

Bio::DB::Primer3.run({:in=>primer_3_input, :out=>primer_3_output})

#5. Pick the best primer and make the primer3 output
kasp_container=Bio::DB::Primer3::KASPContainer.new
kasp_container.line_1=snp_in
kasp_container.line_2=original_name

snps.each do |snp|
  kasp_container.add_snp(snp) 
end

kasp_container.add_primers_file(primer_3_output)
header = "Marker,SNP,RegionSize,SNP_type,#{snp_in},#{original_name},common,primer_type,orientation,#{snp_in}_TM,#{original_name}_TM,common_TM,selected_from"
File.open(output_primers, 'w') { |f| f.write("#{header}\n#{kasp_container.print_primers}") }

