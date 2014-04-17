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


#@snp_map=Hash.new

class HomokaryotContainer < Bio::PolyploidTools::ExonContainer
  
  
  def add_snp_file(filename, chromosome, snp_in, original_name)
    flanking_size = 100
     File.open(filename) do | f |
       f.each_line do | line |
         snp = Bio::PolyploidTools::SNP.parse(line)
         snp.flanking_size = flanking_size
         if snp.position > 0
           snp.container = self
           snp.chromosome = chromosome
           snp.snp_in = snp_in
           snp.original_name = original_name
           snp.use_reference = true
           snp.container = self
           @snp_map[snp.gene] = Array.new unless   @snp_map[snp.gene] 
           @snp_map[snp.gene] << snp   
         end
       end
     end
     
     
  end
  
  def print_primer_3_exons (file, target_chromosome , parental )
    @snp_map.each do | gene, snp_array|
      snp_array.each do |snp|
        string = snp.primer_3_string( snp.chromosome, parental )
        file.puts string if string.size > 0

      end 
    end
  end
end

class Bio::PolyploidTools::SNP
  
  @aligned = false
  
  def aligned_snp_position
     return local_position
     
  end
  
  def aligned_sequences
      
    @aligned_sequences = parental_sequences    
    @aligned_sequences["A"][local_position] = original
    @aligned_sequences["B"][local_position] = snp
    return @aligned_sequences
  end
end





snp_file = ARGV[0]
reference_file = ARGV[1]

snp_in="A"
original_name="B"
snps = Array.new

#0. Load the fasta index 
fasta_reference_db = nil
if reference_file
  fasta_reference_db = Bio::DB::Fasta::FastaFile.new({:fasta=>reference_file})
  fasta_reference_db.load_fai_entries
  p "Fasta reference: #{reference_file}"
end
#1. Read all the SNP files 
#All the SNPs should be on the same chromosome as the first SNP. 
chromosome = nil
File.open(snp_file) do | f |
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


container = HomokaryotContainer.new
container.add_parental({:name=>snp_in})
container.add_parental({:name=>original_name})
container.gene_models(reference_file)

output_folder="#{snp_file}_primer_design_#{Time.now.strftime('%Y%m%d-%H%M%S')}/"
Dir.mkdir(output_folder)
primer_3_input="#{output_folder}primer_3_input_temp"
primer_3_output="#{output_folder}primer_3_output_temp"
container.add_snp_file(snp_file, "PST130",  snp_in, original_name)
primer_3_config=File.expand_path(File.dirname(__FILE__) + '/../conf/primer3_config')
output_primers="#{output_folder}primers.csv"

file = File.open(primer_3_input, "w")
file.puts("PRIMER_PRODUCT_SIZE_RANGE=50-150")
file.puts("PRIMER_MAX_SIZE=25")
file.puts("PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS=1")
file.puts("PRIMER_LIBERAL_BASE=1")
file.puts("PRIMER_NUM_RETURN=5")
file.puts("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=#{primer_3_config}/")


container.print_primer_3_exons(file, "PST130",snp_in)

file.close


Bio::DB::Primer3.run({:in=>primer_3_input, :out=>primer_3_output})

#2. Pick the best primer and make the primer3 output
kasp_container=Bio::DB::Primer3::KASPContainer.new
kasp_container.line_1=original_name
kasp_container.line_2=snp_in

snps.each do |snp|
  kasp_container.add_snp(snp) 
end

kasp_container.add_primers_file(primer_3_output)
header = "Marker,SNP,RegionSize,SNP_type,#{snp_in},#{original_name},common,primer_type,orientation,#{snp_in}_TM,#{original_name}_TM,common_TM,selected_from,product_size"
File.open(output_primers, 'w') { |f| f.write("#{header}\n#{kasp_container.print_primers}") }