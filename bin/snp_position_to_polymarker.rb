#!/usr/bin/env ruby

#This This script converts the a file with snps and positions with the header:
#GENE,BASE,POS,SNP,Chromosome
#  snp.gene, snp.original, snp.position, snp.snp, snp.chromosome 
#To the input expected by polymarker
#ID, Chromosome, sequence
#With sequence containing the SNP in the notation "[A/T]"
require 'bio'
require 'optparse'

$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
$: << File.expand_path('.')
path= File.expand_path(File.dirname(__FILE__) + '/../lib/bioruby-polyploid-tools.rb')
require path


def log(msg)
  time=Time.now.strftime("%Y-%m-%d %H:%M:%S.%L")
  puts "#{time}: #{msg}"
end

markers = nil

options = {}
options[:flanking_size] = 100
test_file=''
OptionParser.new do |opts|
  
  opts.banner = "Usage: snp_postion_to_polymarker.rb [options]"

  opts.on("-s", "--snp_file CSV", "CSV file with the following columnns:\nID,Allele_1,position,Allele_1,target_chromosome") do |o|
    options[:snp_file] = o
    test_file = o 
  end
  opts.on("-r", "--reference FASTA", "reference with the genes/contings/marker seuqnece") do |o|
    options[:reference] = o
  end
  opts.on("-o", "--out CSV", "Output file ") do |o|
    options[:output] = o
  end
  opts.on("-f", "--flanking_size INT", "Flanking size around the SNP") do |o|
    options[:flanking_size] = o.to_i
  end

  opts.on("-t", "--mutant_list FILE", "File with the list of positions with mutation and the mutation line. Example: IWGSC_CSS_1AL_scaff_1455974,Kronos2281,127,C,T\n\
    requires --reference to get the sequence using a position") do |o|
    options[:mutant_list] = o
     test_file = o 
  end
  
end.parse!
#reference="/Users/ramirezr/Documents/TGAC/references/Triticum_aestivum.IWGSP1.21.dna_rm.genome.fa"

fasta_reference = options[:reference] if options[:reference]
fasta_reference_db = Bio::DB::Fasta::FastaFile.new({:fasta=>fasta_reference})
fasta_reference_db.load_fai_entries

out = $stdout
lastRegion = nil
lastTemplate = nil
out = File.open(options[:output], "w") if options[:output]
File.open(test_file) do | f |
  f.each_line do | line |
    	snp = nil
      entry = nil
      if options[:snp_file]
    	   snp = Bio::PolyploidTools::SNP.parse(line)
         entry = fasta_reference_db.index.region_for_entry(snp.gene)
      elsif options[:mutant_list]
         snp = Bio::PolyploidTools::SNPMutant.parse(line)
         entry = fasta_reference_db.index.region_for_entry(snp.contig)
      end
    	
    	if entry
       		region = entry.get_full_region
       		if region != lastRegion
             lastTemplate = fasta_reference_db.fetch_sequence(region)
          end
          snp.template_sequence = lastTemplate
          lastRegion = region
       		out.puts "#{snp.gene}_#{snp.snp_id_in_seq},#{snp.chromosome},#{snp.to_polymarker_sequence(options[:flanking_size])}"
    	else
    	   $stderr.puts "ERROR: Unable to find entry for #{snp.gene}"
    	end
	end
end

out.close if options[:output]

