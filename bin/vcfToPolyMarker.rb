#!/usr/bin/env ruby

require 'optparse'

require 'csv'
require 'bio'
require 'bio-samtools'

$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
$: << File.expand_path('.')
path= File.expand_path(File.dirname(__FILE__) + '/../lib/bioruby-polyploid-tools.rb')
require path

options = {}
options[:arm_selection] = Bio::PolyploidTools::ChromosomeArm.getArmSelection("nrgene");


OptionParser.new do |opts|
	opts.banner = "Usage: polymarker.rb [options]"

	opts.on("-c", "--reference FILE", "File with genome reference to use as database") do |o|
		options[:path_to_contigs] = o
	end

	opts.on("-a", "--arm_selection #{Bio::PolyploidTools::ChromosomeArm.getValidFunctions.join('|')}", "Function to decide the chromome arm") do |o|
		tmp_str = o
		arr = o.split(",")
		if arr.size == 2
			options[:arm_selection] = lambda do |contig_name|
				separator, field = arr
				field = field.to_i
				ret = contig_name.split(separator)[field]
				return ret
			end
		else
			options[:arm_selection] = Bio::PolyploidTools::ChromosomeArm.getArmSelection(o)
		end
	end


end.parse!

def parseVCFheader(head_line="")
	##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">

	m=/##INFO=<ID=(.+),Number=(.+),Type=(.+),Description="(.+)">/.match(head_line)
	{:id=>m[1],:number=>m[2],:type=>m[3],:desc=>m[4]}

end




header_info = Hash.new
ref=options[:path_to_contigs]

fasta_reference_db = Bio::DB::Fasta::FastaFile.new({:fasta=>ref})
fasta_reference_db.load_fai_entries 

$stdin.each do |line|  

	h = nil
	h = parseVCFheader(line) if line.start_with? "##INFO"

	header_info[h[:id]] = h[:desc] if h
	#puts header_info.inspect
	next if line.start_with? "##"
	if line.start_with? "#CHROM"
		arr = line.split
		arr = arr.drop(9)
		arr2 = arr.map { |s| [s.clone().prepend('Cov'), s.clone().prepend('Hap') ]}
		#header += arr2.join("\t")
		#puts header
		next
	end
	line.chomp!
	#puts line
	snp = Bio::PolyploidTools::SNP.parseVCF( line , options[:arm_selection])
	#puts snp.inspect
	snp.setTemplateFromFastaFile(fasta_reference_db, flanking_size: 100)
	puts [snp.gene, snp.chromosome ,snp.to_polymarker_sequence(100)].join(",")
end