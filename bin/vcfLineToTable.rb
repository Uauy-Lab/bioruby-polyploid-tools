require 'bio-samtools'
require 'optparse'

$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
$: << File.expand_path('.')
path=File.expand_path(File.dirname(__FILE__) + '/../lib/bioruby-polyploid-tools.rb')




def parseVCFheader(head_line="")
	##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">

	m=/##INFO=<ID=(.+),Number=(.+),Type=(.+),Description="(.+)">/.match(head_line)
	{:id=>m[1],:number=>m[2],:type=>m[3],:desc=>m[4]}

end


header_info = Hash.new
ARGF.each_line do |line|  
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
		
	vcf = Bio::DB::Vcf.new(line, arr)
#	puts arr.join("\t") if vcf.info["TYPE"] == "snp"
#	puts vcf.inspect
	#pus vcf.pos.inspect
	#next if vcf.info["AO"].to_i != 1
	vcf.info.each_pair { |name, val| puts "#{name}\t#{val}\t#{header_info[name]}" }

    arr2 = Array.new
    puts "____"
    i = 0
	vcf.samples.each do |sample|
		#puts sample.inspect
		puts sample[1].keys.join("\t") if i == 0
        puts sample[1].values.join("\t")
        i+=1
    end

end
