#!/usr/bin/env ruby
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
OptionParser.new do |opts|
  
  opts.banner = "Usage: polymarker.rb [options]"

  opts.on("-c", "--chromosome CHR", "chromosome (1A, 3B, etc)") do |o|
    options[:chromosome] = o.upcase
  end
  opts.on("-r", "--reference FASTA", "reference with the contigs") do |o|
    options[:reference] = o
  end
  opts.on("-m", "--map CSV", "File with the map and sequence \n Header: INDEX_90K,SNP_ID,SNP_NAME,CHR,COORDINATES_CHR,MAP_ORDER,CHR_ARM,DISTANCE_CM,SEQUENCE") do |o|
    options[:map] = o
  end
  
end.parse!
#reference="/Users/ramirezr/Documents/TGAC/references/Triticum_aestivum.IWGSP1.21.dna_rm.genome.fa"
reference = options[:reference] if options[:reference]
throw raise Exception.new(), "Reference has to be provided" unless reference

map = Bio::PolyploidTools::ArmMap.new
map.chromosome = options[:chromosome]
map.global_reference(reference)
log "Reading markers file"
Bio::PolyploidTools::Marker.parse(options[:map]) do |marker|
 if options[:chromosome] == marker.chr
    map.markers[marker.snp_name] = marker 
  end
end



fasta_tmp="markers_#{options[:chromosome]}.fa"
contigs_tmp="contigs_#{options[:chromosome]}.fa"
aln_tmp="align_#{options[:chromosome]}.psl"
contigs_map="contigs_map_#{options[:chromosome]}.fa"
map_with_contigs="contigs_map_#{options[:chromosome]}.csv"

#1. Prints the sequences to print according to the chromosome to search
log "Writing markers: #{fasta_tmp}"
map.print_fasta_markers(fasta_tmp)
log "Writing contigs: #{contigs_tmp}"
map.print_fasta_contigs_from_reference(contigs_tmp)
log "Aligning markers #{aln_tmp}"
map.align_markers(aln_tmp)
log "printing contigs with markers #{contigs_map}"
map.print_fasta_contigs_for_markers(contigs_map)
log "printing map with contigs #{map_with_contigs}"
map.print_map_with_contigs(map_with_contigs)
