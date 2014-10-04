#!/usr/bin/env ruby
require 'bio'

def load_blat_alignments (blat_filename, best_aln)
  blat_aln = Bio::Blat::Report.new(Bio::FlatFile.open(blat_filename).to_io)
  blat_aln.each_hit() do |hit|
    current_matches = hit.match 
    current_name = hit.query_id
    current_identity = hit.percent_identity
    current_score = hit.score
    #p current_name

    best = best_aln[current_name]

    if best == nil 
      best_aln[current_name] = hit
    else
      if current_score > best.score
        best_aln[current_name] = hit
      end
    end
  end
end

blat_file=ARGV[0]
best_aln = Hash.new

load_blat_alignments( blat_file,best_aln)
puts "QUERY\tTARGET"
best_aln.each do |k, hit|
  #puts "#{k}\t#{hit.target_id}"
   puts hit.data.join("\t")
end