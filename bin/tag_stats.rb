#!/usr/bin/env ruby
require 'optparse'

require 'csv'
require 'fileutils'
require 'tmpdir'
require 'bio-samtools'
require 'bio'
require 'descriptive_statistics'

class Bio::DB::Tag
  def set(str) 
    @tag   = str[0..1]
    @type  = str[3]
    @value = str[5..-1]
    @value = @value.to_i if @type == "i"
  end
end

$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
$: << File.expand_path('.')
path= File.expand_path(File.dirname(__FILE__) + '/../lib/bioruby-polyploid-tools.rb')
require path
opts = {}
opts[:tag] = "NH"
opts[:bam] = nil
opts[:out] = nil
opts[:ref] = nil

out = $stdout

OptionParser.new do |o|
  o.banner = "Usage: tag_stats.rb [options]"

  o.on("-t", "--tag str", "The tag to extract (default NH)") do |o|
    opts[:tag] = o
  end

  o.on("-b", "--bam FILE" , "BAM file with the alignments ") do |o|
    opts[:bam] = o
  end

  o.on("-o", "--out_file CHAR", "File to save the stats") do |o|
    opts[:out] = o
  end
  
  o.on("-r", "--reference FILE", "Fasta file with the reference") do |o|
    opts[:ref] = o
  end
end.parse!

bam =  Bio::DB::Sam.new(fasta: opts[:ref], bam: opts[:bam])
tag = opts[:tag]

sample = File.basename(opts[:bam], '.sorted.bam')
last_ref = ""
values = []
to_print = [:sum, :min, :max, :mean, :mode, :median, :q1, :q2, :q3]
percentiles = [90, 95, 97.5, 99]
#Add the 90, 95, 97.5 and 99 percentiles.
out = File.open(opts[:out], "w")  if opts[:out]
bam.view do |aln |
  if(last_ref != aln.rname)
    
    desc_stats = values.descriptive_statistics
    to_print.each    { |e| out.puts [sample, last_ref, e      , desc_stats[e]       ].join("\t")  } if(last_ref !=  "")
    percentiles.each { |e| out.puts [sample, last_ref, "P#{e}", values.percentile(e)].join("\t")  } if(last_ref !=  "")
    out.puts [sample, last_ref, "N", values.length].join("\t") if(last_ref !=  "")
    values.clear
    last_ref = aln.rname  
  end
  values << aln.tags[tag].value
end

out.close  if opts[:out]