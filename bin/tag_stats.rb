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

OptionParser.new do |o|
  o.banner = "Usage: tag_stats.rb [options]"

  o.on("-t", "--tag str", "The tag to extract (default NH)") do |o|
    opts[:tag] = o
  end

  o.on("-b", "--bam FILE" , "BAM file with the alignments ") do |o|
    opts[:bam] = o
  end

  o.on("-o", "--out_file CHAR", "File to save the stats") do |o|
    opts[:out_file] = o
  end
  o.on("-r", "--reference FILE", "Fasta file with the reference") do |o|
    opts[:ref] = o
  end
end.parse!

bam =  Bio::DB::Sam.new(fasta: opts[:ref], bam: opts[:bam])
tag = opts[:tag]
puts bam.inspect 

last_ref = ""
values = []
to_print = [:sum, :min, :max, :mean, :mode, :median, :q1, :q2, :q3]
bam.view do |aln |
  if(last_ref != aln.rname)
    
    desc_stats = values.descriptive_statistics
    to_print.each { |e| puts "#{last_ref}\t#{e}\t#{desc_stats[e]}" } if(last_ref !=  "")
   
    values.clear
    last_ref = aln.rname  
  end
  
  values << aln.tags[tag].value
	
end