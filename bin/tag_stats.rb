#!/usr/bin/env ruby
require 'optparse'

require 'csv'
require 'fileutils'
require 'tmpdir'
require 'bio-samtools'
require 'bio'


$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
$: << File.expand_path('.')
path= File.expand_path(File.dirname(__FILE__) + '/../lib/bioruby-polyploid-tools.rb')
require path
opts = {}
opts[:tag] = "NH"
opts[:bam] = nil
opts[:out] = nil

OptionParser.new do |o|
  o.banner = "Usage: tag_stats.rb [options]"

  o.on("-t", "--tag str", "CSV file with the gene triad names in the named columns 'A','B' and 'D' ") do |o|
    opts[:tag] = o
  end

  o.on("-b", "--bam FILE" , "FASTA file containing all the possible peptide sequences. ") do |o|
    opts[:bam] = o
  end

  o.on("-o", "--out_file CHAR", "Character used to split the sequence name. The name will be evarything before this token on the name of the sequences") do |o|
    opts[:out_file] = o
  end
end


