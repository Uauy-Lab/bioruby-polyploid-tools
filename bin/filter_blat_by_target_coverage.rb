#!/usr/bin/env ruby
require 'bio'
require 'optparse'
$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
$: << File.expand_path('.')
path= File.expand_path(File.dirname(__FILE__) + '/../lib/bioruby-polyploid-tools.rb')
require path
module Bio
	class Blat
  		class StreamedReport < Report 

  			def self.each_hit(text = '')
  				flag = false
  				head = []

  				text.each_line do |line|
  					if flag then
  						yield Hit.new(line)
  					else
            		# for headerless data
            			if /^\d/ =~ line then
            				flag = true
            				redo
            			end
            			line = line.chomp
            			if /\A\-+\s*\z/ =~ line
            				flag = true
            			else
            				head << line
            			end
            		end
            	end
            end
        end
    end
end


#blat_file=ARGV[0]

options = {}
options[:identity] = 95
options[:covered] = 60
OptionParser.new do |opts|
  
  opts.banner = "Usage: filter_blat_by_target_coverage.rb [options]"

  opts.on("-p", "--psl FILE", "PSL file") do |o|
    options[:blat_file] = o.upcase
  end
  opts.on("-i", "--identity FLOAT", "Minimum percentage identity") do |o|
    options[:identity] = o.to_f
  end
  opts.on("-c", "--covered FLOAT", "Minimum percentage coverage") do |o|
    options[:covered] = o.to_f
  end
  
end.parse!


blat_file = options[:blat_file]

Bio::Blat::StreamedReport.each_hit(Bio::FlatFile.open(blat_file).to_io) do |hit|
  if hit.percentage_covered >= options[:covered] and hit.percent_identity >= options[:identity]
    puts hit.data.join("\t")
  end
end


