require 'bio'
require 'bio-samtools'
require "set"
require 'systemu'
require 'json'

#require 'bio/db/exonerate'
#require 'bio/db/primer3'
#require 'bio/db/fasta'

puts "Loading all..."

Dir[File.dirname(__FILE__) + "/bio/*.rb"].each {|file| 
 # puts file 
  require_relative file }
Dir[File.dirname(__FILE__) + "/bio/*/*.rb"].each {|file| 
 # puts file
  require_relative file }


#require_relative "bio/BFRTools.rb"
#require_relative "bio/PolyploidTools/ExonContainer.rb"