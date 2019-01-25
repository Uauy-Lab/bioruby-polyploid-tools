require 'bio'
require 'bio-samtools'
require "set"
require 'systemu'
require 'json'

#require 'bio/db/exonerate'
#require 'bio/db/primer3'
#require 'bio/db/fasta'

#puts "Loading all... #{Dir[File.dirname(__FILE__) + "/bio/*/*.rb"]}"
module Bio::PolyploidTools

end
Dir[File.dirname(__FILE__) + "/bio/*.rb"].each {|file| 
 # puts file 
  require_relative file }
Dir[File.dirname(__FILE__) + "/bio/*/*.rb"].each do |file| 
#  $stderr.puts "loading #{file}"
  require_relative file 
end
  
require_relative File.dirname(__FILE__) + "/../conf/defaults.rb"


#require_relative "bio/BFRTools.rb"
#require_relative "bio/PolyploidTools/ExonContainer.rb"

