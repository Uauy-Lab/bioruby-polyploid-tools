
require_relative "SNP"
require 'bio-samtools'
module Bio::PolyploidTools
  class SNPSequenceException < RuntimeError 
  end

  class SNPSequence < SNP

    attr_accessor :sequence_original
    #Format: 
    #snp name,chromsome from contig,microarray sequence
    #BS00068396_51,2AS,CGAAGCGATCCTACTACATTGCGTTCCTTTCCCACTCCCAGGTCCCCCTA[T/C]ATGCAGGATCTTGATTAGTCGTGTGAACAACTGAAATTTGAGCGCCACAA
    def self.parse(reg_str)
      reg_str.chomp!
      snp = SNPSequence.new

      arr = reg_str.split(",")
      
      if arr.size == 3
        snp.gene, snp.chromosome, snp.sequence_original = arr
      elsif arr.size == 2
       snp.gene, snp.sequence_original = arr
       snp.chromosome = ""
     else
       throw SNPSequenceException.new "Need two or three fields to parse, and got #{arr.size} in #{reg_str}"
      end
      #snp.position = snp.position.to_i
      #snp.original.upcase!
      #snp.snp.upcase!  
      snp.chromosome. strip!
      snp.parse_sequence_snp
  
      snp
    end
    
    def parse_sequence_snp
      pos = 0
      match_data = /(?<pre>\w*)\[(?<org>[ACGT])\/(?<snp>[ACGT])\](?<pos>\w*)/.match(sequence_original.strip)
      if match_data
        @position = Regexp.last_match(:pre).size + 1
        @original = Regexp.last_match(:org)
        @snp = Regexp.last_match(:snp)
        
        amb_base = Bio::NucleicAcid.to_IUAPC("#{@original}#{@snp}")
        
        @template_sequence = "#{Regexp.last_match(:pre)}#{amb_base}#{Regexp.last_match(:pos)}"
        
     end 
    end
    
   

  end
end
