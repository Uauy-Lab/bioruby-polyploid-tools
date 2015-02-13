
require_relative "SNP"
require 'bio-samtools'
module Bio::PolyploidTools
  class SNPSequenceException < RuntimeError 
  end

  class SNPMutant < SNPSequence

    attr_accessor :library, :chr
    #Format: 
    #seqid,library,position,wt_base,mut_base
    #IWGSC_CSS_1AL_scaff_1455974,Kronos2281,127,C,T
    def self.parse(reg_str)
      reg_str.chomp!
      snp = SNPMutant.new

      arr = reg_str.split(",")
      
      throw SNPSequenceException.new "Need five fields to parse, and got #{arr.size} in #{reg_str}" unless arr.size == 5
      
      snp.gene, snp.library, snp.position, snp.original, snp.snp = reg_str.split(",")
      snp.position = snp.position.to_i
      snp.chromosome = snp.gene
      begin
        snp.chr = snp.gene.split('_')[2][0,2] #This parses the default from the IWGSC. We may want to make this a lambda  
      rescue Exception => e
        $stderr.puts "WARN: snp.chr couldnt be set, the sequence id to parse was #{snp.gene}"
      end
      
      snp.exon_list = Hash.new()
      snp.flanking_size=100
      snp
    end
    
    def full_sequence=(seq)
      self.template_sequence = seq
      self.sequence_original = self.to_polymarker_sequence(self.flanking_size)
      self.parse_sequence_snp
    end

    def full_sequence()
      self.template_sequence
    end

    def chromosome_group
      chr[0]
    end
    
    def chromosome_genome
      chr[1]
    end
    
    def chromosome_genome
      return chr[3] if chr[3]
      return nil
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