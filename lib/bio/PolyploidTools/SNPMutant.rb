
require_relative "SNPSequence"
require 'bio-samtools'
module Bio::PolyploidTools
  class SNPSequenceException < RuntimeError 
  end

  class SNPMutant < SNPSequence

    attr_accessor :library, :contig, :chr, :parsed_start, :parsed_flanking, :region_size
    #Format: 
    #seqid,library,position,wt_base,mut_base
    #IWGSC_CSS_1AL_scaff_1455974,Kronos2281,127,C,T
    def self.parse(reg_str)
      reg_str.chomp!
      snp = SNPMutant.new

      arr = reg_str.split(",")
      
      throw SNPSequenceException.new "Need five fields to parse, and got #{arr.size} in #{reg_str}" if arr.size < 5
      
      snp.contig, snp.library, snp.position, snp.original, snp.snp, parsed_flanking, region_size = reg_str.split(",")
      snp.position = snp.position.to_i
      snp.gene = "EMPTY"
      begin
        toks = snp.contig.split('_')
        #1AL_1455974_Kronos2281_127C
        #snp.chr = contig.split('_')[2][0,2] #This parses the default from the IWGSC. We may want to make this a lambda  
        #snp.chr = toks[2][0,2]
        name = toks[2] + "_" + toks[4] + "_" + snp.library + "_" + snp.position.to_s 
        snp.gene = name
        snp.chromosome = toks[2][0,2]
        snp.chr = snp.chromosome
        
      rescue Exception => e
        $stderr.puts "WARN: snp.chr couldnt be set, the sequence id to parse was #{snp.contig}. We expect something like: IWGSC_CSS_1AL_scaff_1455974"
        snp.gene = "Error"
        $stderr.puts e
      end
      
      snp.exon_list = Hash.new()
      snp.flanking_size=100
      snp.region_size = region_size.to_i if region_size
      snp.flanking_size = parsed_flanking.to_i if parsed_flanking
      snp
    end
    
    def full_sequence=(seq)
      self.template_sequence = seq
#       puts self.inspect
       puts self.contig
#       puts self.region_size

      self.sequence_original = self.to_polymarker_sequence(self.flanking_size, total:region_size)
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