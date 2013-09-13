module Bio::PolyploidTools
  class SNPSequenceException < RuntimeError 
  end

  class SNPSequence < SNP

    attr_accessor :sequence_original, :template_sequence
    #Format: 
    #snp name,chromsome from contig,microarray sequence
    #BS00068396_51,2AS,CGAAGCGATCCTACTACATTGCGTTCCTTTCCCACTCCCAGGTCCCCCTA[T/C]ATGCAGGATCTTGATTAGTCGTGTGAACAACTGAAATTTGAGCGCCACAA
    def self.parse(reg_str)
      reg_str.chomp!
      snp = SNPSequence.new


      snp.gene, snp.chromosome, snp.sequence_original = reg_str.split(",")
      #snp.position = snp.position.to_i
      #snp.original.upcase!
      #snp.snp.upcase!  
      snp.parse_sequence_snp
      snp.exon_list = Hash.new()
      snp
    end
    
    def parse_snp
      
    end

    def parse_sequence_snp
      pos = 0
      match_data = /(?<pre>\w*)\[(?<org>[ACGT])\/(?<snp>[ACGT])\](?<pos>\w*)/.match(sequence_original)
      @position = Regexp.last_match(:pre).size + 1
      @original = Regexp.last_match(:org)
      @snp = Regexp.last_match(:snp)
      amb_base = Bio::NucleicAcid.to_IUAPC("#{@original}#{@snp}")
      
      @template_sequence = "#{Regexp.last_match(:pre)}#{amb_base}#{Regexp.last_match(:pos)}"
    end
    
    

  end
end