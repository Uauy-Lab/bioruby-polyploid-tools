module Bio::PolyploidTools
  class PrimerRegion
    attr_accessor :snp_pos, :almost_chromosome_specific_in_mask 
    attr_accessor :chromosome_specific_in_mask, :sequence 
    attr_accessor :chromosome_specific, :almost_chromosome_specific
    attr_accessor :crhomosome_specific_intron , :almost_crhomosome_specific_intron
    attr_accessor :homoeologous, :position_in_mask_from_template

    def initialize

      @chromosome_specific = Array.new
      @almost_chromosome_specific = Array.new
      @crhomosome_specific_intron  = Array.new
      @almost_crhomosome_specific_intron = Array.new
      #For deletions
      @chromosome_specific_in_mask = Array.new
      @almost_chromosome_specific_in_mask = Array.new
      @position_in_mask_from_template = Hash.new
    end

    def tail_candidates
      @chromosome_specific.size + @almost_chromosome_specific.size
    end

    def to_fasta
      ">Primer_#{snp_pos}_#{chromosome_specific.to_s}_#{almost_chromosome_specific.to_s}_#{crhomosome_specific_intron.to_s}_#{almost_crhomosome_specific_intron.to_s}\n#{sequence}\n"
    end

  end
end