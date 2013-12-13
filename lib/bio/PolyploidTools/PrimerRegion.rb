module Bio::PolyploidTools
  class PrimerRegion
    attr_accessor :snp_pos, :sequence, :chromosome_specific, :almost_chromosome_specific, :crhomosome_specific_intron , :almost_crhomosome_specific_intron, :homeologous

    def initialize

      @chromosome_specific = Array.new
      @almost_chromosome_specific = Array.new
      @crhomosome_specific_intron  = Array.new
      @almost_crhomosome_specific_intron = Array.new
    end

    def tail_candidates
      @chromosome_specific.size + @almost_chromosome_specific.size
    end

    def to_fasta
      ">Primer_#{snp_pos}_#{chromosome_specific.to_s}_#{almost_chromosome_specific.to_s}_#{crhomosome_specific_intron.to_s}_#{almost_crhomosome_specific_intron.to_s}\n#{sequence}\n"
    end

  end
end