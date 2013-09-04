module Bio::PolyploidTools
  
  class ChromosomeArm
    attr_accessor :name
    attr_reader :genes
    attr_reader :fasta_db

    def initialize(name, path_to_fasta)
      @name = name
      @fasta_db = Bio::DB::Fasta::FastaFile.new(path_to_fasta)
      @genes = Hash.new
    end
  end

end