module Bio::PolyploidTools
  
  class ChromosomeArm
    attr_accessor :name
    attr_reader :genes
    attr_reader :loaded_entries
    attr_reader :fasta_db

    def initialize(name, path_to_fasta)
      @name = name
      @fasta_db = Bio::DB::Fasta::FastaFile.new({:fasta=>path_to_fasta})
      #$stderr.puts "Loading entries for #{name}"
      
      @genes = Hash.new
    end
    
    def fetch_contig(contig_id)
      
      @fasta_db.load_fai_entries unless @loaded_entries
      @loaded_entries = true
      entry = fasta_db.index.region_for_entry(contig_id)
     # puts entry
      @fasta_db.fetch_sequence(entry.get_full_region)
    end

    #Loads all the chromosome arms in a folder. 
    #The current version requires that all the references end with .fa, and start with XXX_*.fa
    #Where XXX is the chromosome name
    def self.load_from_folder(path_to_contigs)
      chromosomeArms = Hash.new
      
      Dir.foreach(path_to_contigs) do |filename |
        if  File.fnmatch("*.fa", filename)
          
          parsed = /^(?<arm>\d\w+)/.match(filename)
          target="#{path_to_contigs}/#{filename}"
          #fasta_file = Bio::DB::Fasta::FastaFile.new(target)
          #fasta_file.load_fai_entries
          arm = ChromosomeArm.new(parsed[:arm], target)
          chromosomeArms[arm.name] = arm
        end
      end
      return chromosomeArms
    end
    
  end

end