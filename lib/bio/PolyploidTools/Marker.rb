module Bio::PolyploidTools
  class Marker
    include Comparable
    #include Virgola
    attr_reader :template_sequence, :original, :snp
    attr_accessor :best_hit
    attr_accessor :index_90k
    attr_accessor :snp_id
    attr_accessor :snp_name
    attr_accessor :chr
    attr_accessor :coordinates_chr
    attr_accessor :map_order
    attr_accessor :chr_arm
    attr_accessor :distance_cm
    attr_accessor :sequence
    attr_writer  :contig



    #after_map :parse_sequence_snp

    def to_fasta
      ">#{self.snp_name}\n#{self.template_sequence}"
    end

    def contig
      @contig = best_hit.target_id.chomp if best_hit
      @contig
    end

    def to_csv
      "#{index_90k},#{snp_id},#{snp_name},#{chr},#{coordinates_chr},#{map_order},#{chr_arm},#{distance_cm},#{sequence},#{contig}"
    end
    
    def <=>(anOter)
     return 0 if anOter.snp_name == @snp_name 
     return @chr_arm <=> anOter.chr_arm  if anOter.chr_arm != @chr_arm
     return @snp_name  <=> anOter.snp_name if anOter.coordinates_chr == @coordinates_chr
     return @coordinates_chr <=> anOter.coordinates_chr    
    end

    def initialize(line)
      line.chomp!
      @template_sequence = nil
      #INDEX_90K,SNP_ID,SNP_NAME,CHR,COORDINATES_CHR,MAP_ORDER,CHR_ARM,DISTANCE_CM,SEQUENCE
      @index_90k, @snp_id, @snp_name, @chr, @coordinates_chr, @map_order, @chr_arm, @distance_cm, @sequence, @contig = line.split(',')
      parse_sequence_snp
    end

    def self.parse(filename)
      f = File.open(filename, "r").read
      f.each_line do |line|
        m = Marker.new(line)
        yield m if m.template_sequence

      end
    end

    protected
    def parse_sequence_snp
      pos = 0
      @chr.upcase!
      match_data = /(?<pre>\w*)\[(?<org>[ACGT])\/(?<snp>[ACGT])\](?<pos>\w*)/.match(sequence)
      if match_data
        @position = Regexp.last_match(:pre).size + 1
        @original = Regexp.last_match(:org)
        @snp = Regexp.last_match(:snp)
        amb_base = Bio::NucleicAcid.to_IUAPC("#{@original}#{@snp}")
        @template_sequence = "#{Regexp.last_match(:pre)}#{amb_base}#{Regexp.last_match(:pos)}"
        return @template_sequence
      end
      return nil
    end
  end


  #The map hast to come sorted. 
  class ArmMap
    attr_reader :markers , :global_reference, :reference
    attr_accessor :chromosome
    def initialize
      @markers = Hash.new
    end

    def align_markers(output)
      Bio::Blat.align(@reference.fasta_path, @fasta_markers, output) do |hit|
        marker = markers[hit.query_id]
        best = marker.best_hit
        unless marker.best_hit 
          markers[hit.query_id].best_hit = hit
        else
          marker.best_hit = hit if hit.score > marker.best_hit.score
        end
      end
    end

    def print_fasta_contigs_for_markers(contigs_file)
      
      contigs = Set.new
      markers.each do |k, marker|

        if marker.best_hit
          contigs << marker.best_hit.target_id
        end
      end
      
      fasta=File.open(contigs_file, "w")
        contigs.each do |contig_id|
             reg = @reference.index.region_for_entry(contig_id)
           fasta.puts ">#{contig_id}\n#{@reference.fetch_sequence(reg.get_full_region)}"
        end
      fasta.close
    end




    def print_fasta_markers(filename)
      @fasta_markers = filename
      fasta=File.open(filename, "w")

      markers.each do |k, marker|
        fasta.puts marker.to_fasta 
      end
      fasta.close
    end

    def global_reference(reference)
      @global_reference = Bio::DB::Fasta::FastaFile.new({:fasta=>reference})
      @global_reference.load_fai_entries
    end

    def reference(reference)
      @reference = Bio::DB::Fasta::FastaFile.new({:fasta=>reference})
      @reference.load_fai_entries
    end

    def print_fasta_contigs_from_reference(filename)
      if File.exist?(filename) 
        reference(filename)
        return
      end

      #puts "loaded"

      fasta=File.open(filename, "w")

      Bio::FlatFile.auto( @global_reference.fasta_path) do |ff|
        ff.each do |f|
          chr_reg = arm_selection_embl(f.entry_id)
          if chr_reg == chromosome
            fasta.puts f.entry
          end
        end
      end
      fasta.close
      reference(filename)
    end


    def print_map_with_contigs(filename)
      file = File.open(filename, "w")
      markers.values.sort { |x,y| x.map_order <=> y.map_order }.each do | marker |
        file.puts marker.to_csv
      end 
      file.close
    end

    protected
    def arm_selection_embl(contig_name)
      ret = contig_name.split('_')[2][0,2]
      return ret
    end
  end
end