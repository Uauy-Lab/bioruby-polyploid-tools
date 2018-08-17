#puts "Loading ExonCointainer..."
module Bio::PolyploidTools
  class ExonContainer
    attr_reader :parental_1_sam,  :parental_2_sam  
    attr_reader :parental_1_name, :parental_2_name, :gene_models_db
    attr_reader :chromosomes, :snp_map
    attr_reader :parents
    attr_accessor :flanking_size , :primer_3_min_seq_length, :max_hits

    BASES = [:A, :C, :G, :T]
    #Sets the reference file for the gene models

    def initialize
      @parents=Hash.new
      @snp_map = Hash.new 
      @snp_contigs
      @primer_3_min_seq_length = 50
      @max_hits = 10
    end

    def gene_models(path)
      @gene_models_db = Bio::DB::Fasta::FastaFile.new({:fasta=>path})
      @gene_models_db.index
      @gene_models_path = path
    end

    #Returns the sequence for a region in the gene models (exon)
    def gene_model_sequence(region)
      #puts "Region: "
      #puts region.inspect
      target_reg = @gene_models_db.index.region_for_entry(region.entry)
      #puts target_reg.inspect
      region.end = target_reg.length if region.end > target_reg.length
      #entries[region.entry]

      seq=@gene_models_db.fetch_sequence(region)
      #puts "sequence: "
      #This is a patch that we need to fix in biosamtools: 
      #puts seq
      index = seq.index('>')
      if(index )
        index -= 1
        #puts "Index: #{index}"
        seq = seq.slice(0..index) 
      end
      #puts seq
      seq
    end

    #Sets the reference file for the gene models
    def chromosomes(path)
      @chromosomes_db = Bio::DB::Fasta::FastaFile.new({:fasta=>path})
      @chromosomes_path = path
    end

    #Retunrs the sequence for a region in the gene models (exon)
    def chromosome_sequence(region)
      left_pad = 0
      #TODO: Padd if it goes to the right
      if(region.start < 1)
        left_pad = region.start * -1
        left_pad += 1
        region.start = 1
      end
      str = "-" * left_pad << @chromosomes_db.fetch_sequence(region)
      #str << "n" * (region.size - str.size + 1) if region.size > str.size
      str
    end 

   
    def add_chromosome_arm(opts)
      @chromosomes = Hash.new unless @chromosomes
      name = opts[:name]
      path = opts[:reference_path]
      path = opts[:alig_path]
      chromosomes[name] = Bio::DB::Fasta::FastaFile.new({:fasta=>path})
    end
    
    def add_snp(snp)
      #TODO: add to the snp the maximum number of hits? 
      snp.max_hits = self.max_hits
      @snp_map[snp.gene] = Array.new unless   @snp_map[snp.gene] 
      @snp_map[snp.gene] << snp

    end

    def add_snp_file(filename, chromosome, snp_in, original_name)

      File.open(filename) do | f |
        f.each_line do | line |
          snp = SNP.parse(line)
          snp.flanking_size = flanking_size
          if snp.position > 0
            snp.container = self
            snp.chromosome = chromosome
            snp.snp_in = snp_in
            snp.original_name = original_name
            @snp_map[snp.gene] = Array.new unless   @snp_map[snp.gene] 
            @snp_map[snp.gene] << snp   
          end

        end
      end
    end

   

    def fasta_string_for_snp(snp)
      gene_region = snp.covered_region
      local_pos_in_gene = snp.local_position
      ret_str = ""
      @parents.each  do |name, bam|
        ret_str << ">#{gene_region.id}_SNP-#{snp.position}_#{name} Overlapping_exons:#{gene_region.to_s} localSNPpo:#{local_pos_in_gene+1}\n" 
        to_print =  bam.consensus_with_ambiguities({:region=>gene_region}).to_s
        to_print[local_pos_in_gene] = to_print[local_pos_in_gene].upcase
        ret_str << to_print << "\n"
      end

      snp.exon_list.each do | chromosome,  exon |
        target_region = exon.target_region
        exon_start_offset = exon.query_region.start - gene_region.start
        chr_local_pos=local_pos_in_gene + target_region.start + 1
        ret_str  << ">#{chromosome}_SNP-#{chr_local_pos} #{exon.to_s} #{target_region.orientation}\n"
        to_print =  "-" * exon_start_offset 
        chr_seq  =  chromosome_sequence(exon.target_region).to_s
        l_pos    =  exon_start_offset + local_pos_in_gene
        to_print << chr_seq
        to_print[local_pos_in_gene] = to_print[local_pos_in_gene].upcase
        ret_str  << to_print
      end
      ret_str
    end

    def print_fasta_snp_exones (file)
      @missing_exons = Set.new unless @missing_exons
      @snp_map.each do | gene, snp_array|
        snp_array.each do |snp|
          #file.puts snp.primer_fasta_string 
          #puts "In print_fast_np_exones"
          #puts snp.inspect

          begin 
            file.puts snp.aligned_sequences_fasta
          rescue Exception=>e
            @missing_exons << snp.to_s
            $stderr.puts "print_fasta_snp_exones:" + snp.to_s + ":" + e.to_s
            $stderr.puts "Local position: #{snp.local_position}"
            $stderr.puts "Local position: #{snp.parental_sequences.to_s}"
            $stderr.puts e.backtrace
          end
        end
      end
    end

    def print_primer_3_exons (file, target_chromosome , parental )
      added = 0
      
      @snp_map.each do | gene, snp_array|
        snp_array.each do |snp|
          string = ""
          begin 
            primer_3_min_seq_length
            string = snp.primer_3_string( snp.chromosome, parental )
            #TODO: add tan error to the SNP this snp has more than max_hits. Or maybe inside the SNP file. 
            #puts "print_primer_3_exons: #{string.size}"
            if string.size > 0
              file.puts string
              added += 1
            end 
           rescue Exception=>e
              @missing_exons << snp.to_s
             # $stderr.puts ""

              $stderr.puts "print_primer_3_exons: #{e.to_s} : snp.to_s"
              $stderr.puts e.backtrace
            end
        end 
      end
      return added
    end

    def add_alignments(opts=Hash.new) 
      opts = { :min_identity=>90, filter_best:false }.merge!(opts)
      exonerate_filename = opts[:exonerate_file]
      arm_selection = opts[:arm_selection]
      filter_best = opts[:filter_best]

      unless arm_selection
        arm_selection = lambda do | contig_name |
          ret = contig_name[0,3]       
          return ret
        end
      end


      File.open(exonerate_filename) do |f|
        f.each_line do | line |
          record = Bio::DB::Exonerate::Alignment.parse_custom(line)
          if  record and record.identity >= opts[:min_identity]
            snp_array = @snp_map[record.query_id]
            if snp_array != nil
              snp_array.each do |snp|                            
                if snp != nil and snp.position.between?( (record.query_start + 1) , record.query_end)
                  begin
                    exon = record.exon_on_gene_position(snp.position)
                    snp.add_exon(exon, arm_selection.call(record.target_id), filter_best:filter_best)
                  rescue Bio::DB::Exonerate::ExonerateException
                    $stderr.puts "Failed for the range #{record.query_start}-#{record.query_end} for position #{snp.position}"
                  end
                end
              end
            end
          end
        end
      end
      remove_alignments_over_max_hits
    end

    def remove_alignments_over_max_hits
      @snp_map.each_pair do | gene, snp_array| 
        snp_array.each do |snp|
          total_hits = snp.exon_list.size
          snp.hit_count = total_hits
          if total_hits > max_hits
            snp.exon_list = {} 
            snp.repetitive = true
            snp.errors << "The marker is in a repetitive region (#{total_hits} hits to reference)"
          end
        end
      end
    end

    def add_parental(opts=Hash.new) 
      # opts = { :name=>opts[:path]}.merge!(opts)
      sam = nil
      name = opts[:name] ? opts[:name] : "Unknown"
      if opts[:path] 
        path = opts[:path]
        name = opts[:name] ? opts[:name] : path.basename(".bam")
        sam =  Bio::DB::Sam.new({:fasta=>@gene_models_path, :bam=>opts[:path]})
      end
      @parents[name] = sam
    end
  end

end