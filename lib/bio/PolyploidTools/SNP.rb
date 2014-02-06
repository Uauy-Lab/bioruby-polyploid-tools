require 'bio'
module Bio::PolyploidTools
  class SNPException < RuntimeError 
   end
  class SNP

    #GENE,ORIGINAL,POS,SNP
    attr_accessor :gene, :original, :position, :snp, :snp_in, :original_name
    attr_accessor :exon_list
    attr_accessor :container
    attr_accessor :flanking_size, :ideal_min, :ideal_max
    attr_accessor :template_sequence
    attr_accessor :use_reference
    attr_accessor :genomes_count

    attr_reader :chromosome

    #Format: 
    #Gene_name,Original,SNP_Pos,pos,chromosome
    #A_comp0_c0_seq1,C,519,A
    def self.parse(reg_str)
      reg_str.chomp!
      snp = SNP.new
      snp.gene, snp.original, snp.position, snp.snp, snp.chromosome = reg_str.split(",")
      snp.position = snp.position.to_i
      snp.original.upcase!
      snp.snp.upcase!  
      snp.exon_list = Hash.new()
      snp.use_reference = false
      snp
    end
    
    def initialize
       @genomes_count = 3 #TODO: if we want to use this with other polyploids, me need to set this as a variable in the main script. 
    end
    
    #We Only want the chromosome, we drop the arm. 
    def chromosome= (chr)
      @chromosome = chr[0,2]
    end
    
    def chromosome_group
      chromosome[0]
    end
    
    def chromosome_genome
      chromosome[1]
    end
    
     def to_fasta
        return ">#{self.gene}\n#{self.template_sequence}\n"
      end

    def add_exon(exon, arm)
      @exon_list[arm] = exon unless @exon_list[arm]
      @exon_list[arm] = exon if exon.record.score > @exon_list[arm].record.score
    end

    def covered_region
      return @covered_region if @covered_region
      if self.use_reference
         reg = Bio::DB::Fasta::Region.new()
         reg.entry = gene
         reg.orientation = :forward
         reg.start = self.position - self.flanking_size
         reg.end = self.position + self.flanking_size
         
         reg.start = 1 if reg.start < 1
         
         return reg
      end
      
      min = @position
      max = @position
     # puts "Calculating covered region for #{self.inspect}"
    #  puts "#{@exon_list.inspect}"
      raise SNPException.new "Exons haven't been loaded for #{self.to_s}" if @exon_list.size == 0
      
      @exon_list.each do | chromosome, exon |
       # puts exon.inspect
        reg = exon.query_region
        min = reg.start if reg.start < min
        max = reg.end if reg.end > max
      end
      reg = Bio::DB::Fasta::Region.new()
      reg.entry = gene
      reg.orientation = :forward
      reg.start = min
      reg.end = max

      @covered_region = reg
      @covered_region
    end

    def left_padding
      flanking_size - self.local_position + 1
      #  primer_region.start - covered_region.start 
      # 0
    end

    def right_padding
      ret =  (2*flanking_size) - (left_padding  + self.covered_region.size ) 
      ret = 0 if ret < 0
      ret  
    end

    def local_position
#      puts "local_position #{self.position} #{self.covered_region.start}"
      self.position - self.covered_region.start
    end

    def padded_position (pos)
      pos + left_padding
    end

    def exon_fasta_string
      gene_region = self.covered_region
      local_pos_in_gene = self.local_position
      ret_str = ""
      container.parents.each  do |name, bam|
        ret_str << ">#{gene_region.entry}-#{self.position}_#{name} Overlapping_exons:#{gene_region.to_s} localSNPpo:#{local_pos_in_gene+1}\n" 
        to_print = parental_sequences[name]
        ret_str << to_print << "\n"
      end
      self.exon_sequences.each do | chromosome, exon_seq | 
        ret_str << ">#{chromosome}\n#{exon_seq}\n"
      end
      mask = masked_chromosomal_snps("1BS", flanking_size)
      ret_str << ">Mask\n#{mask}\n"
      ret_str
    end


    def primer_fasta_string
      gene_region = self.covered_region
      local_pos_in_gene = self.local_position
      ret_str = ""

      surrounding_parental_sequences.each do |name, seq|
        ret_str << ">#{gene_region.entry}-#{self.position}_#{name}\n" 
        ret_str << "#{seq}\n"
      end

      #self.exon_sequences.each do | chromosome, exon_seq | 
      #  ex_seq = cut_sequence_to_primer_region(exon_seq)
      #  ret_str << ">#{chromosome}\n#{ex_seq}\n"
      #end
      self.surrounding_exon_sequences.each do |chromosome, exon_seq|
        ret_str << ">#{chromosome}\n#{exon_seq}\n"
      end

      mask = surrounding_masked_chromosomal_snps(chromosome)
      ret_str << ">Mask\n#{mask}\n"

      pr = primer_region(chromosome, snp_in )
      ret_str << pr.to_fasta
      ret_str
    end

    def primer_region(target_chromosome, parental )
      parental = aligned_sequences[parental].downcase
      chromosome_seq = aligned_sequences[target_chromosome]
      chromosome_seq = "-" * parental.size unless chromosome_seq
      chromosome_seq = chromosome_seq.downcase
      mask = mask_aligned_chromosomal_snp(target_chromosome)
      #puts "'#{mask}'"

      pr = PrimerRegion.new
      position_in_region = 0
      (0..parental.size-1).each do |i|

        if chromosome_seq[i] != '-' or parental[i] != '-'
          case   
          when mask[i] == '&'
            #This is the SNP we take the parental
            pr.snp_pos = position_in_region
            pr.homeologous = false
          when mask[i] == ':'
            #This is the SNP we take the parental
            pr.snp_pos = position_in_region
            pr.homeologous = true
          when mask[i] == '-'
            #When the mask doesnt detect a SNP, so we take the parental
            parental[i] = chromosome_seq[i] unless Bio::NucleicAcid::is_unambiguous(parental[i])

          when /[[:upper:]]/.match(mask[i])
            #This is a good candidate for marking a SNP
            #We validate that the consensus from the sam file accepts the variation from the chromosomal sequence
            if parental[i] == '-'
              parental[i] = mask[i]
              pr.crhomosome_specific_intron << position_in_region
            elsif Bio::NucleicAcid.is_valid(parental[i], mask[i])
              parental[i] = mask[i]
              pr.chromosome_specific << position_in_region
            end
          when /[[:lower:]]/.match(mask[i])
            #this is not that good candidate, but sitll gives specificity

            if parental[i] == '-'
              parental[i] = mask[i]
              pr.almost_crhomosome_specific_intron << position_in_region
            elsif Bio::NucleicAcid.is_valid(parental[i], mask[i])
              parental[i] = mask[i].upcase
              pr.almost_chromosome_specific << position_in_region
            end
          end #Case closes
          position_in_region += 1
        end #Closes region with bases 
      end         

      pr.sequence=parental.gsub('-','')
      pr
    end

    def reverse_complement_string(sequenc_str)
      complement = sequenc_str.tr('atgcrymkdhvbswnATGCRYMKDHVBSWN', 'tacgyrkmhdbvswnTACGYRKMHDBVSWN')
      complement.reverse!
    end

    def return_primer_3_string(opts={})

      left = opts[:left_pos]
      right = opts[:right_pos]
      sequence =  opts[:sequence]
      orientation = "forward"
      if opts[:right_pos]
        orientation = "forward"
        if left > right
          left = sequence.size - left - 1
          right = sequence.size - right - 1
          sequence = reverse_complement_string(sequence)
          orientation = "reverse"
        end
      end

      str = "SEQUENCE_ID=#{opts[:name]} #{orientation}\n"
      str << "SEQUENCE_FORCE_LEFT_END=#{left}\n"
      str << "SEQUENCE_FORCE_RIGHT_END=#{right}\n" if opts[:right_pos]
      str << "SEQUENCE_TEMPLATE=#{sequence}\n"
      str << "=\n"


      #In case that we don't have a right primer, wi do both orientation
      unless opts[:right_pos]
        sequence =  opts[:sequence]    
        left = sequence.size - left - 1
        orientation = "reverse"
        sequence = reverse_complement_string(sequence)
        str << "SEQUENCE_ID=#{opts[:name]} #{orientation}\n"
        str << "SEQUENCE_FORCE_LEFT_END=#{left}\n"
        str << "SEQUENCE_TEMPLATE=#{sequence}\n"
        str << "=\n"
      else

      end

      str
    end


    def primer_3_all_strings(target_chromosome, parental) 
      pr = primer_region(target_chromosome, parental )
      primer_3_propertes = Array.new

      seq_original = String.new(pr.sequence)
      seq_original[pr.snp_pos] = self.original
      seq_original_reverse = reverse_complement_string(seq_original)

      seq_snp =  String.new(pr.sequence)
      seq_snp[pr.snp_pos] =  self.snp
      seq_snp_reverse = reverse_complement_string(seq_snp)

      rev_pos = seq_snp.size - position

      if pr.homeologous
        snp_type = "homeologous"
      else
        snp_type = "non-homeologous"
      end

      pr.chromosome_specific.each do |pos|

        args = {:name =>"#{gene}:#{original}#{position}#{snp} #{original_name} chromosome_specific exon #{snp_type} #{chromosome}", :left_pos => pr.snp_pos, :right_pos => pos, :sequence=>seq_original}
        primer_3_propertes << return_primer_3_string(args)
        args[:name] = "#{gene}:#{original}#{position}#{snp} #{snp_in} chromosome_specific exon #{snp_type} #{chromosome}"
        args[:sequence] = seq_snp
        primer_3_propertes << return_primer_3_string(args)
      end

      pr.almost_chromosome_specific.each do |pos|
        args = {:name =>"#{gene}:#{original}#{position}#{snp} #{original_name} chromosome_semispecific exon #{snp_type} #{chromosome}", :left_pos => pr.snp_pos, :right_pos => pos, :sequence=>seq_original}
        primer_3_propertes << return_primer_3_string(args)
        args[:name] = "#{gene}:#{original}#{position}#{snp} #{snp_in} chromosome_semispecific exon #{snp_type} #{chromosome}"
        args[:sequence] = seq_snp
        primer_3_propertes << return_primer_3_string(args)

      end

      pr.crhomosome_specific_intron.each do |pos|

        args = {:name =>"#{gene}:#{original}#{position}#{snp} #{original_name} chromosome_specific intron #{snp_type} #{chromosome}", :left_pos => pr.snp_pos, :right_pos => pos, :sequence=>seq_original}
        primer_3_propertes << return_primer_3_string(args)
        args[:name] = "#{gene}:#{original}#{position}#{snp} #{snp_in} chromosome_specific exon #{snp_type} #{chromosome}"
        args[:sequence] = seq_snp
        primer_3_propertes << return_primer_3_string(args)
      end

      pr.almost_crhomosome_specific_intron.each do |pos|
        args = {:name =>"#{gene}:#{original}#{position}#{snp} #{original_name} chromosome_semispecific intron #{snp_type} #{chromosome}", :left_pos => pr.snp_pos, :right_pos => pos, :sequence=>seq_original}
        primer_3_propertes << return_primer_3_string(args)
        args[:name] = "#{gene}:#{original}#{position}#{snp} #{snp_in} chromosome_semispecific exon #{snp_type} #{chromosome}"
        args[:sequence] = seq_snp
        primer_3_propertes << return_primer_3_string(args)

      end


      args = {:name =>"#{gene}:#{original}#{position}#{snp} #{original_name} chromosome_nonspecific all #{snp_type} #{chromosome}", :left_pos => pr.snp_pos, :sequence=>seq_original}
      primer_3_propertes << return_primer_3_string(args)
      args[:name] = "#{gene}:#{original}#{position}#{snp} #{snp_in} chromosome_nonspecific all #{snp_type} #{chromosome}"
      args[:sequence] = seq_snp
      primer_3_propertes << return_primer_3_string(args)


      primer_3_propertes
    end

    def to_s
      "#{gene}:#{original}#{position}#{snp}#{chromosome}"
    end

    def primer_3_string(target_chromosome, parental) 
      strings = primer_3_all_strings(target_chromosome, parental) 
      strings.join
    end

    def exon_for_chromosome (chromosome)
      selected_exon=exon_list[chromosome]
      puts "No exon with chromosome #{chromosome} for #{gene}"  unless selected_exon
      selected_exon
    end

    def parental_sequences
      return @parental_sequences if @parental_sequences
      gene_region = self.covered_region
      local_pos_in_gene = self.local_position

      @parental_sequences = Bio::Alignment::SequenceHash.new
      container.parents.each  do |name, bam|
        seq = nil
        if bam
          seq = bam.consensus_with_ambiguities({:region=>gene_region}).to_s
        else
          seq = container.gene_model_sequence(gene_region)
          unless name == self.snp_in
            seq[local_pos_in_gene] = self.original 
          end
        end
        seq[local_pos_in_gene] = seq[local_pos_in_gene].upcase
        
        seq[local_pos_in_gene] = self.snp if name == self.snp_in    
        @parental_sequences [name] = seq
        #puts name
        #puts seq
      end
      @parental_sequences
    end

    def surrounding_parental_sequences
      return @surrounding_parental_sequences if @surrounding_parental_sequences
      gene_region = self.covered_region
      local_pos_in_gene = self.local_position

      @surrounding_parental_sequences = Bio::Alignment::SequenceHash.new
      container.parents.each  do |name, bam|
        seq = nil
        if bam
          seq =  bam.consensus_with_ambiguities({:region=>gene_region}).to_s
        else
          seq = container.gene_model_sequence(gene_region)
           unless name == self.snp_in
             # puts "Modiging original: #{name} #{self.original}"  
              seq[local_pos_in_gene] = self.original 
            end
        end
        #puts "local_pos_in_gene #{local_pos_in_gene}"
        #puts "'#{name}' compared to '#{self.snp_in}'"
        #puts seq
        seq[local_pos_in_gene] = seq[local_pos_in_gene].upcase
        seq[local_pos_in_gene] = self.snp if name == self.snp_in  
        #puts seq
        #puts "__"
        @surrounding_parental_sequences [name] = cut_and_pad_sequence_to_primer_region(seq)
      end
     # puts "&&&&\n#{surrounding_parental_sequences['A']}\n#{surrounding_parental_sequences['B']}\n&&&&"
      @surrounding_parental_sequences
    end

    def cut_sequence_to_primer_region(sequence)
      ideal_min = self.local_position - flanking_size 
      ideal_max = self.local_position + flanking_size
      ideal_min = 0 if ideal_min < 0
      ideal_max = sequence.size - 1 if ideal_max > sequence.size
      # len = ideal_max - ideal_min
      sequence[ideal_min..ideal_max]
    end

    def cut_and_pad_sequence_to_primer_region(sequence)
     # p "cut_and_pad_sequence_to_primer_region #{local_position} #{flanking_size}" 
      ideal_min = self.local_position - flanking_size 
      ideal_max = self.local_position + flanking_size
      left_pad = 0
      right_pad=0
      if ideal_min < 0
        left_pad = ideal_min * -1
        ideal_min = 0 
      end
      if ideal_max > sequence.size
        right_pad = ideal_max - sequence.size
        ideal_max = sequence.size - 1 
      end  
      ret = "-" * left_pad  << sequence[ideal_min..ideal_max] <<  "-" * right_pad 
      ret
    end

    def sequences_to_align
      @sequences_to_align = surrounding_parental_sequences.merge(surrounding_exon_sequences) unless @sequences_to_align
    #  p "sequences_to_align"
     # p @sequences_to_align.inspect
      @sequences_to_align
    end

    def aligned_sequences
     
      return @aligned_sequences if @aligned_sequences
      options = ['--maxiterate', '1000', '--localpair', '--quiet']
      mafft = Bio::MAFFT.new( "mafft" , options)
      #puts "Before MAFT:#{sequences_to_align.inspect}"
      report = mafft.query_align(sequences_to_align)
      @aligned_sequences = report.alignment
      #puts "MAFFT: #{report.alignment.inspect}" 
      @aligned_sequences
    end

    def aligned_sequences_fasta
      ret_str = ""
      aligned_sequences.each_pair do |name, seq|
        ret_str << ">#{self.to_s}-#{name}\n#{seq}\n" 
      end
      ret_str << ">MASK #{chromosome}\n#{mask_aligned_chromosomal_snp(chromosome)}\n"

      pr = primer_region(chromosome, snp_in )
      ret_str << pr.to_fasta
      ret_str
      ret_str 
    end

    def aligned_snp_position
      return @aligned_snp_position if @aligned_snp_position
      pos = -1
      parental_strings = Array.new
      parental_sequences.keys.each do | par |
        
        parental_strings << aligned_sequences[par]
      end
      template_sequence = nil
      aligned_sequences.keys.each do |temp |
        template_sequence = aligned_sequences[ temp ] if  aligned_sequences[ temp ][0] != "-"
      end
      $stderr.puts "WARN: #{self.to_s} #{parental_sequences.keys} is not of size 2 (#{parental_strings.size})" if parental_strings.size != 2

      i = 0
      differences = 0
      local_pos_in_gene = flanking_size
      local_pos = 0
      started = false
#TODO: Validate the cases when the alignment has padding on the left on all the chromosomes

      while i < parental_strings[0].size  do
        if local_pos_in_gene == local_pos
          pos = i
          if parental_strings[0][i] == parental_strings[1][i]
            $stderr.puts "WARN: #{self.to_s} doesn't have a SNP in the marked place (#{i})! \n#{parental_strings[0]}\n#{parental_strings[1]}"
          end 
    
        end

        started = true if template_sequence[i] != "-" 
        if started == false or template_sequence[i] != "-" 
          local_pos += 1
        end
        i += 1
      end
      @aligned_snp_position = pos
      return pos
    end

    def mask_aligned_chromosomal_snp(chromosome)
      names = exon_sequences.keys
      parentals =  parental_sequences.keys

      local_pos_in_gene = aligned_snp_position
      masked_snps = aligned_sequences[chromosome].downcase if aligned_sequences[chromosome]
      masked_snps = "-" * aligned_sequences.values[0].size  unless aligned_sequences[chromosome]
      #TODO: Make this chromosome specific, even when we have more than one alignment going to the region we want.
      i = 0
      while i < masked_snps.size
        different = 0
        cov = 0
        from_group = 0
        names.each do | chr |
          if aligned_sequences[chr] and aligned_sequences[chr][i]  != "-"
            cov += 1 

            from_group += 1 if chr[0] == chromosome_group
            #puts "Comparing #{chromosome_group} and #{chr[0]} as chromosomes"
            if chr != chromosome 
              $stderr.puts "WARN: No base for #{masked_snps} : ##{i}" unless masked_snps[i].upcase
              $stderr.puts "WARN: No base for #{aligned_sequences[chr]} : ##{i}" unless masked_snps[i].upcase
              different += 1  if masked_snps[i].upcase != aligned_sequences[chr][i].upcase 
            end
          end
        end
        masked_snps[i] = "-" if different == 0
        masked_snps[i] = "-" if cov == 1
        masked_snps[i] = "*" if cov == 0
        expected_snps = names.size - 1
       # puts "Diferences: #{different} to expected: #{ expected_snps } [#{i}] Genome count (#{from_group} == #{genomes_count})"
        
        masked_snps[i] = masked_snps[i].upcase if different == expected_snps and from_group == genomes_count

        if i == local_pos_in_gene
          masked_snps[i] = "&"
          bases = ""
          names.each do | chr |
            bases << aligned_sequences[chr][i]  if aligned_sequences[chr] and aligned_sequences[chr][i]  != "-"
          end
          
          code_reference = "n"
          code_reference = Bio::NucleicAcid.to_IUAPC(bases) unless bases == ""

          if Bio::NucleicAcid.is_valid(code_reference,   original) and Bio::NucleicAcid.is_valid(code_reference,   snp)
            masked_snps[i] = ":"
          end 

        end
        i += 1
      end
      masked_snps
    end

    def masked_chromosomal_snps(chromosome)
      chromosomes = exon_sequences
      names = chromosomes.keys
      masked_snps = chromosomes[chromosome].tr("-","+") if chromosomes[chromosome]
      masked_snps = "-" * covered_region.size unless chromosomes[chromosome]
      local_pos_in_gene = self.local_position
      ideal_min = local_pos_in_gene - flanking_size
      ideal_max = local_pos_in_gene + flanking_size
      i = 0
      while i < masked_snps.size  do
        if i > ideal_min and i <= ideal_max

          different = 0
          cov = 0
          names.each do | chr |
            if chromosomes[chr][i]  != "-"
              cov += 1 
              if chr != chromosome and masked_snps[i] != "+"
                different += 1  if masked_snps[i] != chromosomes[chr][i] 
              end
            end

          end
          masked_snps[i] = "-" if different == 0 and  masked_snps[i] != "+"
          masked_snps[i] = "-" if cov < 2
          masked_snps[i] = masked_snps[i].upcase if different > 1   

        else
          masked_snps[i] = "*"
        end   
        if i == local_pos_in_gene
          masked_snps[i] = "&"
        end
        i += 1
      end
      masked_snps
    end

    def surrounding_masked_chromosomal_snps(chromosome)

      chromosomes = surrounding_exon_sequences
      names = chromosomes.keys
      masked_snps = chromosomes[chromosome].tr("-","+") if chromosomes[chromosome]
      masked_snps = "-" * (flanking_size * 2 ) unless chromosomes[chromosome]
      local_pos_in_gene = flanking_size 
      # ideal_min = local_pos_in_gene - flanking_size
      #ideal_max = local_pos_in_gene + flanking_size
      i = 0
      while i < masked_snps.size  do


        different = 0
        cov = 0
        names.each do | chr |
          if chromosomes[chr][i]  != "-" and chromosomes[chr][i]. != 'N' and chromosomes[chr][i]. != 'n'
            cov += 1 
            if chr != chromosome and masked_snps[i] != "+"
              different += 1  if masked_snps[i] != chromosomes[chr][i] 
            end
          end

        end
        masked_snps[i] = "-" if different == 0 and  masked_snps[i] != "+"
        masked_snps[i] = "-" if cov < 2
        masked_snps[i] = masked_snps[i].upcase if different > 1   


        if i == local_pos_in_gene
          masked_snps[i] = "&"
        end
        i += 1
      end
      masked_snps
    end

    def surrounding_exon_sequences
      return @surrounding_exon_sequences if @surrounding_exon_sequences
      @surrounding_exon_sequences =  Bio::Alignment::SequenceHash.new
      self.exon_list.each do |chromosome, exon| 
        #puts "surrounding_exon_sequences #{flanking_size}"
       #puts chromosome
        #puts exon
        flanquing_region  = exon.target_flanking_region_from_position(position,flanking_size)
        #TODO: Padd when the exon goes over the regions... 
        
        #Ignoring when the exon is in a gap
        unless exon.snp_in_gap 
          exon_seq = container.chromosome_sequence(flanquing_region)
          @surrounding_exon_sequences[chromosome] = exon_seq
        end
      end
      @surrounding_exon_sequences
    end  


    def exon_sequences
      return @exon_sequences if @exon_sequences
      gene_region = self.covered_region
      local_pos_in_gene = self.local_position
      @exon_sequences = Bio::Alignment::SequenceHash.new
      self.exon_list.each do |chromosome, exon| 
        exon_start_offset = exon.query_region.start - gene_region.start
        exon_seq  = "-" * exon_start_offset 
        exon_seq << container.chromosome_sequence(exon.target_region).to_s
        #puts exon_seq
       # l_pos = exon_start_offset + local_pos_in_gene
        unless exon.snp_in_gap
          #puts "local position: #{local_pos_in_gene}"
          #puts "Exon_seq: #{exon_seq}"
          exon_seq[local_pos_in_gene] = exon_seq[local_pos_in_gene].upcase
          exon_seq << "-" * (gene_region.size - exon_seq.size + 1)
          @exon_sequences[chromosome] = exon_seq
        end
      end
      @exon_sequences[@chromosome] = "-" * gene_region.size unless @exon_sequences[@chromosome]
      @exon_sequences
    end
  end
end