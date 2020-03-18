
require_relative "SNP"
require 'bio-samtools'
module Bio::PolyploidTools
  class SNPSequenceException < RuntimeError 
  end

  class NoSNPSequence < SNP

    attr_accessor :sequence_original
    #Format: 
    #snp name,chromsome from contig,microarray sequence
    #BS00068396_51,2AS,CGAAGCGATCCTACTACATTGCGTTCCTTTCCCACTCCCAGGTCCCCCTA[T/C]ATGCAGGATCTTGATTAGTCGTGTGAACAACTGAAATTTGAGCGCCACAA
    def self.parse(reg_str)
      reg_str.chomp!
      snp = NoSNPSequence.new

      arr = reg_str.split(",")
      
      if arr.size == 3
        snp.gene, snp.chromosome, snp.sequence_original = reg_str.split(",")
      elsif arr.size == 2
       snp.gene, snp.sequence_original = arr
     else
       throw SNPSequenceException.new "Need two or three fields to parse, and got #{arr.size} in #{reg_str}"
      end
      #snp.position = snp.position.to_i
      #snp.original.upcase!
      #snp.snp.upcase!  
      snp.chromosome. strip!
      snp.snp_in = snp.chromosome
      snp.parse_sequence_snp
      snp.exon_list = Hash.new()
      snp
    end
    
    def parse_snp
       
    end

    def parse_sequence_snp
       @position = (sequence_original.length / 2).to_i 
       @original = sequence_original[@position]
       @snp = @original
    end

    def to_s
      "#{gene}:#{chromosome}"
    end

    def sequences_to_align
      @sequences_to_align = surrounding_exon_sequences unless @sequences_to_align
      @sequences_to_align
    end

     def mask_aligned_chromosomal_snp(chromosome)
      return nil if  aligned_sequences.values.size == 0
      names = aligned_sequences.keys
      parentals =  parental_sequences.keys
      names = names - parentals


      best_target = get_target_sequence(names, chromosome)
      masked_snps = aligned_sequences[best_target].downcase if aligned_sequences[best_target]
      masked_snps = "-" * aligned_sequences.values[0].size  unless aligned_sequences[best_target]

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
       #puts "Diferences: #{different} to expected: #{ expected_snps } [#{i}] Genome count (#{from_group} == #{genomes_count})"
        
        masked_snps[i] = masked_snps[i].upcase if different == expected_snps and from_group == genomes_count

        i += 1
      end
      masked_snps
    end

    def count_deletions_around(position,target_chromosome)
      first_aligned = aligned_sequences[target_chromosome]

      pos_start = position - flanking_size
      pos_end = position + flanking_size
      pos_start = 0 if pos_start < 0
      pos_end = first_aligned.size - 1 if pos_end >= first_aligned.size
      count = 0
      for i in pos_start..pos_end
        has_del = false

        aligned_sequences.each_pair do |name, val|  
          has_del = true if val[i] == '-'
          #print "#{val[i]}\t"
        end
        count += 1 if has_del
        #print "#{count}\n"
      end
      return count
    end

    def primer_region(target_chromosome, parental_chr )
      chromosome_seq = aligned_sequences[target_chromosome]
      names = aligned_sequences.keys
      target_chromosome = get_target_sequence(names, target_chromosome)
      chromosome_seq = aligned_sequences[target_chromosome]
      chromosome_seq = surrounding_exon_sequences[target_chromosome ]if aligned_sequences.size == 0
      chromosome_seq = "-" * sequence_original.size unless chromosome_seq
      chromosome_seq = chromosome_seq.downcase
      #puts chromosome_seq
      mask = mask_aligned_chromosomal_snp(target_chromosome)
    
      pr = PrimerRegion.new
      pr.homoeologous = false
      position_in_region = 0
      parental = chromosome_seq.clone
      (0..chromosome_seq.size-1).each do |i|

        if chromosome_seq[i] != '-'
          case   
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
              pr.chromosome_specific << position_in_region #if count_deletions_around(1,target_chromosome) < 3
              pr.chromosome_specific_in_mask << i
            end

          when /[[:lower:]]/.match(mask[i])
            #this is not that good candidate, but sitll gives specificity
            if parental[i] == '-'
              parental[i] = mask[i]
              pr.almost_crhomosome_specific_intron << position_in_region
            elsif Bio::NucleicAcid.is_valid(parental[i], mask[i])
              parental[i] = mask[i].upcase
              pr.almost_chromosome_specific << position_in_region
              pr.almost_chromosome_specific_in_mask << i
            end
          end #Case closes
          pr.position_in_mask_from_template[position_in_region] = i
          position_in_region += 1
        end #Closes region with bases 
      end         
      pr.sequence=parental.gsub('-','')
      pr
    end

    def return_primer_3_string(opts={})
      #puts "return_primer_3_string #{opts.inspect}"
      left = opts[:left_pos]
      right = opts[:right_pos]
      sequence =  opts[:sequence].clone
      orientation = "forward"
      if opts[:right_pos]
        orientation = "forward"
        if left > right
          left = sequence.size - left - 1
          right = sequence.size - right - 1
          sequence = reverse_complement_string(sequence)
          orientation = "reverse"
        end
        if @variation_free_region > 0
          check_str = sequence[right+1, @variation_free_region]
          return nil if check_str != check_str.downcase
        end

      end


      str = "SEQUENCE_ID=#{opts[:name]} #{orientation}\n"
      str << "SEQUENCE_FORCE_LEFT_END=#{left}\n"
      str << "SEQUENCE_FORCE_RIGHT_END=#{right}\n" if opts[:right_pos]
      str << "SEQUENCE_TEMPLATE=#{sequence}\n"
      str << "=\n"


      #In case that we don't have a right primer, we do both orientations
      unless opts[:right_pos]
        sequence =  opts[:sequence].clone    
        left = sequence.size - left - 1
        orientation = "reverse"
        sequence = reverse_complement_string(sequence)
        str << "SEQUENCE_ID=#{opts[:name]} #{orientation}\n"
        str << "SEQUENCE_FORCE_LEFT_END=#{left}\n"
        str << "SEQUENCE_TEMPLATE=#{sequence}\n"
        str << "=\n"
      end

      str
    end

    def get_base_in_different_chromosome(position, target_chromosome)

        aligned_sequences.each_pair do |name, val|  
          next if target_chromosome == name
          return val[position]
        end
    end

    def primer_3_all_strings(target_chromosome, parental, max_specific_primers: nil) 
      #puts "primer_3_all_strings: #{target_chromosome} #{parental}"
      pr = primer_region(target_chromosome, parental )
      #puts pr.inspect
      primer_3_propertes = Array.new

      seq_original = String.new(pr.sequence)
      #puts seq_original.size.to_s << "-" << primer_3_min_seq_length.to_s
      return primer_3_propertes if seq_original.size < primer_3_min_seq_length

      if pr.homoeologous
        snp_type = "homoeologous"
      else
        snp_type = "non-homoeologous"
      end

      pr.chromosome_specific.each_with_index do |pos , i|
        seq_snp =  seq_original.clone
        #original_base = seq_snp[pos]
        #puts "___"
        #puts aligned_sequences.keys.inspect
        #puts target_chromosome
        t_chr =  get_target_sequence(aligned_sequences.keys, target_chromosome)
        other_chromosome_base = get_base_in_different_chromosome(pr.chromosome_specific_in_mask[i], t_chr)

        args = {
          :name =>"#{gene} A chromosome_specific exon #{snp_type} #{chromosome}", 
          :left_pos => pos,  
          :sequence=>seq_snp
        }
        
        seq_snp =  seq_original.clone
        primer_3_propertes << return_primer_3_string(args)
        
        args[:name] = "#{gene} B chromosome_specific exon #{snp_type} #{chromosome}"
        seq_snp[pos] =  other_chromosome_base.upcase
        args[:sequence] = seq_snp
        
        
        primer_3_propertes << return_primer_3_string(args)
      end

  
      primer_3_propertes
    end

    def aligned_sequences
     
      return @aligned_sequences if @aligned_sequences
      if sequences_to_align.size <= 1
        @aligned_sequences = sequences_to_align
        return @aligned_sequences
      end
      options = ['--maxiterate', '1000', '--localpair', '--quiet']
      mafft = Bio::MAFFT.new( "mafft" , options)
    #  puts "Before MAFT:#{sequences_to_align.inspect}"
      report = mafft.query_align(sequences_to_align)
      @aligned_sequences = report.alignment
   #   puts "MAFFT: #{report.alignment.inspect}" 
      @aligned_sequences
    end
    


   

  end
end