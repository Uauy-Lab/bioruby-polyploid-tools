require 'rubygems'
#require 'extensions/all'
#require 'bio-samtools'
#require 'bio/db/pileup'
#require 'bio/db/vcf'
require 'pathname'
#require_relative 'BIOExtensions.rb'


require 'bio'
require 'bio-samtools'

require "set"
require 'systemu'
require 'json'
#require 'strmask'


module Bio::BFRTools
  
  
  
  class BFRToolsException < StandardError; end


  class Container

    attr_reader :putative_snps, :processed_regions, :total_length, :parental_1_sam,  :parental_2_sam,  :bulk_1_sam,  :bulk_2_sam
    attr_reader :parental_1_name, :parental_2_name, :bulk_1_name, :bulk_2_name, :reference_db

    BASES = [:A, :C, :G, :T]
    #Sets the reference file
    def reference(path)
      @reference_db = Bio::DB::Fasta::FastaFile.new(path)
      @reference_path = path
    end

    def reference_sequence(region)
      @reference_db.fetch_sequence(region)
    end

    #Sets the sorted BAM file of the first parental
    #It accepts the following arguments
    #:name=>A name for thie parental 1 (optional). If not provided, 
    #:path=>
    def parental_1(opts)
      raise BFRToolsException.new("Missing path for parental 1") if opts[:path] == nil
      path = Pathname.new(opts[:path])
      raise BFRToolsException.new("Unable to open #{path}") unless path.readable? or path.directory?        

      @parental_1_name = opts[:name] ? opts[:name] : path.basename(".bam").to_s
      @parental_1_sam =  Bio::DB::Sam.new({:fasta=>@reference_path, :bam=>path.realpath.to_s})
      @parental_1_path = path

    end 

    #Sets the sorted BAM file of the second parental
    def parental_2(opts)
      raise BFRToolsException.new("Missing path for parental 2") if opts[:path] == nil
      path = Pathname.new(opts[:path])
      raise BFRToolsException.new("Unable to open #{path}") unless path.readable? or path.directory?        

      @parental_2_name = @name = opts[:name] ? opts[:name] : path.basename(".bam").to_s
      @parental_2_sam =  Bio::DB::Sam.new({:fasta=>@reference_path, :bam=>path.realpath.to_s})
      @parental_2_path = path
    end

    #Sets the sorted BAM file of the first bulk
    def bulk_1(opts)
      raise BFRToolsException.new("Missing path for bulk 1") if opts[:path] == nil
      path = Pathname.new(opts[:path])
      raise BFRToolsException.new("Unable to open #{path}") unless path.readable? or path.directory?        

      @bulk_1_name =  opts[:name] ? opts[:name] : path.basename(".bam").to_s
      @bulk_1_sam =  Bio::DB::Sam.new({:fasta=>@reference_path, :bam=>path.realpath.to_s})
      @bulk_1_path = path
    end 

    #Sets the sorted BAM file of the second bulk
    def bulk_2(opts)
      raise BFRToolsException.new("Missing path for bulk 2") if opts[:path] == nil
      path = Pathname.new(opts[:path])
      raise BFRToolsException.new("Unable to open #{path}") unless path.readable? or path.directory?        

      @bulk_2_name =  opts[:name] ? opts[:name] : path.basename(".bam").to_s
      @bulk_2_sam =  Bio::DB::Sam.new({:fasta=>@reference_path, :bam=>path.realpath.to_s})
      @bulk_2_path = path
    end     


    def self.snps_between(seq1, seq2)
      snps=0
      for i in (0..seq1.size)
        snps += 1 if seq1[i] != seq2[i]
      end
      snps
    end
  end

  class BFRLine
    attr_reader :original_base, :variation_base, :position, :bulk_1_ratio, :bulk_2_ratio, :bfr


  end

  class BFRRegion < Bio::DB::Fasta::Region
    BASES = [:A, :C, :G, :T]
    attr_reader :parental_1_sequence, :parental_2_sequence, :bulk_1_sequence, :bulk_2_sequence, :snp_count
    attr_reader :ratios_bulk_1, :ratios_bulk_2, :avg_cov_bulk_1, :avg_cov_bulk_2, :coverages_1, :coverages_2, :bases_bulk_1, :bases_bulk_2

    def initialize(opts)
      opts = { :min_cov=>20, :max_snp_1kbp => 5 }.merge!(opts)
      reg = Bio::DB::Fasta::Region.parse_region(opts[:region])
      self.entry = reg.entry
      self.start = reg.start
      self.end   = reg.end

      @container = opts[:container]

      parental_1_sam = @container.parental_1_sam
      parental_2_sam = @container.parental_2_sam
      bulk_1_sam = @container.bulk_1_sam
      bulk_2_sam = @container.bulk_2_sam

      @parental_1_sequence = parental_1_sam.consensus_with_ambiguities(opts)
      @parental_2_sequence = parental_2_sam.consensus_with_ambiguities(opts)

      @bulk_1_sequence = bulk_1_sam.consensus_with_ambiguities(opts)
      @bulk_2_sequence = bulk_2_sam.consensus_with_ambiguities(opts)

      @snp_count = Container.snps_between( @parental_1_sequence , @parental_2_sequence )

      @ratios_bulk_1 = bulk_1_sam.base_ratios_in_region(opts)
      @ratios_bulk_2 = bulk_2_sam.base_ratios_in_region(opts)

      @bases_bulk_1 = bulk_1_sam.bases_in_region(opts)
      @bases_bulk_2 = bulk_2_sam.bases_in_region(opts)

      @avg_cov_bulk_1 = bulk_1_sam.average_coverage_from_pileup(opts)
      @avg_cov_bulk_2 = bulk_2_sam.average_coverage_from_pileup(opts)

      @coverages_1 =  bulk_1_sam.coverages_from_pileup(opts)
      @coverages_2 =  bulk_2_sam.coverages_from_pileup(opts)

    end

    def get_bfr_lines(opts = {})
     
      opts = { :min_cov=>20, :max_snp_1kbp => 5 }.merge!(opts)
      p opts.inspect
      region = self
      line  = String.new
      info = Array.new
      
      for i in (0..region.size-1)
        
        if region.coverages_1[i] > opts[:min_cov] and region.coverages_2[i] > opts[:min_cov]
          BASES.each do |base|

            info.clear
            if  Bio::NucleicAcid.is_valid( region.parental_1_sequence[i],  base.to_s  ) and 
              not  Bio::NucleicAcid.is_valid( region.parental_2_sequence[i],  base.to_s  )
              info << :first
            end

            if   Bio::NucleicAcid.is_valid( region.parental_2_sequence[i],  base.to_s  ) and 
              not Bio::NucleicAcid.is_valid( region.parental_1_sequence[i],  base.to_s  )
              info << :second
            end


            for informative in info
              l = region.get_bfr_line(i, base, informative)
              puts l << "\n"
              line << l << "\n"

              #     output.print  line , "\n"
            end
          end
        end
      end
      line
    end


    def snp_1kbp
      @snp_count.to_f * 1000 / self.size.to_f
    end

    def bfrs
      return @BFRs if @BFRs
      @BFRs = Hash.new
      
      [:first, :second].each do | reference |
        @BFRs[reference] = Hash.new
        BASES.each do |base|
          @BFRs[reference][base] = Array.new
        end
      end
      

      for i in (0..self.size-1)
        ratios_1 = @ratios_bulk_1[i]
        ratios_2 = @ratios_bulk_2[i]
        BASES.each do |base|
          
          if ratios_1[base] == 0 and ratios_2[base] == 0
            bfr1 = 0
            bfr2  = 0
          elsif ratios_1[base] == 0
           bfr1  = 0
           bfr2 = Float::INFINITY
          elsif ratios_2[base] == 0
           bfr1 = Float::INFINITY
            bfr2 = 0
            #bfr = Float::INFINITY
          else
            bfr1  =  ratios_1[base] / ratios_2[base]
            bfr2 =  ratios_2[base] / ratios_1[base]
          end
          @BFRs[:first][base] << bfr1
          @BFRs[:second][base] << bfr2
        end
      end
      @BFRs
    end

    def get_bfr_line(position, base, reference)
      if(reference == :first)
        informative = @container.parental_1_name
        ref_base = @parental_2_sequence[position]
      elsif(reference == :second )
        informative = @container.parental_2_name
        ref_base = @parental_1_sequence[position]
      else
        raise BFRToolsException.new ("The reference for the line should be :first or :second, but was " + reference.to_s )
      end
      
      relative_position = self.start +  position + 1
      
      bfr = bfrs[reference][base][position]
      cov_1 = @coverages_1[position]
      cov_2 = @coverages_2[position]
      ratios_1 = @ratios_bulk_1[position][base]
      ratios_2 = @ratios_bulk_2[position][base]
      base_1_count = @bases_bulk_1[position][base.to_sym]
      base_2_count = @bases_bulk_2[position][base.to_sym]
      #puts "bases_1 #{@bases_bulk_1[position].to_s}"
      #puts "bases_2 #{@bases_bulk_2[position].to_s}"
      line = String.new
      line << @container.parental_1_name  << "\t" << @container.parental_2_name << "\t" <<  @container.bulk_1_name << "\t" << @container.bulk_2_name << "\t" << self.entry << "\t"
      line << ref_base  << "\t" << relative_position.to_s 
      line << "\t" << base.to_s << "\t" 
      line << bfr.round(2).to_s << "\t"  
      line << cov_1.to_s << "\t" << cov_2.to_s  << "\t" 
      line << informative 
      line << "\t" << ratios_1.round(2).to_s << "\t" << ratios_2.round(2).to_s
      line << "\t" << base_1_count.to_s << "\t" << base_2_count.to_s
      line
    end

    def to_multi_fasta
      fasta_string = String.new
      fasta_string << ">"<< self.to_s << ":" << @container.parental_1_name << "\n" << @parental_1_sequence << "\n"
      fasta_string << ">"<< self.to_s << ":" << @container.parental_2_name << "\n" << @parental_2_sequence << "\n"
      fasta_string << ">"<< self.to_s << ":" << @container.bulk_1_name << "\n" << @bulk_1_sequence << "\n"
      fasta_string << ">"<< self.to_s << ":" << @container.bulk_2_name << "\n" << @bulk_2_sequence << "\n"
      fasta_string 
    end

    def to_json (opts)
      #      puts JSON.dump self
      #      JSON.dump self
      #{}"{\"firstName\": \"John\"}"
      out = String.new
      out << "{" 
      out << "\"Parental_1\" : \"" << @container.parental_1_name << "\"\n"
      out << "\"Parental 2\" : \"" << @container.parental_2_name << "\"\n"
      out << "\"Bulk 1\" :  \"" <<  @container.bulk_1_name << "\"\n"
      out << "\"Bulk 2\" : \"" << @container.bulk_2_name << "\"\n"
      out << "\"Positions\" : " << (1..self.size).to_a.to_json << "\n" #TODO: Make this for any subsection, so we can subquery in case we are working on something bigger
      out << "\"Parental_1_consensus\":" << @parental_1_sequence .split(//).to_json << "\n"
      out << "\"Parental_2_consensus\":" << @parental_2_sequence .split(//).to_json << "\n"
      out << "\"Bulk_1_consensus\":" << @bulk_1_sequence .split(//).to_json << "\n"
      out << "\"Bulk_1_coverage\":" << @coverages_1.to_json << "\n"
      #  puts BASES
      
      BASES.each do |base|
        out << "\"Bases_Bulk_1"  << base.to_s << "\":" <<  base_count_for_base(base, @bases_bulk_1).join(",") << "\n"
        out << "\"Ratios_Bulk_1" << base.to_s << "\":" << base_ratios_for_base(base, @ratios_bulk_1).join(",")  << "\n"
      end
      out << "\"Bulk_2_consensus\":" << @bulk_2_sequence .split(//).join(",") << "\n"
      out << "\"Bulk_2_coverage\":" << @coverages_2.join(",") << "\n"

      BASES.each do |base|
        out << "\"Bases_Bulk_2"<< base.to_s << "\":" <<  base_count_for_base(base, @bases_bulk_2).join(",") << "\n"
        out << "\"Ratios_Bulk_2" << base.to_s << "\":" << base_ratios_for_base(base, @ratios_bulk_2).join(",")  << "\n"
      end
      BASES.each do |base|
        out << "\"BFR" << base.to_s << "\":"  << bfrs[:first][base].join(",") << "\n"
      end
      #     << "\t" << @container.bulk_2_name << "\t" << self.entry << "\t"
      out << "}"
      out

    end

    def to_csv
      out = String.new
      out << "Parental 1," << @container.parental_1_name << "\n"
      out << "Parental 2," << @container.parental_2_name << "\n"
      out << "Bulk 1, " <<  @container.bulk_1_name << "\n"
      out << "Bulk 2," << @container.bulk_2_name << "\n"
      out << "Positions," << (1..self.size).to_a.join(",") << "\n"
      out << "Parental 1 consensus," << @parental_1_sequence .split(//).join(",") << "\n"
      out << "Parental 2 consensus," << @parental_2_sequence .split(//).join(",") << "\n"
      out << "Bulk 1 consensus," << @bulk_1_sequence .split(//).join(",") << "\n"
      out << "Bulk 1 coverage," << @coverages_1.join(",") << "\n"
      #  puts BASES
      BASES.each do |base|
        out << "Bases Bulk 1"<< base.to_s << "," <<  base_count_for_base(base, @bases_bulk_1).join(",") << "\n"
        out << "Ratios Bulk 1 " << base.to_s << "," << base_ratios_for_base(base, @ratios_bulk_1).join(",")  << "\n"
      end
      out << "Bulk 2 consensus," << @bulk_2_sequence .split(//).join(",") << "\n"
      out << "Bulk 2 coverage," << @coverages_2.join(",") << "\n"

      BASES.each do |base|
        out << "Bases Bulk 2   "<< base.to_s << "," <<  base_count_for_base(base, @bases_bulk_2).join(",") << "\n"
        out << "Ratios Bulk 2 " << base.to_s << "," << base_ratios_for_base(base, @ratios_bulk_2).join(",")  << "\n"
      end
      BASES.each do |base|
        out << "BFRs" << base.to_s << ","  << bfrs[:first][base].join(",") << "\n"
      end
      #     << "\t" << @container.bulk_2_name << "\t" << self.entry << "\t"
      out
    end

    def base_ratios_for_base(base, ratios_matrix)
      ratios = Array.new
      for i in (0..ratios_matrix.size-1)
        ratios << ratios_matrix[i][base]
      end
      ratios
    end

    def base_count_for_base(base, base_matrix)
      bases = Array.new
      for i in (0..base_matrix.size-1)
        bases << base_matrix[i][base]
      end
      bases
    end

  end


  class BFRContainer < Container

    def init_counters
      @putative_snps      = 0
      @proccesed_regions  = 0
      @not_enogh_coverage = 0
      @total_avg_coverage_bulk_1 = 0.0
      @total_avg_coverage_bulk_2 = 0.0
      @total_snp_1kbp = 0.0
      @no_snps = 0
      @too_many_snps = 0

    end
    def print_header(opts={}) 
      output = opts[:output_file_stats] ? opts[:output_file_stats] : $stderr
      output.print "#bulk_1\tbulk_2\tProcessed_regions\tputative_snps\tno_snps\ttoo_many_snps\tno_enough_coverage\tavg_cov_bulk_1\tavg_cov_bulk_2\tavg_snp_1kbp\n"
    end

    def print_stats(opts={}) 
      output = opts[:output_file_stats] ? opts[:output_file_stats] : $stderr
      output.print @bulk_1_name, "\t", @bulk_2_name, "\t"
      output.print @proccesed_regions, "\t", @putative_snps, "\t", @no_snps, "\t", @too_many_snps,"\t", @not_enogh_coverage, "\t"
      output.print @total_avg_coverage_bulk_1/@proccesed_regions, "\t",@total_avg_coverage_bulk_2/@proccesed_regions, "\t" 
      output.print @total_snp_1kbp / @proccesed_regions,"\n"
    end

    def get_region(opts={})
      opts[:container] = self
      region = BFRRegion.new(opts)
    end

    def process_region(opts={})        
      opts = { :min_cov=>20, :max_snp_1kbp => 10 }.merge!(opts)

      @proccesed_regions += 1
      output = opts[:output_file] ? opts[:output_file] : $stdout
      print_output = opts[:output_file] ? true : false
      opts[:container] = self

      region = BFRRegion.new(opts)

      #puts region.to_multi_fasta

      @total_snp_1kbp += region.snp_1kbp 
     # puts "SNPS: #{region.snp_1kbp}"
      if region.snp_count == 0
        @no_snps += 1 
        print_output = false 
      end

      if region.snp_1kbp  > opts[:max_snp_1kbp]
        @too_many_snps += 1
        print_output = false
      end



      @total_avg_coverage_bulk_2 += region.avg_cov_bulk_2
      @total_avg_coverage_bulk_1 += region.avg_cov_bulk_1

      if region.avg_cov_bulk_2 < opts[:min_cov] or region.avg_cov_bulk_1 < opts[:min_cov]
        @not_enogh_coverage += 1
        print_output = false
      end

      info = Array.new

      if print_output
        for i in (0..region.size-1)
          if region.coverages_1[i] > opts[:min_cov] and region.coverages_2[i] > opts[:min_cov]
            BASES.each do |base|

              info.clear
              if  Bio::NucleicAcid.is_valid( region.parental_1_sequence[i],  base.to_s  ) and 
                not  Bio::NucleicAcid.is_valid( region.parental_2_sequence[i],  base.to_s  )
                info << :first
              end

              if   Bio::NucleicAcid.is_valid( region.parental_2_sequence[i],  base.to_s  ) and 
                not Bio::NucleicAcid.is_valid( region.parental_1_sequence[i],  base.to_s  )
                info << :second
              end


              for informative in info
                line = region.get_bfr_line(i+1, base, informative)
                output.print  line , "\n"
              end
            end
          end
        end
      end


      @parental_1_sam.mpileup_clear_cache region
      @parental_2_sam.mpileup_clear_cache region
      @bulk_2_sam.mpileup_clear_cache region
      @bulk_1_sam.mpileup_clear_cache region
      return region
    end

  end


end

