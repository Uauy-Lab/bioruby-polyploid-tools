require 'pp'
module Bio::DB::Primer3
  class Primer3Exception < RuntimeError 
  end

  def self.read_primer_preferences(file, defaults)

    File.open(file) do |f|
      f.each_line do | line | 
        line.chomp!
        arr = line.split("=")
        defaults[arr[0].downcase.to_sym] = arr[1];
      end
    end

    return defaults
  end

  def self.prepare_input_file(file, opts2={})
   opts = {
    :primer_product_size_range => "50-150" ,
    :primer_max_size => 25 , 
    :primer_lib_ambiguity_codes_consensus => 1,
    :primer_liberal_base => 1, 
    :primer_num_return => 5,
    :primer_explain_flag => 1,
    :primer_thermodynamic_parameters_path => File.expand_path(File.dirname(__FILE__) + '../../../../conf/primer3_config/') + '/'
    }.merge(opts2)
    
    opts.each do |key,value|
      file.puts "#{key.to_s.upcase}=#{value}\n"
    end
    # file.puts "="
  end

  def self.run(opts={})
    puts "Primer3.run running..."

    f_in=opts[:in]
    f_out=opts[:out]
    opts.delete(:in)
    opts.delete(:out)
    primer_3_in = File.read(f_in)
    status = systemu "primer3_core", 0=>primer_3_in, 1=>stdout='', 2=>stderr=''
    # $stderr.puts cmdline
    if status.exitstatus == 0
      File.open(f_out, 'w') { |f| f.write(stdout) }
    else
      raise Primer3Exception.new(), "Error running primer3. Command line was 'primer3_core'\nPrimer3 STDERR was:\n#{stderr}"
    end
  end

  class SNP

    attr_accessor :gene, :original, :position, :snp, :chromosome, :line_1, :line_2
    attr_accessor :primer3_line_1, :primer3_line_2, :template_length
    attr_accessor :primers_line_1, :primers_line_2
    attr_accessor :used_contigs
    attr_accessor :snp_from
    attr_accessor :regions
    attr_accessor :primer3_errors
    
    def line_1_name
      "#{gene}:#{position}#{original}>#{snp} #{line_1}}"
    end

    def initialize
      @primers_line_1 = SortedSet.new
      @primers_line_2 = SortedSet.new
      @regions = SortedSet.new
      @primer3_errors = Set.new
    end

    def line_2_name
      "#{gene}:#{position}#{original}>#{snp} #{line_2}}"
    end

    def to_s
      "#{gene}:#{original}#{position}#{snp}:#{snp_from.chromosome}"
    end

    def find_left_primer_temp(primer)
      primers_line_1.each do |pr|
        return pr.find_left_tm(primer) if pr.find_left_tm(primer)
      end
      primers_line_2.each do |pr|
        return pr.find_left_tm(primer) if pr.find_left_tm(primer)
      end
      return "NA"
    end


    def find_primer_pair_first
      primers_line_1.each do |pr|
        primer = pr.left_primer_snp(self)
        return pr if find_left_primer_temp(primer) != "NA"
      end
      nil
    end

    def find_primer_pair_second
      primers_line_2.each do |pr|
        primer = pr.left_primer_snp(self)
        return pr if find_left_primer_temp(primer) != "NA"
      end
      nil
    end

    def values
      return @values if @values
      left_start = 0
      left_end = 0
      right_start = 0
      right_end = 0
      total_columns_before_messages=17
      #puts "Values in primer3"
      #puts snp_from.inspect
      @values = Array.new
      #@values << "#{gene},,#{template_length},"
      @values << gene
      @values << "#{original}#{position}#{snp}"
      @values << template_length
      @values << snp_from.chromosome
      @values << regions.size
      @values << regions.join("|")
      if primer3_line_1 and primer3_line_2
        @values <<  primer3_line_1.polymorphism

        #Block that searches both if both pairs have a TM
        primer_2 = primer3_line_2.left_primer_with_coordinates(primer3_line_1.left_coordinates, primer3_line_1.orientation)
        primer_2_tm = find_left_primer_temp(primer_2)
        primer_1 = primer3_line_1.left_primer_with_coordinates(primer3_line_2.left_coordinates, primer3_line_2.orientation) 
        primer_1_tm = find_left_primer_temp(primer_1)
        #  $stderr.puts primer_1
        #  $stderr.puts primer_2
        if primer3_line_1 < primer3_line_2 and primer_2_tm != "NA"
          @values << primer3_line_1.left_primer
          @values << primer_2
          @values << primer3_line_1.right_primer 
          @values << primer3_line_1.type.to_s 
          @values << primer3_line_1.orientation.to_s 
          @values << primer3_line_1.best_pair.left.tm 
          @values << primer_2_tm
          @values << primer3_line_1.best_pair.right.tm
          @values << "first" 
          @values << primer3_line_1.best_pair.product_size
        elsif  primer_1_tm != "NA"
          @values << primer_1
          @values << primer3_line_2.left_primer
          @values << primer3_line_2.right_primer
          @values << primer3_line_2.type.to_s
          @values << primer3_line_2.orientation.to_s
          @values << primer_1_tm
          @values << primer3_line_2.best_pair.left.tm
          @values << primer3_line_2.best_pair.right.tm
          @values << "second"
          @values << primer3_line_2.best_pair.product_size
        else

          first_candidate = find_primer_pair_first
          second_candidate = find_primer_pair_second

          if first_candidate
            primer_2 = primer3_line_2.left_primer_with_coordinates(first_candidate.left_coordinates, first_candidate.orientation)
            primer_2_tm = find_left_primer_temp(primer_2)
          end
          if second_candidate
            #puts "input to search #{first_candidate.left_coordinates}"
            primer_1 = primer3_line_1.left_primer_with_coordinates(second_candidate.left_coordinates, second_candidate.orientation) 
            primer_1_tm = find_left_primer_temp(primer_1)
            #puts "In the other funky if #{primer_2}"
          end

          if first_candidate and second_candidate and first_candidate < second_candidate 
            #puts "A"
            @values << first_candidate.left_primer
            @values << primer_2
            @values << first_candidate.right_primer 
            @values << first_candidate.type.to_s 
            @values << first_candidate.orientation.to_s 
            @values << first_candidate.best_pair.left.tm 
            @values << primer_2_tm
            @values << first_candidate.best_pair.right.tm
            @values << "first" 
            @values << first_candidate.best_pair.product_size
          elsif  second_candidate 
            #puts "B"
            @values << primer_1
            @values << second_candidate.left_primer
            @values << second_candidate.right_primer
            @values << second_candidate.type.to_s
            @values << second_candidate.orientation.to_s
            @values << primer_1_tm
            @values << second_candidate.best_pair.left.tm
            @values << second_candidate.best_pair.right.tm
            @values << "second"
            @values << second_candidate.best_pair.product_size
          elsif  first_candidate 
            #puts "C"
            @values << first_candidate.left_primer
            @values << primer_2
            @values << first_candidate.right_primer
            @values << first_candidate.type.to_s
            @values << first_candidate.orientation.to_s
            @values << primer_2_tm
            @values << first_candidate.best_pair.left.tm
            @values << first_candidate.best_pair.right.tm
            @values << "first"
            @values << first_candidate.best_pair.product_size
          end
        end

      elsif primer3_line_1 
        @values << primer3_line_1.polymorphism
        @values << primer3_line_1.left_primer
        @values << primer3_line_1.left_primer_snp(self) 
        @values << primer3_line_1.right_primer 
        @values << primer3_line_1.type.to_s 
        @values << primer3_line_1.orientation.to_s      
        @values << primer3_line_1.best_pair.left.tm 
        @values << "NA"
        @values << primer3_line_1.best_pair.right.tm

        @values << "first+"
        @values << primer3_line_1.best_pair.product_size
      elsif primer3_line_2 
        @values << primer3_line_2.polymorphism
        @values << primer3_line_2.left_primer_snp(self) 
        @values << primer3_line_2.left_primer
        @values << primer3_line_2.right_primer
        @values << primer3_line_2.type.to_s
        @values << primer3_line_2.orientation.to_s
        @values << "NA"
        @values << primer3_line_2.best_pair.left.tm
        @values << primer3_line_2.best_pair.right.tm
        @values << "second+"
        @values << primer3_line_2.best_pair.product_size

      end 
      if @values.size < total_columns_before_messages
        @values[total_columns_before_messages] = primer3_errors.to_a.join("|")
      else
        @values << nil
      end
      return @values
    end

    def print_primers
      self.values.join(",")
    end

    def found_primers?
      return self.values[7] && self.values[7] != nil
    end

    def first_primer
      return self.values[7] if self.values[7] && self.values[7] != nil
      return ""
    end

    def second_primer
      return self.values[8] if self.values[8] && self.values[8] != nil
      return ""
    end

    def common_primer 
      return self.values[9] if self.values[9] && self.values[9] != nil
      return ""
    end

    def product_size
      return self.values[16].to_i if self.values[16]&& self.values[16] != nil 
      return 0
    end

    def orientation
      puts "insideOrientation: #{self.values[11]}"
      return self.values[11] if self.values[11]&& self.values[11] != nil 
      return 'unknown'
    end


    def first_product
      left = first_primer
      right = common_primer
      nlen = product_size - left.size - right.size
      product = left + ('n' * nlen) + Primer3Record.reverse_complement_string(right)
      #puts "orientation: #{orientation}"
      
      product = Primer3Record.reverse_complement_string(product)  if orientation == 'reverse'
    
      product
    end

     def second_product
      left = second_primer
      right = common_primer
      nlen = product_size - left.size - right.size
      product = left + ('n' * nlen) + Primer3Record.reverse_complement_string(right)
      product = Primer3Record.reverse_complement_string(product) if orientation == 'reverse'
      

      product
    end


    def realigned_primers_fasta
      ret_str = ""
      realigned_primers.each_pair do |name, seq|
        ret_str << ">#{self.to_s}-#{name}\n#{seq}\n" 
      end
      ret_str 
    end


    def realigned_primers
     
      return @realigned_primers if @realigned_primers
      sequences_to_align = Hash.new
      sequences_to_align["first_product"]  = first_product
      sequences_to_align["second_product"] = second_product
      sequences_to_align.merge!(snp_from.surrounding_exon_sequences)
      if sequences_to_align.size == 1
        @realigned_primers = sequences_to_align
        return @realigned_primers
      end
      options = ['--maxiterate', '1000', '--localpair', '--quiet']
      mafft = Bio::MAFFT.new( "mafft" , options)
      #puts "Before MAFT:#{sequences_to_align.inspect}"
      report = mafft.query_align(sequences_to_align)
      @realigned_primers = report.alignment
      #puts "MAFFT: #{report.alignment.inspect}" 
      @realigned_primers
    end 

    def self.parse(reg_str)
      reg_str.chomp!
      snp = SNP.new
      snp.gene, snp.original, snp.position, snp.snp = reg_str.split(",")
      snp.position = snp.position.to_i
      snp.original.upcase!
      snp.snp.upcase!  
      snp
    end

    def self.parse_file(filename)
      File.open(filename) do | f |
        f.each_line do | line |
          snp = SNP.parse(line)
          if snp.position > 0
            yield snp
          end
        end
      end
    end

    
    def add_record(primer3record)
      @primer3_errors = Set.new unless @primer3_errors
      @template_length = primer3record.sequence_template.size
      if primer3record.primer_error != nil 
        primer3_errors << primer3record.primer_error
        return
      end
      case

      when primer3record.line == @line_1

        @line_1_template = primer3record.sequence_template

      when primer3record.line == @line_2
        @line_2_template = primer3record.sequence_template
      else
        raise Primer3Exception.new "#{primer3record.line} is not recognized (#{line_1}, #{line_2})"
      end

      if  primer3record.primer_left_num_returned.to_i > 0 
        case
        when primer3record.line == @line_1
          primers_line_1 << primer3record
          #puts primer3record.inspect
          @primer3_line_1 = primer3record if not @primer3_line_1  or @primer3_line_1 > primer3record
        when primer3record.line == @line_2
          primers_line_2 << primer3record
          @primer3_line_2 = primer3record if not @primer3_line_2 or @primer3_line_2 > primer3record
        else
          raise Primer3Exception.new "#{primer3record.line} is not recognized (#{line_1}, #{line_2})"
        end
      else
        primer3_errors << "#{primer3record.line}(#{primer3record.orientation}):#{primer3record.primer_left_explain.gsub!(',',';')}"
        primer3_errors << "common(#{primer3record.orientation}#{primer3record.type}):#{primer3record.primer_right_explain.gsub!(',',';')}"
        primer3_errors << "pair(#{primer3record.orientation}#{primer3record.type}):#{primer3record.primer_pair_explain.gsub!(',',';')}"
      end
    end
  end

  class Primer3Record
    include Comparable
    attr_accessor :properties, :polymorphism
    attr_accessor :scores


    def best_pair
      return @best_pair if @best_pair
      @best_pair = nil
      @primerPairs.each do | primer |
        @best_pair = primer if @best_pair.nil?
        @best_pair = primer if primer.size < @best_pair.size
      end
      #@best_pair = @primerPairs.min
      @best_pair
    end

    def primer_error
      return @properties[:primer_error] if @properties[:primer_error]
      return nil
    end
    
    def method_missing(method_name, *args)
      return @properties[method_name] if @properties[method_name] 
      $stderr.puts "Missing #{method_name}"
      $stderr.puts @properties.inspect
      return "" #if a property is missing, return blank. 
      raise NoMethodError.new() 
    end

    def find_left_tm(primer)
      last = size - 1
      (0..last).each do | i |
        seq_prop = "primer_left_#{i}_sequence".to_sym
        #        $stderr.puts seq_prop
        temp_property = "primer_left_#{i}_tm".to_sym  
        #       $stderr.puts "comparing  #{@properties[seq_prop] } == #{primer}"
        return @properties[temp_property]  if @properties[seq_prop] == primer

      end
      return nil
    end

    def score
      ret = 0
      ret += @scores[type]
      ret += @scores[:exon] if exon?
      ret -= product_length
      ret
    end

    def <=>(anOther)
      return  anOther.score <=> score
    end

    def parse_coordinates(str)
      coords = str.split(',')
      coords[0] = coords[0].to_i
      coords[1] = coords[1].to_i
      coords
    end


    def left_coordinates
      #@left_coordinates = parse_coordinates(self.primer_left_0) unless @left_coordinates 
      @left_coordinates = best_pair.left.coordinates
      @left_coordinates 
    end

    def right_coordinates
      unless @right_coordinates 
        @right_coordinates = best_pair.right.coordinates
        @right_coordinates[0] = @right_coordinates[0] - @right_coordinates[1] + 1
      end
      @right_coordinates 
    end

    def left_primer
      #@left_primer = self.sequence_template[left_coordinates[0],left_coordinates[1]] unless @left_primer
      @left_primer = best_pair.left.sequence
      @left_primer
    end

    def left_primer_snp(snp)
      tmp_primer = String.new(left_primer)
      if self.orientation == :forward
        base_original = snp.original 
        base_snp = snp.snp
      elsif self.orientation == :reverse
        #puts self.inspect
        base_original =Primer3Record.reverse_complement_string(snp.original )
        base_snp = Primer3Record.reverse_complement_string(snp.snp)
      else
        raise Primer3Exception.new "#{self.orientation} is not a valid orientation"
      end

      #puts "#{snp.to_s} #{self.orientation} #{tmp_primer[-1] } #{base_original} #{base_snp}"
      if tmp_primer[-1] == base_original
        tmp_primer[-1] = base_snp
      elsif tmp_primer[-1] == base_snp
        tmp_primer[-1] = base_original  
      else
        raise Primer3Exception.new "#{tmp_primer} doesnt end in a base in the SNP #{snp.to_s}"
      end
      #puts "tmp_primer: #{tmp_primer}"
      return tmp_primer
    end

    def left_primer_with_coordinates(coordinates, other_orientation)

      seq = self.sequence_template
      #puts "Left coordinates: #{seq}"
      seq = Primer3Record.reverse_complement_string(seq) if self.orientation != other_orientation

      seq[coordinates[0],coordinates[1]] 
    end

    def self.reverse_complement_string(sequenc_str)
      complement = sequenc_str.tr('atgcrymkdhvbswnATGCRYMKDHVBSWN', 'tacgyrkmhdbvswnTACGYRKMHDBVSWN')
      complement.reverse!
    end

    def right_primer_delete
      @right_primer = self.sequence_template[right_coordinates[0],right_coordinates[1]] unless @right_primer
      @right_primer = Primer3Record.reverse_complement_string(@right_primer)
      @right_primer
    end
    
    def right_primer
      return best_pair.right.sequence
    end

    def product_length
      return best_pair.size
    end

    def initialize
      @properties = Hash.new
      @scores = Hash.new
      @scores[:chromosome_specific] = 1000
      @scores[:chromosome_semispecific] = 100
      @scores[:chromosome_nonspecific] = 0
      @scores[:exon] = 50

    end

    def snp
      return @snp if @snp
      parse_header
      @snp
    end

    #CL3339Contig1:T509C AvocetS chromosome_specific exon 4D forward 
    def parse_header
      #puts "Parsing header: '#{self.sequence_id}'"
      @snp, @line, @type, @in, @polymorphism, @chromosome, @orientation   = self.sequence_id.split(" ")  
      @type = @type.to_sym
      if @in
        @in = @in.to_sym == :exon 
      else
        @exon = false
      end

      if @polymorphism.to_sym == :homoeologous
        @homoeologous = true
      else
        @homoeologous = false
      end
      @parsed = true
      @orientation = @orientation.to_sym
    end

    def orientation
      return @orientation if @parsed
      parse_header
      @orientation
    end

    def chromosome
      return @chromosome if @parsed
      parse_header
      @chromosome
    end
    
    def homoeologous?
      return @homoeologous if @parsed
      parse_header
      @homoeologous
    end

    def type
      return @type if @parsed
      parse_header
      @type
    end

    def exon?
      return @exon if @parsed
      parse_header
      @exon
    end

    def line
      return @line if @parsed
      parse_header
      @line
    end

    def size
      @properties[:primer_pair_num_returned].to_i
    end

    def parse_blocks
      total_blocks = size - 1 
      @primerPairs = Array.new
      for i in 0..total_blocks
        @primerPairs << PrimerPair.new(self, i)
      end

    end

    def self.parse_file(filename, scores: nil)
      File.open(filename) do | f |
        record = Primer3Record.new
        record.scores = scores if scores
        f.each_line do | line |
          line.chomp!
          if line == "="

            record.parse_blocks
            yield record
            record = Primer3Record.new
            record.scores = scores if scores
          else
            tokens = line.split("=")
            i = 0
            reg = ""
            #TODO: Look if there is a join function or something similar to go around this... 
            tokens.each do |tok|
              if i > 0
                if i > 1
                  reg << "="
                end
                reg << tok
              end
              i+=1
            end
            record.properties[tokens[0].downcase.to_sym] = reg
          end
        end
      end
    end
  end


  class Primer
    attr_accessor :pair
    def initialize
      @values = Hash.new
    end

    def method_missing(m, *args, &block)  

      return @values[m.to_s] if @values[m.to_s] != nil
      raise NoMethodError.new(), "There's no method called #{m}, available: #{@values.keys.to_s}."  
    end

    def set_value(key, value)
      @values[key] = value
    end



  end

  class PrimerPair
    include Comparable
    attr_reader :record
    attr_reader :left, :right
    
    def parse_coordinates(str)
      coords = str.split(',')
      coords[0] = coords[0].to_i
      coords[1] = coords[1].to_i
      coords
    end

    def size
      return product_size.to_i
    end

    def <=>(anOther)
      penalty.to_f <=> anOther.penalty.to_f
    end

    def initialize(record, index)
      raise Primer3Exception.new(), "Index #{index} is greater than the number of records" unless index < record.size
      @record = record
      @left = Primer.new
      @right = Primer.new
      @values = Hash.new
      

      @left.set_value("added", false)
      @right.set_value("added", false)
      @left.pair = self
      @right.pair = self
      index_s = index.to_s
      record.properties.each do |key, value|
        tokens = key.to_s.split("_")
        if tokens.size > 2 and tokens[2] == index_s
          primer = nil
          primer = @right if tokens[1] == "right"
          primer = @left if tokens[1] == "left"
          if primer != nil
            primer.set_value("added", true)
            if tokens.size == 3
              primer.set_value("coordinates", parse_coordinates(value) )
            else

              to_add = value
              to_add = value.to_f unless tokens[3]=="sequence"
              n_key = tokens[3..6].join("_")
              primer.set_value(n_key, to_add)
            end
          else
            n_key = tokens[3..6].join("_")
            @values[n_key] = value  
          end

        end
      end

      raise Primer3Exception.new(), "The pair is not complete (l:#{left.added}, r:#{right.added})" if @left.added == false or @right.added == false

    end

    def method_missing(m, *args, &block)  

      return @values[m.to_s] if @values[m.to_s]
      raise NoMethodError.new(), "There's no method called #{m}. Available methods: #{@values.keys.to_s}"
    end
  end

  class KASPContainer

    attr_accessor :line_1, :line_2
    attr_accessor :snp_hash
    attr_accessor :scores

    def add_snp_file(filename)
      @snp_hash=Hash.new unless @snp_hash
      SNP.parse_file(filename) do |snp|
        @snp_hash[snp.to_s] = snp
        snp.line_1 = @line_1
        snp.line_2 = @line_2
      end
    end

    def add_snp(snp_in) 
      #TODO: Here we need to also copy the errors that will be printed. 
      @snp_hash=Hash.new unless @snp_hash
      snp = SNP.new
      snp.gene = snp_in.gene
      snp.original = snp_in.original

      snp.position = snp_in.position
      snp.snp = snp_in.snp

      snp.line_1 = @line_1
      snp.line_2 = @line_2 
      snp.snp_from = snp_in
      snp.regions = snp_in.exon_list.values.collect { |x|  x.collect {|y| y.target_region.to_s }}
      @snp_hash[snp.to_s] = snp
      snp
    end

    def add_primers_file(filename)
      #primer3record.scores = @scores if @scores
      Primer3Record.parse_file(filename, scores: @scores) do | primer3record |
        current_snp = @snp_hash["#{primer3record.snp.to_s}:#{primer3record.chromosome}"]
        current_snp.add_record(primer3record)
      end
    end

    def print_primers
      str = ""
      snp_hash.each do |k, snp|
        str << snp.print_primers << "\n"
      end
      return str
    end

    def print_primers_with_tails(tail_a: "GAAGGTCGGAGTCAACGGATT", tail_b: "GAAGGTGACCAAGTTCATGCT")
      str = ""
      snp_hash.each do |k, snp|  
        if snp.found_primers?
          str << snp.gene << snp.original << "\t" << tail_a << snp.first_primer << "\n"
          str << snp.gene << snp.snp      << "\t" << tail_b << snp.second_primer << "\n"
          str << snp.gene                 << "\t"           << snp.common_primer << "\n"
        end
      end
      return str
    end

  end
end

