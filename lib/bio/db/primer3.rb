
module Bio::DB::Primer3
  class Primer3Exception < RuntimeError 
  end
  class SNP

    attr_accessor :gene, :original, :position, :snp, :chromosome, :line_1, :line_2
    attr_accessor :primer3_line_1, :primer3_line_2, :template_length
    attr_accessor :primers_line_1, :primers_line_2
    def line_1_name
      "#{gene}:#{position}#{original}>#{snp} #{line_1}}"
    end

    def initialize
      @primers_line_1 = SortedSet.new
      @primers_line_2 = SortedSet.new
    end

    def line_2_name
      "#{gene}:#{position}#{original}>#{snp} #{line_2}}"
    end

    def to_s
      "#{gene}:#{original}#{position}#{snp}"
    end

    def print_primers_old
      count = 0
      count += 1 if @primer3_line_1    
      count += 1 if @primer3_line_2

      str = "#{gene},#{position}#{original}>#{snp},#{template_length},#{count}," 
      if @primer3_line_1 
        str << primer3_line_1.primer_left_0_sequence << "," << primer3_line_1.primer_right_0_sequence << "," << primer3_line_1.primer_pair_0_product_size << "," <<  primer3_line_1.type.to_s << ","
      else
        str << ",,,,"
      end
      if @primer3_line_2 
        str << primer3_line_2.primer_left_0_sequence << "," <<  primer3_line_2.primer_right_0_sequence << "," << primer3_line_2.primer_pair_0_product_size << "," <<  primer3_line_2.type.to_s 
      else
        str << ",,,"
      end  
      str

    end
    
    def find_left_primer_temp(primer)
      
      primers_line_1.each do |pr|
        #$stderr.puts pr
        #$stderr.puts pr.find_left_tm(primer)
      
        return pr.find_left_tm(primer) if pr.find_left_tm(primer)
      end
      primers_line_2.each do |pr|
        # $stderr.puts pr
        #  $stderr.puts pr.find_left_tm(primer)
        return pr.find_left_tm(primer) if pr.find_left_tm(primer)
      end
      return "NA"
    end
    

    def find_primer_pair_first

      primers_line_1.each do |pr|
        primer = pr.left_primer_snp(self)
        #puts pr.left_primer
        #puts primer
        #puts  find_left_primer_temp(primer)
        return pr if find_left_primer_temp(primer) != "NA"
      end
      nil
    end

    def find_primer_pair_second
     
      primers_line_2.each do |pr|
        primer = pr.left_primer_snp(self)
        # puts pr.left_primer
        #  puts primer
        #  puts  find_left_primer_temp(primer)
        return pr if find_left_primer_temp(primer) != "NA"
      end
      nil
    end


    def print_primers
      
      left_start = 0
      left_end = 0
      right_start = 0
      right_end = 0
      values = Array.new
      #values << "#{gene},,#{template_length},"
      values << gene
      values << "#{original}#{position}#{snp}"
      values << template_length

      if primer3_line_1 and primer3_line_2
        values <<  primer3_line_1.polymorphism

        #Block that searches both if both pairs have a TM
        primer_2 = primer3_line_2.left_primer_with_coordinates(primer3_line_1.left_coordinates, primer3_line_1.orientation)
        primer_2_tm = find_left_primer_temp(primer_2)
        primer_1 = primer3_line_1.left_primer_with_coordinates(primer3_line_2.left_coordinates, primer3_line_2.orientation) 
        primer_1_tm = find_left_primer_temp(primer_1)
      #  $stderr.puts primer_1
      #  $stderr.puts primer_2
        if primer3_line_1 < primer3_line_2 and primer_2_tm != "NA"
          values << primer3_line_1.left_primer
          values << primer_2
          values << primer3_line_1.right_primer 
          values << primer3_line_1.type.to_s 
          values << primer3_line_1.orientation.to_s 
          values << primer3_line_1.primer_left_0_tm 
          values << primer_2_tm
          values << primer3_line_1.primer_right_0_tm
          values << "first" 
        elsif  primer_1_tm != "NA"
          values << primer_1
          values << primer3_line_2.left_primer
          values << primer3_line_2.right_primer
          values << primer3_line_2.type.to_s
          values << primer3_line_2.orientation.to_s
          values << primer_1_tm
          values << primer3_line_2.primer_left_0_tm
          values << primer3_line_2.primer_right_0_tm
          values << "second"
        else
          first_candidate = find_primer_pair_first
          second_candidate = find_primer_pair_second
          
          if first_candidate
            primer_2 = primer3_line_2.left_primer_with_coordinates(first_candidate.left_coordinates, first_candidate.orientation)
            primer_2_tm = find_left_primer_temp(primer_2)
          end
          if second_candidate
            primer_1 = primer3_line_1.left_primer_with_coordinates(second_candidate.left_coordinates, second_candidate.orientation) 
            primer_1_tm = find_left_primer_temp(primer_1)
          end
          
          if first_candidate and second_candidate and first_candidate < second_candidate 
            values << first_candidate.left_primer
            values << primer_2
            values << first_candidate.right_primer 
            values << first_candidate.type.to_s 
            values << first_candidate.orientation.to_s 
            values << first_candidate.primer_left_0_tm 
            values << primer_2_tm
            values << first_candidate.primer_right_0_tm
            values << "first" 
          elsif  second_candidate 
            values << primer_1
            values << second_candidate.left_primer
            values << second_candidate.right_primer
            values << second_candidate.type.to_s
            values << second_candidate.orientation.to_s
            values << primer_1_tm
            values << second_candidate.primer_left_0_tm
            values << second_candidate.primer_right_0_tm
            values << "second"
          elsif  first_candidate 
            values << primer_2
            values << first_candidate.left_primer
            values << first_candidate.right_primer
            values << first_candidate.type.to_s
            values << first_candidate.orientation.to_s
            values << primer_2_tm
            values << first_candidate.primer_left_0_tm
            values << first_candidate.primer_right_0_tm
            values << "first"
          else
            values << "Pair not found"
          end
          
        end
        
      elsif primer3_line_1 
        values << primer3_line_1.polymorphism
        values << primer3_line_1.left_primer
        values << primer3_line_1.left_primer_snp(self) 
        values << primer3_line_1.right_primer 
        values << primer3_line_1.type.to_s 
        values << primer3_line_1.orientation.to_s      
        values << primer3_line_1.primer_left_0_tm 
        values << "NA"
        values << primer3_line_1.primer_right_0_tm
       
        values << "first+"

      elsif primer3_line_2 
        values << primer3_line_2.polymorphism
        values << primer3_line_2.left_primer_snp(self) 
        values << primer3_line_2.left_primer
        values << primer3_line_2.right_primer
        values << primer3_line_2.type.to_s
        values << primer3_line_2.orientation.to_s
        values << "NA"
        values << primer3_line_2.primer_left_0_tm
        values << primer3_line_2.primer_right_0_tm
        values << "second+"

      end 
      values.join(",")
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

    #TODO: make this pick the best primer. At the minute, it just prints the first in the list. 
    def add_record(primer3record)
      @template_length = primer3record.sequence_template.size
      
      case
      when primer3record.line == @line_1
        @line_1_template = primer3record.sequence_template
      when primer3record.line == @line_2
        @line_2_template = primer3record.sequence_template
      else
        raise Primer3Exception.new "#{primer3record.line} is not recognized (#{line_1}, #{line_2})"
      end
      
      if primer3record.primer_left_num_returned.to_i > 0
        case
        when primer3record.line == @line_1
          primers_line_1 << primer3record
          @primer3_line_1 = primer3record if not @primer3_line_1  or @primer3_line_1 > primer3record
        when primer3record.line == @line_2
          primers_line_1 << primer3record
          @primer3_line_2 = primer3record if not @primer3_line_2 or @primer3_line_2 > primer3record
        else
          raise Primer3Exception.new "#{primer3record.line} is not recognized (#{line_1}, #{line_2})"
        end
      end
    end
  end

  class Primer3Record
    include Comparable
    attr_accessor :properties, :polymorphism

    def method_missing(method_name, *args)

      return @properties[method_name] if @properties[method_name] 
      $stderr.puts "Missing #{method_name}"
      puts @properties.inspect
      raise NoMethodError 
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

    def <=>(anOther)
      ret = snp <=> anOther.snp
      return ret if ret != 0


      #Sorting by the types. 
      if type == :chromosome_specific 
        if anOther.type != :chromosome_specific
          return -1
        end
      elsif type == :chromosome_semispecific
        if anOther.type == :chromosome_specific
          return 1
        else anOther.type == :chromosome_nonspecific
          return -1
        end
      elsif type == :chromosome_nonspecific
        if anOther.type != :chromosome_nonspecific
          return 1
        end
      end

      #Sorting if it is in intron or not This will give priority 
      #to the cases when we know for sure the sequence from the line
      #and reduce the chances of getting messed with a short indel
      if self.exon?
        unless anOther.exon? 
          return -1
        end
      else
        if anOther.exon?
          return 1
        end
      end

      #Sorting for how long the product is, the shorter, the better 
      return  product_length <=> anOther.product_length

    end

    def parse_coordinates(str)
      coords = str.split(',')
      coords[0] = coords[0].to_i
      coords[1] = coords[1].to_i
      coords
    end
    
    
    def left_coordinates
      @left_coordinates = parse_coordinates(self.primer_left_0) unless @left_coordinates 
      @left_coordinates 
    end
    
    def right_coordinates
      unless @right_coordinates 
        @right_coordinates = parse_coordinates(self.primer_right_0) 
        @right_coordinates[0] = @right_coordinates[0] - @right_coordinates[1] + 1
      end
     @right_coordinates 
    end
    
    def left_primer
      @left_primer = self.sequence_template[left_coordinates[0],left_coordinates[1]] unless @left_primer
      @left_primer
    end
    
    def left_primer_snp(snp)
      tmp_primer = String.new(left_primer)
      if self.orientation == :forward
        base_original = snp.original 
        base_snp = snp.snp
      elsif self.orientation == :reverse
        base_original = reverse_complement_string(snp.original )
        base_snp = reverse_complement_string(snp.snp)
      else
        raise Primer3Exception.new "#{self.orientation} is not a valid orientation"
      end
      
     # puts "#{snp.to_s} #{self.orientation} #{tmp_primer[-1] } #{base_original} #{base_snp}"
      if tmp_primer[-1] == base_original
        tmp_primer[-1] = base_snp
      elsif tmp_primer[-1] == base_snp
        tmp_primer[-1] = base_original  
      else
         raise Primer3Exception.new "#{tmp_primer} doesnt end in a base in the SNP #{snp.to_s}"
      end
      return tmp_primer
    end
    
    def left_primer_with_coordinates(coordinates, other_orientation)
      
      seq = self.sequence_template
      
       seq = reverse_complement_string(seq) if self.orientation != other_orientation
      
      seq[coordinates[0],coordinates[1]] 
    end
    
    def reverse_complement_string(sequenc_str)
      complement = sequenc_str.tr('atgcrymkdhvbswnATGCRYMKDHVBSWN', 'tacgyrkmhdbvswnTACGYRKMHDBVSWN')
      complement.reverse!
    end
    
    def right_primer
      @right_primer = self.sequence_template[right_coordinates[0],right_coordinates[1]] unless @right_primer
      @right_primer = reverse_complement_string(@right_primer)
     # puts self.sequence_template
      #puts right_coordinates.inspect
      #puts @right_primer.inspect
      @right_primer
    end

    def product_length
      #TODO: Pick the shortest? Distance between primers? 
      return self.primer_pair_0_product_size
    end

    def initialize
      @properties = Hash.new
    end

    def snp
      return @snp if @snp
      parse_header
      @snp
    end
    
    #CL3339Contig1:T509C AvocetS chromosome_specific exon forward
    def parse_header
      @snp, @line, @type, @in, @polymorphism, @orientation  = self.sequence_id.split(" ")  
      @type = @type.to_sym
      if @in
        @in = @in.to_sym == :exon 
      else
        @exon = false
      end
      
      if @polymorphism.to_sym == :homeologous
        @homeologous = true
      else
        @homeologous = false
      end
      @parsed = true
      @orientation = @orientation.to_sym
    end
    
    def orientation
      return @orientation if @parsed
      parse_header
      @orientation
    end

    def homeologous?
      return @homeologous if @parsed
      parse_header
      @homeologous
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

    def self.parse_file(filename)
      File.open(filename) do | f |
        record = Primer3Record.new
        f.each_line do | line |
          line.chomp!
          if line == "="
            yield record

            record = Primer3Record.new
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

  class KASPContainer

    attr_accessor :line_1, :line_2
    attr_accessor :snp_hash

    def add_snp_file(filename)
      @snp_hash=Hash.new unless @snp_hash
      SNP.parse_file(filename) do |snp|
        @snp_hash[snp.to_s] = snp
        snp.line_1 = @line_1
        snp.line_2 = @line_2
      end

    end

    def add_primers_file(filename)
      Primer3Record.parse_file(filename) do | primer3record |
        current_snp = @snp_hash[primer3record.snp]
        current_snp.add_record(primer3record)
        #puts current_snp.inspect
      end
    end

    def print_primers
      snp_hash.each do |k, snp|
        puts snp.print_primers
      end

    end

  end

end

