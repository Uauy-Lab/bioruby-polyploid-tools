# RYO %S\t%pi\t%ql\t%tl\t%g\t%V\n 


module Bio::DB::Exonerate


  #TODO: Make a proper object with generic parser
  def self.align(opts={})
    opts = {
      :model => 'affine:local' ,
      :ryo => "RESULT:\\t%S\\t%pi\\t%ql\\t%tl\\t%g\\t%V\\n" , 
      :bestn => 20,
      :percentage => 50
    }
    .merge(opts)

    target=opts[:target]
    query=opts[:query]

    cmdline = "exonerate --verbose 0 --showalignment no --bestn #{opts[:bestn]} --showvulgar no  --model #{opts[:model]}   --ryo '#{opts[:ryo]}' #{query} #{target}"
    status, stdout, stderr = systemu cmdline
    #$stderr.puts cmdline
    if status.exitstatus == 0
      alns = Array.new unless block_given?
      stdout.each_line do |line|
        aln = Alignment.parse_custom(line) 
        if aln
          if block_given?
            yield aln
          else
            alns << aln
          end
        end
      end
      return alns unless block_given?
    else
      raise ExonerateException.new(), "Error running exonerate. Command line was '#{cmdline}'\nExonerate STDERR was:\n#{stderr}"
    end
  end
  
  
  class ExonerateException < RuntimeError 
  end

  class Alignment
    attr_accessor   :query_id, :query_start, :query_end, :query_strand
    attr_accessor :target_id, :target_start, :target_end, :target_strand, :score
    attr_accessor :vulgar_block, :pi, :ql, :tl, :g
    attr_accessor :line

    #This one day may grow to work with complex ryo....
    def self.parse_custom(line)
      fields=line.split(/\t/)
      if fields[0] == "RESULT:"
        al = Bio::DB::Exonerate::Alignment.new()      
        al.parse_sugar(fields[1])
        al.pi = fields[2].to_f 
        al.ql = fields[3].to_i
        al.tl = fields[4].to_i
        al.g  = fields[5]
        al.parse_vulgar(fields[6])
        al.line = line
        return al
      else
        return nil
      end
    end

    def query
     unless @query
       @query = Bio::DB::Fasta::Region.new()
       @query.entry = query_id 
       @query.start = query_start + 1
       @query.end = query_end
       @query.orientation = query_strand
       if @query.orientation == :reverse
        @query.end = query_start 
        @query.start = query_end + 1
       end
       @query
     end
     @query
   end

   def target
    unless @target
      @target = Bio::DB::Fasta::Region.new()
      @target.entry = target_id 
      @target.start = target_start + 1
      @target.end = target_end
      @target.orientation = target_strand
      if @target.orientation == :reverse
        @target.end = target_start 
        @target.start = target_end + 1
      end
    end
    @target
  end

  def identity
    @pi
  end
  def query_length
    @ql
  end
  def query_coverage
    total_m = 0
    vulgar_block.each do |v|  
      #p v.label
      if v.label == :M
        total_m += v.query_length
      end
    end
   #puts "Total m #{total_m}"
   #puts "ql #{query_length}"
    return 100.00 * total_m.to_f / query_length.to_f 
  end

  def parse_sugar(sugar_str)
    @query_id, @query_start, @query_end, @query_strand, @target_id, @target_start, @target_end, @target_strand, @score = sugar_str.split(/\s+/)

    @query_start  = @query_start.to_i
    @query_end = @query_end.to_i
    @target_start = @target_start.to_i
    @target_end = @target_end.to_i
    @score = @score.to_f

    if @target_strand == "+"
      @target_strand = :forward
    elsif @target_strand == "-"
      @target_strand = :reverse
    else
      raise  ExonerateException.new(), "Ivalid target orientation #{@target_strand} for line:\n#{sugar_str}"
    end


    if @query_strand == "+"
      @query_strand = :forward
    elsif @query_strand == "-"
      @query_strand = :reverse
    else
      raise  ExonerateException.new(), "Ivalid query orientation #{@query_strand} for line:\n#{sugar_str}"
    end

    raise  ExonerateException.new(), "Inconsistent orientation (forward, query)" if @query_strand == :forward and @query_start > @query_end
    raise  ExonerateException.new(), "Inconsistent orientation (reverse, query)" if @query_strand == :reverse and @query_start < @query_end
    raise  ExonerateException.new(), "Inconsistent orientation (forward, target)" if @target_strand == :forward and @target_start > @target_end
    raise  ExonerateException.new(), "Inconsistent orientation (reverse, target)" if @target_strand == :reverse and @target_start < @target_end


    self
  end


    #The vulgar has to be parsed AFTER the sugar, otherwise it is impossible to determine the orientations
    def parse_vulgar(vulgar_str)

      tarcurrent = @target_start
      query_current = @query_start
      target_multiply = 1
      query_multiply = 1

      if @target_strand == :reverse
        target_multiply = -1
      end

      if @query_strand == :reverse
        query_multiply = -1
      end

      @vulgar_block = Array.new
      #p "VULGAR #{vulgar_str}"
      vulgar_str.split(/\s/).each_slice(3) do | block |
        #p block
        vulgar = Vulgar.new(block[0].to_sym, block[1].to_i, block[2].to_i, tarcurrent, target_multiply, query_current, query_multiply, self)
        query_current = vulgar.query_end
        tarcurrent = vulgar.target_end
        vulgar_block << vulgar
      end
      self
    end

    #This assumes that the gene is the query and the chromosome is the target
    def exon_on_gene_position(position)
      @vulgar_block.each do |vulgar|
        if position.between?(vulgar.query_start, vulgar.query_end)
          return vulgar
        end
      end
      nil
    end

    def query_position_on_target(position, base:0)
      vulgar = exon_on_gene_position(position)
      qr = vulgar.query_region
      tr = vulgar.target_region
      
      offset = qr.orientation == :forward ? position - qr.start + 1 : qr.end - position

      #puts vulgar.to_s
      #puts "SNP position: #{position}"
      #puts vulgar.query_region
      #puts vulgar.query_region.orientation
      #puts "Offset query: #{offset}"
      #puts vulgar.target_region
      #puts vulgar.target_region.orientation

      new_pos = tr.orientation == :forward ? offset + tr.start - 1 :  tr.end - offset + 1

      return new_pos
    end

    def tarpostion_from_query_position(position)
      ret = nil
      vulgar_block = exon_on_gene_position(position)
      ret
    end

    def print_features
      out = String.new

      @vulgar_block.each do |  vulgar |
        out << vulgar.to_s << "\n"
      end
      out
    end
  end

  class Vulgar
    attr_reader :label, :query_length, :target_length, :query_start, :query_end, :target_start, :target_end, :record, :snp_in_gap
    def initialize(label, ql, tl, target_start, target_multiply, query_start, query_multiply, record)
      @label = label
      @query_length = ql
      @target_length = tl
      @query_start = query_start
      @query_end = query_start + (query_multiply   * query_length)
      @target_start = target_start
      @target_end = target_start + (target_multiply * target_length)
      @record = record
      @snp_in_gap = false
    end    

    def to_s
      out = String.new
      out << @label.to_s << "\t" << @query_length.to_s << "\t" << @target_length.to_s << "\t" << @query_start.to_s << "\t" << @query_end.to_s << "\t" << @target_start.to_s << "\t" << @target_end.to_s 
      out
    end

    def query_id
      record.query_id
    end

    def target_id
      record.target_id
    end

    def target_flanking_region_from_position(position, flanking_size)
      reg = reg = Bio::DB::Fasta::Region.new()
      reg.entry = target_id
      target_snp_pos = target_position_from_query(position)
      return nil if snp_in_gap 
      reg.orientation = record.target_strand
      reg.start = target_snp_pos - flanking_size
      reg.end = target_snp_pos + flanking_size
      raise  ExonerateException.new "Target Query out of bounds!" unless position.between?(query_start, query_end)

      reg
    end

    def target_position_from_query(position)
      raise ExonerateException.new(), "Position: #{position} not in range (#{query_start}-#{query_end}) #{self.to_s} " unless position.between?(query_start, query_end) or position.between?(query_end, query_start) 
      offset = 0
      ret = 0
      if record.query_strand == :forward
        offset = position - query_start
      elsif record.query_strand == :reverse
        offset = query_start - position
      else
        raise ExonerateException.new(), "The strand is not forward or reverse (#{record.query_strand}) ! #{self.inspect}"
      end

      if record.target_strand == :forward
        ret = target_start + offset
      elsif record.target_strand == :reverse
        ret = target_start - offset + 1
      else
        raise ExonerateException.new(), "The strand is not forward or reverse! #{self.inspect}"
      end
      #THis is in case the position is on a gap. 
      if @target_length == 0 and label == :G
        @snp_in_gap = true
        ret = target_start
      end
      raise ExonerateException.new(), "Return position #{ret} outside block (#{target_start}-#{target_end}, #{self.inspect})" unless ret.between?(target_start, target_end) or ret.between?(target_end, target_start)
      ret
    end

    def query_region
      reg = Bio::DB::Fasta::Region.new()
      reg.entry = query_id
      reg.orientation = record.query_strand 
      if record.query_strand == :forward
        reg.start = @query_start + 1
        reg.end =  @query_end 
      elsif record.query_strand == :reverse
        reg.start = @query_end + 1
        reg.end =  @query_start 
      else
        raise  ExonerateException.new(), "Ivalid query orientation #{@query_strand}"
      end
      reg
    end

    def target_region
      reg = Bio::DB::Fasta::Region.new()

      reg.entry = target_id
      reg.orientation = record.target_strand 
      if record.target_strand == :forward
        reg.start = @target_start + 1
        reg.end =  @target_end 
      elsif record.target_strand == :reverse
        reg.start = @target_end + 1
        reg.end =  @target_start 
      else
        raise  ExonerateException.new(), "Ivalid target orientation #{@target_strand}"
      end
      reg
    end

  end

end
