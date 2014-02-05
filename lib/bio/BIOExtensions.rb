
module Bio::NucleicAcid::Data
  IUPAC_CODES = {

    'y'	=> 'ct',
    'r'	=> 'ag',
    'w'	=> 'at',
    's'	=> 'cg',
    'k'	=> 'gt',
    'm'	=> 'ac',

    'b'	=> 'cgt',
    'd'	=> 'agt',
    'h'	=> 'act',
    'v'	=> 'acg',

    'n'	=> 'acgt',

    'a'	=> 'a',
    't'	=> 't',
    'g'	=> 'g',
    'c'	=> 'c',
    'u'	=> 'u',

    'ct' => 'y',
    'ag' => 'r',
    'at' => 'w',
    'cg' => 's',
    'gt' => 'k',
    'ac' => 'm',

    'cgt' => 'b',
    'agt' => 'd',
    'act' => 'h',
    'acg' => 'v',

    'acgt' => 'n'
  }


end

class Bio::NucleicAcid

  IUPAC_CODES = {

    'y'	=> 'ct',
    'r'	=> 'ag',
    'w'	=> 'at',
    's'	=> 'cg',
    'k'	=> 'gt',
    'm'	=> 'ac',

    'b'	=> 'cgt',
    'd'	=> 'agt',
    'h'	=> 'act',
    'v'	=> 'acg',

    'n'	=> 'acgt',

    'a'	=> 'a',
    't'	=> 't',
    'g'	=> 'g',
    'c'	=> 'c',
    'u'	=> 'u',

    'ct' => 'y',
    'ag' => 'r',
    'at' => 'w',
    'cg' => 's',
    'gt' => 'k',
    'ac' => 'm',

    'cgt' => 'b',
    'agt' => 'd',
    'act' => 'h',
    'acg' => 'v',

    'acgt' => 'n'
  }

  def self.is_unambiguous(base)
    "acgtACGT".match(base)
  end

  def self.to_IUAPC(bases)    
    base = IUPAC_CODES[bases.to_s.downcase.chars.sort.uniq.join]
    if base == nil
      p "Invalid base! #{base}"
      base = 'n' #This is a patch... as one of the scripts failed here. 
    end
    base.upcase
  end

  def self.is_valid(code, base)
    IUPAC_CODES[code.downcase].chars.include? base.downcase
  end

end

class Bio::Sequence
  def self.snps_between(seq1, seq2)
    snps=0
    for i in (0..seq1.size-1)
      snps += 1 if seq1[i] != seq2[i] 
    end
    snps
  end
end

class  String
  def count_ambiguities
    snps=0

    for i in (0..self.size-1)

      snps += 1 if !Bio::NucleicAcid.is_unambiguous(self[i]) and /[[:upper:]]/.match(self[i])  and /[[:upper:]]/.match(self[i])
    end
    snps
  end
  
  def upper_case_count
    match(/[^A-Z]*/).to_s.size
  end
end

class Bio::Blat
  def self.align(database , query , output)
    cmdline = "blat #{database} #{query}  #{output}"
    puts $stderr.puts cmdline
    status, stdout, stderr = systemu cmdline
    if status.exitstatus == 0
      alns = Array.new unless block_given?
      blat_aln = Bio::Blat::Report.new(Bio::FlatFile.open(output).to_io)
      #p blat_aln
      blat_aln.each_hit() do |hit|
        if block_given?
          yield hit
        else
          alns << hit
        end
      end
      return alns unless block_given?
    else
      raise Exception.new(), "Error running exonerate. Command line was '#{cmdline}'\nBlat STDERR was:\n#{stderr}"
    end
  end
end