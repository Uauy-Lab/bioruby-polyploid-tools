

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

class Bio::Blat::Report::Hit
  
  #Function to parse stuff like: IWGSC_CSS_1AL_scaff_110
  def wheat_chr_arm
    @wheat_chr_arm if @wheat_chr_arm
    @wheat_chr_arm = target_id.split('_')[2]
  end
  
  def wheat_chr
    wheat_chr_arm[0,2]
  end
  
  def wheat_chr_group 
    raise Exception.new(), "No wheat group for #{target_id} #{self.inspect}"  unless wheat_chr
    wheat_chr_arm[0]
  end
  
  def wheat_genome
    wheat_chr_arm[1]
  end
  
  def wheat_arm
    wheat_chr_arm[2]
  end
  
  def percentage_covered
    ( match + mismatch ) * 100.0 / query_len.to_f
  end
  
end


class Hash
  def join(keyvaldelim=$,, entrydelim=$,)
    map {|e| e.join(keyvaldelim) }.join(entrydelim)
  end
end


class Bio::NucleicAcid

  IUPAC_CODES ||= {

    'y' => 'ct',
    'r' => 'ag',
    'w' => 'at',
    's' => 'cg',
    'k' => 'gt',
    'm' => 'ac',

    'b' => 'cgt',
    'd' => 'agt',
    'h' => 'act',
    'v' => 'acg',

    'n' => 'acgt',

    'a' => 'a',
    't' => 't',
    'g' => 'g',
    'c' => 'c',
    'u' => 'u',

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

#Monkey patching to Bio::Sequence to find snps between sequences. It assumes the
#sequences are already aligned and doesn't check if a base on the first sequence is
#valid on the second. 
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
  #Monkey patching to count how many ambiguity codes are present in the string, for Nucleic Acids
  def count_ambiguities
    snps=0

    for i in (0..self.size-1)

      snps += 1 if !Bio::NucleicAcid.is_unambiguous(self[i])
    end
    snps
  end
  
  #Counts how many bases are uppercase
  def upper_case_count
    match(/[^A-Z]*/).to_s.size
  end
end
