
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
    #puts "TADA"    
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
