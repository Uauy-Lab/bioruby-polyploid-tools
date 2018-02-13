class Array
  def sum
    inject(0.0) { |result, el| result + el }
  end

  def mean 
    sum / size
  end
end

module Bio::PolyploidTools::Mask

  def self.get(seqs, target: nil)
    names = seqs.keys
    target = names[0] if target.nil? 
    masked_snps = seqs[target].downcase
    i = 0
    while i < masked_snps.size
      different = 0
      cov = 0
      gap = false
      names.each do | chr |
        if seqs[chr] and seqs[chr][i]  != "-"
          cov += 1 
          if chr != target 
            different += 1  if masked_snps[i].upcase != seqs[chr][i].upcase 
          end
        elsif seqs[chr] and seqs[chr][i]  == "-" and chr == target
            gap = true
        end
      end
      masked_snps[i] = "-" if gap 
      masked_snps[i] = "." if different == 0
      masked_snps[i] = "." if cov == 1
      masked_snps[i] = "*" if cov == 0
      expected_snps = names.size - 1
      masked_snps[i] = masked_snps[i].upcase if different == expected_snps
      i += 1
    end
    masked_snps
  end

  def self.stats(mask)
    specific = []
    semispecific = []
    sp_i = 0
    semi = 0
    i = 0
    mask.to_s.each_char do |e|  
      case e
      when /[[:lower:]]/ then 
        semispecific << semi
        semi = 0
      when /[[:upper:]]/ then 
        specific     << sp_i
        semispecific << semi
        sp_i = 0
        semi = 0
      else
        semi += 1
        sp_i += 1
      end 
    end

     
    {
      semispecific: semispecific.mean, 
      specific: specific.mean
    }
  end
end