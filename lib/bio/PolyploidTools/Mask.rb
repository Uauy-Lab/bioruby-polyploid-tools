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
        if seqs[chr][i]  != "-" and seqs[chr][i]  != "n" and seqs[chr][i]  != "N"
          cov += 1
          if chr != target
            different += 1  if masked_snps[i].upcase != seqs[chr][i].upcase
          end
        end
        if seqs[chr][i]  == "-" and chr == target
            gap = true
        end
      end
      masked_snps[i] = "." if different == 0
      masked_snps[i] = "." if cov == 1
      masked_snps[i] = "*" if cov == 0
      expected_snps  = names.size - 1
      masked_snps[i] = masked_snps[i].upcase if different == expected_snps
      masked_snps[i] = "-" if gap
      i += 1
    end
    masked_snps
  end

  def self.stats(mask, triad, gene, genome, reference)
    specific = []
    semispecific = []
    sp_i = 0
    semi = 0
    i = 0
    mask.to_s.each_char do |e|
      case e
      when "n","N"
        i += 1
      when /[[:lower:]]/ then
        semispecific << semi
        semi = 0
        i += 1
      when /[[:upper:]]/ then
        specific     << sp_i
        semispecific << semi
        sp_i = 0
        semi = 0
        i += 1
      when "." then
        semi += 1
        sp_i += 1
        i += 1
      end
    end
    {
      reference: reference,
      triad: triad,
      genome: genome,
      gene: gene,
      semispecific: semispecific.mean,
      semispecific_bases: semispecific.size,
      semispecific_identity: (1 - (semispecific.size.to_f / i)) * 100 ,
      specific: specific.mean,
      specific_bases: specific.size,
      specific_identity: (1 - (specific.size.to_f / i )) * 100,
      length: i
    }
  end
end
