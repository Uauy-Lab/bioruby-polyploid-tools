module Bio::PolyploidTools

  def get_mask(seqs, target: nil)
    names = seqs.keys
    target = names[0] if target.nil? 
    masked_snps = seqs[target].downcase

    i = 0
    while i < masked_snps.size
      different = 0
      cov = 0
      names.each do | chr |
        if aligned_sequences[chr] and aligned_sequences[chr][i]  != "-"
          cov += 1 
          puts "Comparing #{chromosome_group} and #{chr[0]} as chromosomes"
          if chr != chromosome 
            different += 1  if masked_snps[i].upcase != aligned_sequences[chr][i].upcase 
          end
        end
      end
      masked_snps[i] = "-" if different == 0
      masked_snps[i] = "-" if cov == 1
      masked_snps[i] = "*" if cov == 0
      expected_snps = names.size - 1
      masked_snps[i] = masked_snps[i].upcase if different == expected_snps
      i += 1
    end
    masked_snps
  end
end