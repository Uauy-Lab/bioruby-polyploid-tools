

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


