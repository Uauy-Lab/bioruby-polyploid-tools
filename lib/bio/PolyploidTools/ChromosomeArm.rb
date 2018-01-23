class Bio::PolyploidTools::ChromosomeArm



  @@arm_selection_functions = Hash.new;

  #example format: chr2A
  @@arm_selection_functions[:nrgene] = lambda do | contig_name |
    ret = contig_name[3,2]
    return ret
  end

  @@arm_selection_functions[:first_two] = lambda do | contig_name |
    contig_name.gsub!(/chr/,"")
    ret = contig_name[0,2]       
    return ret
  end

  #Function to parse stuff like: "IWGSC_CSS_1AL_scaff_110"
  #Or the first two characters in the contig name, to deal with 
  #pseudomolecules that start with headers like: "1A"
  #And with the cases when 3B is named with the prefix: v443
  @@arm_selection_functions[:embl] = lambda do | contig_name|

    arr = contig_name.split('_')
    ret = "U"
    ret = arr[2][0,2] if arr.size >= 3
    ret = "3B" if arr.size == 2 and arr[0] == "v443"
    ret = arr[0][0,2] if arr.size == 1   
    return ret
  end

  @@arm_selection_functions[:morex] = lambda do | contig_name |
    ret = contig_name.split(':')[0].split("_")[1];       
    return ret
  end

  @@arm_selection_functions[:scaffold] = lambda do | contig_name |
    ret = contig_name;       
    return ret
  end

  def self.getArmSelection(name)
    @@arm_selection_functions[name.to_sym]
  end

  def self.getValidFunctions
    @@arm_selection_functions.keys.map { |e| e.to_s }
  end

end