require 'rubygems'
require 'pathname'
#require_relative 'db/fasta.rb'
require 'bio'

require_relative 'db/fastadb.rb'

#require "set"
#require 'systemu'
#require 'json'

=begin

Extends the methods to be able to calculate the BFR and a consensus from the pileup

=end

class Bio::DB::Pileup

  #attr_accessor :minumum_ratio_for_iup_consensus
  #@minumum_ratio_for_iup_consensus = 0.20

  #Returns a hash with the count of bases

  def bases
    return @bases if @bases
    @bases = self.non_refs
    #puts self.ref_count
    @bases[self.ref_base.upcase.to_sym] = self.ref_count 
    @bases
  end

  def base_coverage
    total = 0
    @bases.each do |k,v|
      total += v  
    end
    total
  end

  def base_ratios
    return @base_ratios if @base_ratios
    bases = self.bases
    @base_ratios = Hash.new
    bases.each do |k,v| 
      @base_ratios[k] = v.to_f/self.base_coverage.to_f 
    end
    @base_ratios
  end

  # returns the consensus (most frequent) base from the pileup, if there are equally represented bases returns a string of all equally represented bases in alphabetical order   
  def consensus_iuap(minumum_ratio_for_iup_consensus)
    minumum_ratio_for_iup_consensus
    if @consensus_iuap.nil?
      @consensus_iuap = self.ref_base.downcase
      bases = self.bases
      tmp = String.new
      bases.each do |k,v|
        tmp << k[0].to_s if v/self.coverage > minumum_ratio_for_iup_consensus
      end
      if tmp.length > 0
        @consensus_iuap = Bio::NucleicAcid.to_IUAPC(tmp)
      end
    end 
    @consensus_iuap
  end
end

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


#class Bio::DB::Sam::SAMException < RuntimeError

#end

class Bio::DB::Sam


  attr_accessor :minumum_ratio_for_iup_consensus
  attr_reader :cached_regions
  #attr_accessor :pileup_cache
  @minumum_ratio_for_iup_consensus = 0.20


  #Same as mpilup, but it caches the pileup, so if you want several operations on the same set of regions
  #the pile for different operations, it won't execute the mpilup command several times
  #Whenever you finish using a region, call mpileup_clear_cache to free the cache
  #The argument Region is required, as it will be the key for the underlying hash. 
  #We asume that the options are constant. If they are not, the cache mechanism may not be consistent. 
  #
  #TODO: It may be good to load partially the pileup
  def mpileup_cached (opts={})      
    raise SAMException.new(), "A region must be provided" unless opts[:r] or opts[:region]
    @pileup_cache = Hash.new unless @pileup_cache
    @cached_regions = Hash.new unless @cached_regions

    region = opts[:r] ? opts[:r] : opts[:region]
    opts[:r] = "'#{region.to_s}'"
    opts[:region] = "'#{region.to_s}'"
    opts[:A] = true
    #reg = region.class == Bio::DB::Fasta::Region ? region : Bio::DB::Fasta::Region.parse_region(region.to_s)

    unless @cached_regions[region.to_s]
      @cached_regions[region.to_s] =  Bio::DB::Fasta::Region.parse_region(region.to_s)
      tmp = Array.new
      @cached_regions[region.to_s].pileup =  tmp
      #puts "Loading #{region.to_s}"
      mpileup(opts) do | pile | 
        #   puts pile
        tmp << pile 
        yield pile
      end
    else   
      #   puts "Loaded, reruning #{region.to_s}"
      @cached_regions.pileup[region.to_s] .each do | pile |
        yield pile
      end
    end
  end

  #Clears the pileup cache. If a region is passed as argument, just the specified region is removed
  #If no region is passed, the hash is emptied
  def mpileup_clear_cache (region)
    return unless @cached_regions
    if region
      @cached_regions[region.to_s] = nil
    else
      @cached_regions.clear
    end
  end

  #Gets the coverage of a region from a pileup. 
  def average_coverage_from_pileup(opts={})
    opts[:region] =   opts[:region].to_s if opts[:region] .class == Bio::DB::Fasta::Region 
    region = opts[:region]
    calculate_stats_from_pile(opts) if @cached_regions == nil or @cached_regions[region] == nil
    @cached_regions[region].average_coverage
  end

  #
  def coverages_from_pileup(opts={})
    opts[:region] =   opts[:region].to_s if opts[:region] .class == Bio::DB::Fasta::Region 
    region = opts[:region]
    calculate_stats_from_pile(opts) if @cached_regions == nil or @cached_regions[region] == nil
    @cached_regions[region].coverages
  end

  def consensus_with_ambiguities(opts={})
    opts[:region] =   opts[:region].to_s if opts[:region] .class == Bio::DB::Fasta::Region 
    region = opts[:region]
    #   p "consensus with ambiguities for: " << opts[:region] 
    calculate_stats_from_pile(opts) if @cached_regions == nil or @cached_regions[region] == nil
    @cached_regions[region].consensus
  end

  def calculate_stats_from_pile(opts={})
    min_cov = opts[:min_cov] ? opts[:min_cov] : 20  


    opts[:region] = Bio::DB::Fasta::Region.parse_region( opts[:region] .to_s)  unless opts[:region].class == Bio::DB::Fasta::Region
    region = opts[:region]
    
    mark_case = true if opts[:case]
   # puts "Marcase: #{mark_case}"
    reference = self.fetch_reference(region.entry, region.start, region.end).downcase
    #  p "calculationg from pile..." << region.to_s
    base_ratios = Array.new(region.size, BASE_COUNT_ZERO) 
    bases = Array.new(region.size, BASE_COUNT_ZERO) 
    coverages = Array.new(region.size, 0)
    total_cov = 0
    
    self.mpileup_cached(:region=>"#{region.to_s}") do | pile |
      #puts pile
      #puts pile.coverage
      if pile.coverage > min_cov
        base_ratios[pile.pos - region.start ] = pile.base_ratios
        reference[pile.pos - region.start  - 1 ] = pile.consensus_iuap(0.20)
        coverages[pile.pos - region.start   ]  = pile.coverage.to_i
        bases[pile.pos - region.start   ]  = pile.bases
      end
      total_cov += pile.coverage
    end

    #puts ">Ref\n#{reference}"
    region = @cached_regions[region.to_s]
    region.coverages = coverages
    region.base_ratios = base_ratios
    region.consensus = Bio::Sequence.new(reference)
    region.consensus.na
    if region.orientation == :reverse
      region.consensus.reverse_complement!()
    end
    region.average_coverage = total_cov.to_f/region.size.to_f
    region.bases = bases
    region
  end



  BASE_COUNT_ZERO =  {:A => 0, :C => 0, :G => 0,  :T => 0}

  #Gets an array with the proportions of the bases in the region. If there is no coverage, a
  def base_ratios_in_region(opts={})
    opts[:region] =   opts[:region].to_s if opts[:region] .class == Bio::DB::Fasta::Region 
    region = opts[:region]
    calculate_stats_from_pile(opts) if @cached_regions == nil or @cached_regions[region] == nil
    @cached_regions[region].base_ratios 
  end

  #Gets an array with the bsaes count in the region. If there is no coverage, a
  def bases_in_region(opts={})
    opts[:region] =   opts[:region].to_s if opts[:region] .class == Bio::DB::Fasta::Region 
    region = opts[:region]
    calculate_stats_from_pile(opts) if @cached_regions == nil or @cached_regions[region] == nil
    @cached_regions[region].bases 
  end



  def extract_reads(opts={})
    opts[:region] = Bio::DB::Fasta::Region.parse_region( opts[:region] .to_s)  unless opts[:region].class == Bio::DB::Fasta::Region
    fastq_filename = opts[:fastq]
    fastq_file = opts[:fastq_file]

    out = $stdout

    print_fastq = Proc.new do |alignment|
      out.puts "@#{alignment.qname}"
      out.puts "#{alignment.seq}"
      out.puts "+#{alignment.qname}"
      out.puts "#{alignment.qual}"
    end

    fetch_with_function(chromosome, qstart, qstart+len,  print_fastq)


  end

end

class Bio::DB::Fasta::Region
  attr_accessor :pileup, :average_coverage, :snps, :reference, :base_ratios, :consensus, :coverages, :bases

  #TODO: Debug, as it hasnt been tested in the actual code. 
  def base_ratios_for_base(base)
    @all_ratios = Hash.new unless @all_ratios
    unless @all_ratios[base]
      ratios = Array.new
      for i in (0..region.size-1)
        ratios << @base_ratios[i][base]
      end
      @all_ratios[base] = ratios
    end
    @all_ratios[base]
  end

end