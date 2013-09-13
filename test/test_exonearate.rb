$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
$: << File.expand_path('.')
path= File.expand_path(File.dirname(__FILE__) + '/../lib/bioruby-polyploid-tools.rb')

#puts path
require path
require "test/unit"

class TestPolyploidTools < Test::Unit::TestCase
  Query=File.dirname(__FILE__) + '/data/'+"BS00068396_51.fa"
  Target=File.dirname(__FILE__) + '/data/'+"BS00068396_51_contigs.fa"
  #Set up the paths
  def setup
    File.expand_path(File.dirname(__FILE__) + '/data/')
    
    puts "SEETING UP *******************"
  end
  
  def teardown
    
  end
  
  def test_simple_align_array
 #   puts $LOAD_PATH
    alignments = Bio::DB::Exonerate.align({:query=>Query, :target=>Target})
    assert(alignments.size == 4, "The count of alignments should be 4, it was #{alignments.size}")
  end

  def test_parse_alingn_line
    line="RESULT:\tBS00068396_51 0 101 + 2BS_5163353 7425 7323 - 462\t96.04\t101\t11974\t.\tM 69 69 G 0 1 M 32 32"
    aln =  Bio::DB::Exonerate::Alignment.parse_custom(line)
    assert(aln.query_id == "BS00068396_51")
    assert(aln.query_start==0)
    assert(aln.query_end==101)
    assert(aln.query_strand==:forward)
    assert(aln.target_id=="2BS_5163353")
    assert(aln.target_start==7425) 
    assert(aln.target_end==7323) 
    assert(aln.target_strand==:reverse)
    assert(aln.score==462.0) 
    assert(aln.pi==96.04)
    assert(aln.ql==101)
    assert(aln.tl==11974)
    assert(aln.g==".")
    assert(aln.vulgar_block.to_s=="[M	69	69	0	69	7425	7356, G	0	1	69	69	7356	7355, M	32	32	69	101	7355	7323]" )
    assert(aln.line==line)
    
    puts aln.vulgar_block.inspect
    #puts aln.inspect
  end

  
end