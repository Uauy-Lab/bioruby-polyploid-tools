$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
$: << File.expand_path('.')
path= File.expand_path(File.dirname(__FILE__) + '/../lib/bioruby-polyploid-tools.rb')

#puts path
require path
require "test/unit"

class TestPolyploidTools < Test::Unit::TestCase

  #Set up the paths
  def setup
    
  end
  
  def teardown
    
  end
  
  
  def test_default
 #   puts $LOAD_PATH
    assert(true, "Unit test test")
  end
  
  def test_snp_sequence
    snp = Bio::PolyploidTools::SNPSequence.parse("BS00068396_51,2A,CGAAGCGATCCTACTACATTGCGTTCCTTTCCCACTCCCAGGTCCCCCTA[T/C]ATGCAGGATCTTGATTAGTCGTGTGAACAACTGAAATTTGAGCGCCACAA")
    assert(snp.gene == "BS00068396_51" )
    assert(snp.chromosome == "2A")
   
    assert(snp.sequence_original == "CGAAGCGATCCTACTACATTGCGTTCCTTTCCCACTCCCAGGTCCCCCTA[T/C]ATGCAGGATCTTGATTAGTCGTGTGAACAACTGAAATTTGAGCGCCACAA")
    
    assert_equal(snp.position ,  51, "Position isnt parsed #{snp.position}")
    assert_equal(snp.original , "T", "ORiginal base not parsed, is #{snp.original}")
    assert_equal(snp.snp , "C")
    assert(snp.template_sequence == "CGAAGCGATCCTACTACATTGCGTTCCTTTCCCACTCCCAGGTCCCCCTAYATGCAGGATCTTGATTAGTCGTGTGAACAACTGAAATTTGAGCGCCACAA", "#{snp.template_sequence}!=CGAAGCGATCCTACTACATTGCGTTCCTTTCCCACTCCCAGGTCCCCCTAYATGCAGGATCTTGATTAGTCGTGTGAACAACTGAAATTTGAGCGCCACAA")
    #true
  end
  
end