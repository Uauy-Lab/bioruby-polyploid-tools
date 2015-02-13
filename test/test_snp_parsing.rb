$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
$: << File.expand_path('.')
path= File.expand_path(File.dirname(__FILE__) + '/../lib/bioruby-polyploid-tools.rb')

#puts path
require path
require "test/unit"

class TestPolyploidTools < Test::Unit::TestCase

  #Set up the paths
  def setup
    @data = File.expand_path(File.dirname(__FILE__) + "/data")
    
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

  def test_mutant_snp

    ref=@data + "/IWGSC_CSS_1AL_scaff_1455974_aln_contigs.fa"
    
    fasta_reference_db = Bio::DB::Fasta::FastaFile.new({:fasta=>ref})
    fasta_reference_db.load_fai_entries 
    
    snp = Bio::PolyploidTools::SNPMutant.parse("IWGSC_CSS_1AL_scaff_1455974,Kronos2281,127,C,T")
    assert_equal(snp.gene , "IWGSC_CSS_1AL_scaff_1455974", "The original name was not parsed: #{snp.gene}")
    assert_equal(snp.chromosome, "IWGSC_CSS_1AL_scaff_1455974", "The chromosome wasnt parsed: #{snp.chromosome}")
    assert_equal(snp.position, 127, "The position is not parsed: #{snp.position}")
    assert_equal(snp.chr, "1A")
    
    region = fasta_reference_db.index.region_for_entry(snp.gene).get_full_region
    snp.full_sequence = fasta_reference_db.fetch_sequence(region)

    assert_equal(snp.template_sequence, "actcgatcgtcagcacccgctggaacttggggaacgtcttgaacgccgcaagcaccggggcgtcctctgactgtatgagcacgcgctgcttacaggtctcYttgtcgtacccggacttgacaagcgctttggagaccgcatccaccacgtcaaggcttctggctataaggtacgtagcatgctgcactcggtaggtacaaga")
    assert_equal(snp.sequence_original, "actcgatcgtcagcacccgctggaacttggggaacgtcttgaacgccgcaagcaccggggcgtcctctgactgtatgagcacgcgctgcttacaggtctc[C/T]ttgtcgtacccggacttgacaagcgctttggagaccgcatccaccacgtcaaggcttctggctataaggtacgtagcatgctgcactcggtaggtacaaga")
    assert_equal(snp.position, 101)
    assert_equal(snp.original, "C")
    assert_equal(snp.snp, "T")
    
    
  end

  
end