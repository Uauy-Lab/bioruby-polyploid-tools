$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
$: << File.expand_path('.')
path= File.expand_path(File.dirname(__FILE__) + '/../lib/bioruby-polyploid-tools.rb')

#puts path
require path
require "test/unit"

class TestPolyploidTools < Test::Unit::TestCase

  
  #Set up the paths
  def setupre
    data_path= File.expand_path(File.dirname(__FILE__)   + '/data/' )
    
    @ref=data_path + '/S22380157.fa'
    @a=data_path + "/LIB1721.bam"
    @b=data_path + "/LIB1722.bam"
    @f2_a=data_path + "/LIB1716.bam"
    @f2_b=data_path + "/LIB1719.bam"
    @fasta_db = Bio::DB::Fasta::FastaFile.new({:fasta=>@ref})
    @fasta_db.load_fai_entries
    @bam_a =  Bio::DB::Sam.new({:fasta=>@ref, :bam=>@a})
    @bam_b =  Bio::DB::Sam.new({:fasta=>@ref, :bam=>@b})
    @bam_f2_a =  Bio::DB::Sam.new({:fasta=>@ref, :bam=>@f2_a})
    @bam_f2_b =  Bio::DB::Sam.new({:fasta=>@ref, :bam=>@f2_b})
   # puts "SETUP"
  end
  
  def teardown
    
  end
  
  def test_snp_between_consensus
    setupre
    
    reg="gnl|UG|Ta#S22380157"
    region = @fasta_db.index.region_for_entry(reg).to_region
    min_cov=20
    
    ref_seq=@fasta_db.fetch_sequence(region)
    puts region.to_s
    
    #puts @bam_a.methods
    cons_1 = @bam_a.consensus_with_ambiguities({:region=>region, :case=>true, :min_cov=>min_cov})
    cons_2 = @bam_b.consensus_with_ambiguities({:region=>region, :case=>true, :min_cov=>min_cov})
    
    snps_tot = Bio::Sequence.snps_between(cons_1, cons_2)
    snps_to_ref = Bio::Sequence.snps_between(cons_1, ref_seq)
    puts ">ref\n#{ref_seq}"
    puts ">a\n#{cons_1}"
    puts ">b\n#{cons_2}"
    puts "SNPS between: #{snps_tot}"
    puts "SNPS ref: #{snps_to_ref}"
    assert(ref_seq.to_s=="acgcttgaccttaggcctatttaggtgacactatagaacaagtttgtacaaaaaagcaggctggtaccggtccggaattcccgggatatcgtcgacccacgcgtccgcgtccgaccagcacaaacaagactgtactctgggctcctctgactccgtgtcttgctaaaatatctttggtcgactcgttgcgaggttgatcagatggcggaggaagcgaagcaggatgtggcgccacccgcgccggagccgaccgaggacgtcgcggacgagaaggtggcggttccgtcgccggaggagtctaaggccctcgttgtcgccgagaatgacgctgagaagcctgcagctacagggggctcacacgaacgagatgctctgctcacgagggtcgcgaccgagaagaggatttcgctgatcaaggcatgggaggagaacgagaaggccaaagccgagaacaaggccgtgaagttgctggcggacatcacctcgtgggagaactccaaggccgcggaactggaagccgagctcaagaagatgcaagagcagctggagaagaagaaggcgcgctgcgtggagaagctcaagaacagcgccgcgacggtgcacaaagaggcggaangagaagcgtgccgcggcggaagcgcggcacggcgaggagatcgtcgcggcggaggagaccgccgccaagtaccgcgccaagggtgaagcgccgaagaagctgctcttcggcagaagatagatatcgcttcatcttcagcttctctctgtttgaccgnttgcatgtctcctgcccatggcatcacttgtgtatttatctttgggggngatcttagtttgtatggtatcatcaaatgcgtcgtga")
    assert(cons_1.to_s=="acgcttgaccttaggcctatttaggtgacactatagaacaagtttgtacaaaaaagcaggctggtaccggtccggaattcccgggatatcgtcgacccacgcgtccgcgtccgaccagcacaaacaagactgtactctgggctcctctgactccgtgtcttgctaaaatatytttggtcgactcgttgcgaggttgatcagatggcggaggaagcgaagcaggatgtggcgccacccgcgccggagccgaccgaggacgtcgcggacgagaaggcggcggttccgtcgccggaggagtctaaggccctsgttgtcgccgagaatgacgcygagaagcctgcagctacagggggctcacacgaacgagatgctctgctcacgagggtygcgaccgagaagaggatttcgctgatcaaggcatgggaggagaaygagaaggccaaagccgagaacaaggccgtgaagttgctggcggacatcacctcgtgggagaactccaaggccgcggaactggaagccgagctcaagaagatgcaagagcagctggagaagaagaaggcgcgctgcgtggagaagctcaagaacagcgccgcgacggtgcacaaagaggcgraaggagaagcgtgccgcggcggaagygcggcrcggcgaggagatcgtcgcggcggaggagaccgccgccaagtaccgcgccaagggtgaggcgccgaagaagctgctcttcggcagaggatagatatcgcttcatcttcagcttctctctgtttgaccgnttgcatgtctcctgcccatggcatcacttgtgtatttatctttgggggngatcttagtttgtatggtatcatcaaatgcgtcgtga")  
    assert(cons_2.to_s=="acgcttgaccttaggcctatttaggtgacactatagaacaagtttgtacaaaaaagcaggctggtaccggtccggaattcccgggatatcgtcgacccacgcgtccgcgtccgaccagcacaaacaagactgtactctgggctcctctgactccgtgtcttgctaaaatatytttggtcgactcgttgcgaggttgatcagatggcggasgaagcgaagcaggatgtggcgccacccgcgccggagccgaccgaggacgtcgcggacgagaaggcggcggttccgtcgccggaggartcyaaggccctsgttgtcgccgagaatgacgcygagaagcctgcagctacagggggctcacacgaacgagatgctctgctcacgagggtygcgaccgagaagaggatttcgctgatcaaggcatgggaggagaaygagaaggccaaagccgagaacaaggccgtgaagttgctggcggacatcacctcgtgggagaactccaaggccgcggaactggaagccgagctcaagaagatgcaagagcagctggagaagaagaaggcgcgctgcgtggagaagctcaagaacagcgccgcgacggtgcacaaagaggcgraaggagaagcgtgccgcggcggaagygcggcgcggcgaggagatcgtcgcggcggaggagrccgccgccaagtaccgcgccaagggtgaggcgccgaagaagctgctcttcggcagaagatagatatcgcttcatcttcagcttctctctgtttgaccgnttgcatgtctcctgcccatggcatcacttgtgtatttatctttgggggngatcttagtttgtatggtatcatcaaatgcgtcgtga")
    assert(snps_tot == 6)
    assert(snps_to_ref == 12)
     

  end
  
end