
$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
$: << File.expand_path('.')
path= File.expand_path(File.dirname(__FILE__) + '/../lib/bioruby-polyploid-tools.rb')

tmp_verb = $VERBOSE
$VERBOSE=nil
#puts path
require path
require 'bio-samtools'
require "test/unit"
$VERBOSE=tmp_verb


class TestPolyploidTools < Test::Unit::TestCase

  
  #Set up the paths
  def setupre
    data_path= File.expand_path(File.dirname(__FILE__)   + '/data/' )
    
    @ref=data_path + '/S22380157.fa'
    @a=data_path + "/LIB1721.bam"
    @b=data_path + "/LIB1722.bam"
    @f2_a=data_path + "/LIB1716.bam"
    @f2_b=data_path + "/LIB1719.bam"
    
    @bfr_path=data_path + "/bfr_out_test.csv"
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
    #puts region.to_s
    
    #puts @bam_a.methods
    ref_seq=@fasta_db.fetch_sequence(region)
    reg_a = @bam_a.fetch_region({:region=>region,  :min_cov=>min_cov, :A=>1})
    reg_b = @bam_b.fetch_region({:region=>region,  :min_cov=>min_cov, :A=>1})
    cons_1 = reg_a.consensus
    cons_2 = reg_b.consensus
    
    snps_1 = cons_1.count_ambiguities
    snps_2 = cons_2.count_ambiguities
    
    called_1 = reg_a.called
    called_2 = reg_b.called
    
    snps_tot = Bio::Sequence.snps_between(cons_1, cons_2)
    #block_size = 1000
    #snps_per_1k_1   = (block_size * snps_1.to_f   ) / region.size
    #snps_per_1k_2   = (block_size * snps_2.to_f   ) / region.size
    #snps_per_1k_tot = (block_size * snps_tot.to_f ) / region.size
        
    snps_tot = Bio::Sequence.snps_between(cons_1, cons_2)
    snps_to_ref = Bio::Sequence.snps_between(cons_1, ref_seq)

    assert_equal(ref_seq.to_s, "acgcttgaccttaggcctatttaggtgacactatagaacaagtttgtacaaaaaagcaggctggtaccggtccggaattcccgggatatcgtcgacccacgcgtccgcgtccgaccagcacaaacaagactgtactctgggctcctctgactccgtgtcttgctaaaatatctttggtcgactcgttgcgaggttgatcagatggcggaggaagcgaagcaggatgtggcgccacccgcgccggagccgaccgaggacgtcgcggacgagaaggtggcggttccgtcgccggaggagtctaaggccctcgttgtcgccgagaatgacgctgagaagcctgcagctacagggggctcacacgaacgagatgctctgctcacgagggtcgcgaccgagaagaggatttcgctgatcaaggcatgggaggagaacgagaaggccaaagccgagaacaaggccgtgaagttgctggcggacatcacctcgtgggagaactccaaggccgcggaactggaagccgagctcaagaagatgcaagagcagctggagaagaagaaggcgcgctgcgtggagaagctcaagaacagcgccgcgacggtgcacaaagaggcggaangagaagcgtgccgcggcggaagcgcggcacggcgaggagatcgtcgcggcggaggagaccgccgccaagtaccgcgccaagggtgaagcgccgaagaagctgctcttcggcagaagatagatatcgcttcatcttcagcttctctctgtttgaccgnttgcatgtctcctgcccatggcatcacttgtgtatttatctttgggggngatcttagtttgtatggtatcatcaaatgcgtcgtga")
    assert_equal(cons_1.to_s , "acgcttgaccttaggcctatttaggtgacactatagaacaagtttgtacaaaaaagcaggctggtaccggtccggaattcccgggatatcgtcgacccacgcgtccgcgtccgaccagcacaaacaagactgtactctgggctcctctgactccgtgtcttgctaaaatatytttggtcgactcgttgcgaggttgatcagatggcggaggaagcgaagcaggatgtggcgccacccgcgccggagccgaccgaggacgtcgcggacgagaaggcggcggttccgtcgccggaggagtctaaggccctsgttgtcgccgagaatgacgcygagaagcctgcagctacagggggctcacacgaacgagatgctctgctcacgagggtygcgaccgagaagaggatttcgctgatcaaggcatgggaggagaaygagaaggccaaagccgagaacaaggccgtgaagttgctggcggacatcacctcgtgggagaactccaaggccgcggaactggaagccgagctcaagaagatgcaagagcagctggagaagaagaaggcgcgctgcgtggagaagctcaagaacagcgccgcgacggtgcacaaagaggcgraaggagaagcgtgccgcggcggaagygcggcrcggcgaggagatcgtcgcggcggaggagaccgccgccaagtaccgcgccaagggtgaggcgccgaagaagctgctcttcggcagaggatagatatcgcttcatcttcagcttctctctgtttgaccgnttgcatgtctcctgcccatggcatcacttgtgtatttatctttgggggngatcttagtttgtatggtatcatcaaatgcgtcgtga")  
    assert_equal(cons_2.to_s , "acgcttgaccttaggcctatttaggtgacactatagaacaagtttgtacaaaaaagcaggctggtaccggtccggaattcccgggatatcgtcgacccacgcgtccgcgtccgaccagcacaaacaagactgtactctgggctcctctgactccgtgtcttgctaaaatatytttggtcgactcgttgcgaggttgatcagatggcggasgaagcgaagcaggatgtggcgccacccgcgccggagccgaccgaggacgtcgcggacgagaaggcggcggttccgtcgccggaggartcyaaggccctsgttgtcgccgagaatgacgcygagaagcctgcagctacagggggctcacacgaacgagatgctctgctcacgagggtygcgaccgagaagaggatttcgctgatcaaggcatgggaggagaaygagaaggccaaagccgagaacaaggccgtgaagttgctggcggacatcacctcgtgggagaactccaaggccgcggaactggaagccgagctcaagaagatgcaagagcagctggagaagaagaaggcgcgctgcgtggagaagctcaagaacagcgccgcgacggtgcacaaagaggcgraaggagaagcgtgccgcggcggaagygcggcgcggcgaggagatcgtcgcggcggaggagrccgccgccaagtaccgcgccaagggtgaggcgccgaagaagctgctcttcggcagaagatagatatcgcttcatcttcagcttctctctgtttgaccgnttgcatgtctcctgcccatggcatcacttgtgtatttatctttgggggngatcttagtttgtatggtatcatcaaatgcgtcgtga")
    assert_equal(snps_tot , 6)
    assert_equal(snps_to_ref , 12)
    assert_equal(snps_1,10)
    assert_equal(snps_2,13)
    assert_equal(called_1,617)
    assert_equal(called_2,612)
  end
  
  def test_bfr
    setupre
    container = Bio::BFRTools::BFRContainer.new
    
    container.reference @ref
    container.parental_1  ( {:path => @a } )
    container.parental_2  ( {:path => @b } )
    container.bulk_1 ( {:path => @f2_a  })
    container.bulk_2 ( {:path => @f2_b  })

    i = -1

    container.init_counters
    output_file =  File.open(@bfr_path, "w")
  #  puts "Range: #{min}:#{max}"
    assert_equal(@fasta_db.index.entries.size,1)
    reg = nil
    @fasta_db.index.entries.each do | r |
      i = i  + 1
       
      reg = container.process_region({:A => true, :q => 37, :region => r.get_full_region.to_s,:output_file => output_file , :min_cov => 5} )
      #puts reg.inspect
    end
    
    with_bfr = [210, 297, 300, 645, 674]
    
    bases_1 = Array.new
    bases_2 = Array.new
    bases_1  << {:A=>0, :C=>24, :G=>147, :T=>0}
    bases_2  << {:A=>0, :C=>33, :G=>142, :T=>0}
    bases_1  << {:A=>35, :C=>0, :G=>142, :T=>0}
    bases_2  << {:A=>25, :C=>0, :G=>149, :T=>0}
    bases_1  << {:A=>0, :C=>32, :G=>0, :T=>146}
    bases_2  << {:A=>0, :C=>24, :G=>0, :T=>152}
    bases_1  << {:A=>22, :C=>0, :G=>38, :T=>0}
    bases_2  << {:A=>49, :C=>0, :G=>20, :T=>0}
    bases_1  << {:A=>26, :C=>0, :G=>16, :T=>0}
    bases_2  << {:A=>36, :C=>0, :G=>10, :T=>0}
    i = 0
    with_bfr.each do | pos |
      #puts pos
      assert_equal(bases_1[i], reg.bases_bulk_1[pos - 1 ] )
      assert_equal(bases_2[i], reg.bases_bulk_2[pos - 1 ] )
      i += 1
    end
   
   
    
    output_file.close
    
  end
  
end