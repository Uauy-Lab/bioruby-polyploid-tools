$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
$: << File.expand_path('.')
path= File.expand_path(File.dirname(__FILE__) + '/../lib/bioruby-polyploid-tools.rb')

#puts path
require path
require "test/unit"

class TestSNPparsing < Test::Unit::TestCase

  #Set up the paths
  def setup
    @data = File.expand_path(File.dirname(__FILE__) + "/data")    
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
  end

  def test_mutant_snp

    ref=@data + "/IWGSC_CSS_1AL_scaff_1455974_aln_contigs.fa"
    fasta_reference_db = Bio::DB::Fasta::FastaFile.new({:fasta=>ref})
    fasta_reference_db.index
    fasta_reference_db.load_fai_entries 
    
    snp = Bio::PolyploidTools::SNPMutant.parse("IWGSC_CSS_1AL_scaff_1455974,Kronos2281,127,C,T")
    assert_equal(snp.gene , "1AL_1455974_Kronos2281_127", "The original name was not parsed: #{snp.gene}")
    assert_equal(snp.contig, "IWGSC_CSS_1AL_scaff_1455974")
    assert_equal(snp.chromosome, "1A", "The chromosome wasnt parsed: #{snp.chromosome}")
    assert_equal(snp.position, 127, "The position is not parsed: #{snp.position}")
    #snp.setTemplateFromFastaFile(fasta_reference_db, flanking_size: 100)
    region = fasta_reference_db.index.region_for_entry(snp.contig).get_full_region
    snp.full_sequence = fasta_reference_db.fetch_sequence(region)

    assert_equal(snp.template_sequence, "actcgatcgtcagcacccgctggaacttggggaacgtcttgaacgccgcaagcaccggggcgtcctctgactgtatgagcacgcgctgcttacaggtctcYttgtcgtacccggacttgacaagcgctttggagaccgcatccaccacgtcaaggcttctggctataaggtacgtagcatgctgcactcggtaggtacaag".upcase)
    assert_equal(snp.sequence_original, "actcgatcgtcagcacccgctggaacttggggaacgtcttgaacgccgcaagcaccggggcgtcctctgactgtatgagcacgcgctgcttacaggtctc[C/T]ttgtcgtacccggacttgacaagcgctttggagaccgcatccaccacgtcaaggcttctggctataaggtacgtagcatgctgcactcggtaggtacaag".upcase)
    assert_equal(snp.position, 101)
    assert_equal(snp.original, "C")
    assert_equal(snp.snp, "T")
  end

  def test_vcf_line
    ref=@data + "/IWGSC_CSS_1AL_scaff_1455974_aln_contigs.fa"
    fasta_reference_db = Bio::DB::Fasta::FastaFile.new({:fasta=>ref})
    
    fasta_reference_db.load_fai_entries 
    vcf="IWGSC_CSS_1AL_scaff_1455974	127	test_snp	C	T	135.03	.	"

    chr_arm_parser = Bio::PolyploidTools::ChromosomeArm.getArmSelection("embl");
    snp = Bio::PolyploidTools::SNP.parseVCF(vcf, chr_arm_parser: chr_arm_parser)
    assert_equal(snp.gene , "test_snp", "The original name was not parsed: #{snp.gene}")
    assert_equal("1A", snp.chromosome, "The chromosome wasnt parsed: #{snp.chromosome}")
    assert_equal(127, snp.position, "The position is not parsed: #{snp.position}")
    snp.setTemplateFromFastaFile(fasta_reference_db, flanking_size:  100)
    assert_equal("actcgatcgtcagcacccgctggaacttggggaacgtcttgaacgccgcaagcaccggggcgtcctctgactgtatgagcacgcgctgcttacaggtctcYttgtcgtacccggacttgacaagcgctttggagaccgcatccaccacgtcaaggcttctggctataaggtacgtagcatgctgcactcggtaggtacaaga",     snp.template_sequence)
    assert_equal("actcgatcgtcagcacccgctggaacttggggaacgtcttgaacgccgcaagcaccggggcgtcctctgactgtatgagcacgcgctgcttacaggtctc[C/T]ttgtcgtacccggacttgacaagcgctttggagaccgcatccaccacgtcaaggcttctggctataaggtacgtagcatgctgcactcggtaggtacaag".upcase, snp.to_polymarker_sequence(100))
    assert_equal(101,snp.position)
    assert_equal("C",snp.original)
    assert_equal("T",snp.snp)

    vcf="IWGSC_CSS_1AL_scaff_1455974\t127\ttest_snp\tC\tT\t135.03\t.\tOR=reverse"

    chr_arm_parser = Bio::PolyploidTools::ChromosomeArm.getArmSelection("embl");
    snp = Bio::PolyploidTools::SNP.parseVCF(vcf, chr_arm_parser: chr_arm_parser)
    assert_equal(snp.gene , "test_snp", "The original name was not parsed: #{snp.gene}")
    assert_equal("1A", snp.chromosome, "The chromosome wasnt parsed: #{snp.chromosome}")
    assert_equal(127, snp.position, "The position is not parsed: #{snp.position}")
    snp.setTemplateFromFastaFile(fasta_reference_db, flanking_size:  100)
    assert_equal("actcgatcgtcagcacccgctggaacttggggaacgtcttgaacgccgcaagcaccggggcgtcctctgactgtatgagcacgcgctgcttacaggtctcYttgtcgtacccggacttgacaagcgctttggagaccgcatccaccacgtcaaggcttctggctataaggtacgtagcatgctgcactcggtaggtacaaga",     snp.template_sequence)
    assert_equal("TCTTGTACCTACCGAGTGCAGCATGCTACGTACCTTATAGCCAGAAGCCTTGACGTGGTGGATGCGGTCTCCAAAGCGCTTGTCAAGTCCGGGTACGACAA[G/A]GAGACCTGTAAGCAGCGCGTGCTCATACAGTCAGAGGACGCCCCGGTGCTTGCGGCGTTCAAGACGTTCCCCAAGTTCCAGCGGGTGCTGACGATCGAG", snp.to_polymarker_sequence(100))
    assert_equal(101,snp.position)
    assert_equal("C",snp.original)
    assert_equal("T",snp.snp)

  end

  def test_reference_snp

    ref=@data + "/IWGSC_CSS_1AL_scaff_1455974_aln_contigs.fa"
    fasta_reference_db = Bio::DB::Fasta::FastaFile.new({:fasta=>ref})
    
    fasta_reference_db.load_fai_entries 
    
    snp = Bio::PolyploidTools::SNP.parse("IWGSC_CSS_1AL_scaff_1455974,C,127,T,1A")
    assert_equal(snp.gene , "IWGSC_CSS_1AL_scaff_1455974", "The original name was not parsed: #{snp.gene}")
    assert_equal("1A", snp.chromosome, "The chromosome wasnt parsed: #{snp.chromosome}")
    assert_equal(127, snp.position, "The position is not parsed: #{snp.position}")
    snp.setTemplateFromFastaFile(fasta_reference_db, flanking_size: 100)
    assert_equal("actcgatcgtcagcacccgctggaacttggggaacgtcttgaacgccgcaagcaccggggcgtcctctgactgtatgagcacgcgctgcttacaggtctcYttgtcgtacccggacttgacaagcgctttggagaccgcatccaccacgtcaaggcttctggctataaggtacgtagcatgctgcactcggtaggtacaaga",     snp.template_sequence)
    assert_equal("actcgatcgtcagcacccgctggaacttggggaacgtcttgaacgccgcaagcaccggggcgtcctctgactgtatgagcacgcgctgcttacaggtctc[C/T]ttgtcgtacccggacttgacaagcgctttggagaccgcatccaccacgtcaaggcttctggctataaggtacgtagcatgctgcactcggtaggtacaag".upcase, snp.to_polymarker_sequence(100))
    assert_equal(101,snp.position)
    assert_equal("C",snp.original)
    assert_equal("T",snp.snp)
    
    flanking_size = 3

    snp = Bio::PolyploidTools::SNP.parse("IWGSC_CSS_1DL_scaff_2258883,A,12498,C,1D")
    snp.setTemplateFromFastaFile(fasta_reference_db, flanking_size: flanking_size)
    assert_equal(4,snp.position)
    assert_equal("A",snp.original)
    assert_equal("C",snp.snp)
    assert_equal("gatM", snp.template_sequence)

    snp = Bio::PolyploidTools::SNP.parse("IWGSC_CSS_1BL_scaff_3810460,G,1,T,1B")
    snp.setTemplateFromFastaFile(fasta_reference_db, flanking_size: flanking_size)
    assert_equal(1,snp.position)
    assert_equal("G",snp.original)
    assert_equal("T",snp.snp)
    assert_equal("Kaatt", snp.template_sequence)    
  end

  
end
