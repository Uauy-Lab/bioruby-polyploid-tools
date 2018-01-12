$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
$: << File.expand_path('.')
path= File.expand_path(File.dirname(__FILE__) + '/../lib/bioruby-polyploid-tools.rb')

#puts path
require path
require "test/unit"

class TestBlast < Test::Unit::TestCase
  Query  = File.dirname(__FILE__) + '/data/' + "BS00068396_51.fa"
  Target = File.dirname(__FILE__) + '/data/' + "BS00068396_51_contigs.fa"
  Blast_file = File.dirname(__FILE__) + '/data/' + "BS00068396_51_blast.tab"
  #Set up the paths
  def setup
    File.expand_path(File.dirname(__FILE__) + '/data/')
  end
  
  
  
  def test_blast_to_exo
    lines = File.readlines(Blast_file)
    expected = [
      "BS00068396_51 0 101 + 2AS_5222932 3015 2914 - 99",
    "BS00068396_51 0 101 + 2DS_5334799 6812 6913 + 99",
    "BS00068396_51 0 101 + 2BS_5245544 4549 4651 + 87",
    "BS00068396_51 101 0 - 2BS_5163353 7425 7323 - 87"]

    lines.each_with_index do |line , i|
      tmp =  Bio::DB::Blast.to_sugar(line)
      assert(tmp==expected[i])
    end
    
  end
  
end