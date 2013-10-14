$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
$: << File.expand_path('.')
path= File.expand_path(File.dirname(__FILE__) + '/../lib/bioruby-polyploid-tools.rb')

#puts path
require path
require "test/unit"

class TestExonContainer < Test::Unit::TestCase
  Query=File.dirname(__FILE__) + '/data/'+"BS00068396_51.fa"
  Target=File.dirname(__FILE__) + '/data/'+"BS00068396_51_contigs.fa"
  
  
   def test_simple_container_test
     
   end
end