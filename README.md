bio-polyploid-tools
===================

Introduction
-------------
This tools are designed to deal with polyploid wheat. The first tool is to design KASP primers, 
making them as specific as possible. 


Installation
------------
'gem install bio-polyploid-tools'

You need to have in your $PATH the following programs:
* [MAFFT]{http://mafft.cbrc.jp/alignment/software/}
* [primer3]{http://primer3.sourceforge.net/releases.php}
* [exonerate]{http://www.ebi.ac.uk/~guy/exonerate/} 

The code has been developed on ruby 2.1.0, but it should work on 1.9.3 and above. 


Polymarker
----------

To run poolymerker with the CSS wheat contigs, you need to unzip the 
reference file [Triticum_aestivum.IWGSP1.22.dna_rm.genome.fa.gz{ftp://ftp.ensemblgenomes.org/pub/release-22/plants/fasta/triticum_aestivum/dna/}. 

polymarker.rb --contigs Triticum_aestivum.IWGSP1.22.dna_rm.genome.fa --marker_list snp_list.csv --output output_folder

The snp_list file must follow the convention
<ID>,<Chromosome>,<SEQUENCE>
with the SNP inside the sequence in the format [A/T]. As a reference, look at test/data/short_primer_design_test.csv

Notes
-----

* If the SNP is in a gap in the alignment to the chromosomes, it is ignored. 

BUG: Blocks with NNNs are picked and treated as semi-specific. 
BUG: If the name of the reference have space, the ID is not chopped. ">gene_1 (G12A)" shouls be treated as ">gene_1". 
TODO: If reading from a reference file, only get one reference to align when the region is queried several times
TODO: Add a parameter file to configure the alignments. 
TODO: Produce primers for products of different sizes


