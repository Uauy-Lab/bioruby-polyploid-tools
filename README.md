#bio-polyploid-tools

##Introduction
This tools are designed to deal with polyploid wheat. The first tool is to design KASP primers, 
making them as specific as possible. 


##Installation
```sh
gem install bio-polyploid-tools
```
You need to have in your ```$PATH``` the following programs:

* [MAFFT](http://mafft.cbrc.jp/alignment/software/)
* [primer3](http://primer3.sourceforge.net/releases.php)
* [exonerate](http://www.ebi.ac.uk/~guy/exonerate/) 

The code has been developed on ruby 2.1.0, but it should work on 1.9.3 and above. 


#PolyMarker

To run poolymerker with the CSS wheat contigs, you need to unzip the 
(reference file)[ftp://ftp.ensemblgenomes.org/pub/release-25/plants/fasta/triticum_aestivum/dna/Triticum_aestivum.IWGSC2.25.dna.genome.fa.gz
].



```sh
polymarker.rb --contigs Triticum_aestivum.IWGSC2.25.dna.genome.fa --marker_list snp_list.csv --output output_folder
```

The snp_list file must follow the convention
<ID>,<Chromosome>,<SEQUENCE>
with the SNP inside the sequence in the format [A/T]. As a reference, look at test/data/short_primer_design_test.csv

If you want to use the web interface, visit the [PolyMarker webservice at TGAC](http://polymarker.tgac.ac.uk)

##Release Notes

###0.6.1


* polymarker.rb now validates that all the files exist. 
* BUGFIX: A reference was required even when it was not used to generate contigs. 

#Notes


* If the SNP is in a gap in the alignment to the chromosomes, it is ignored. 

BUG: Blocks with NNNs are picked and treated as semi-specific. 
BUG: If the name of the reference have space, the ID is not chopped. ">gene_1 (G12A)" shouls be treated as ">gene_1". 
TODO: Add a parameter file to configure the alignments. 
TODO: Produce primers for products of different sizes




