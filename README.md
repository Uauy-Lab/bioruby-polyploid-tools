#bio-polyploid-tools

##Introduction
This tools are designed to deal with polyploid wheat. The first tool is to design KASP primers, making them as specific as possible. 


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

To run poolymerker with the CSS wheat contigs, you need to unzip the reference file from [ensembl](ftp://ftp.ensemblgenomes.org/pub/release-25/plants/fasta/triticum_aestivum/dna/Triticum_aestivum.IWGSC2.25.dna.genome.fa.gz).



```sh
polymarker.rb --contigs Triticum_aestivum.IWGSC2.25.dna.genome.fa --marker_list snp_list.csv --output output_folder
```

The ```snp_list``` file must follow the convention ```ID,Chromosome,SEQUENCE``` with the SNP inside the sequence in the format [A/T]. As a reference, look at test/data/short_primer_design_test.csv

If you want to use the web interface, visit the [PolyMarker webservice at TGAC](http://polymarker.tgac.ac.uk)

The available command line arguments are: 

```
Usage: polymarker.rb [options]
    -c, --contigs FILE               File with contigs to use as database
    -m, --marker_list FILE           File with the list of markers to search from
    -g, --genomes_count INT          Number of genomes (default 3, for hexaploid)
    -s, --snp_list FILE              File with the list of snps to search from, requires --reference to get the sequence using a position
    -t, --mutant_list FILE           File with the list of positions with mutation and the mutation line.
    requires --reference to get the sequence using a position
    -r, --reference FILE             Fasta file with the sequence for the markers (to complement --snp_list)
    -o, --output FOLDER              Output folder
    -e, --exonerate_model MODEL      Model to be used in exonerate to search for the contigs
    -a, --arm_selection arm_selection_embl|arm_selection_morex|arm_selection_first_two
                    Function to decide the chromome arm
    -p, --primer_3_preferences FILE  file with preferences to be sent to primer3
    -v, --variation_free_region INT  If present, avoid generating the common primer if there are homoeologous SNPs within the specified distance (not tested)
    -x, --extract_found_contigs      If present, save in a separate file the contigs with matches. Useful to debug.
    -P, --primers_to_order			 If present, saves a file named primers_to_order which contains the KASP tails
```

###Custom reference sequences. 
By default, the contigs and pseudomolecules from [ensembl](ftp://ftp.ensemblgenomes.org/pub/release-25/plants/fasta/triticum_aestivum/dna/Triticum_aestivum.IWGSC2.25.dna.genome.fa.gz
) are used. However, it is possible to use a custom reference. To define the chromosome where each contig belongs the argument ```arm_selection``` is used.  The defailt uses ids like: ```IWGSC_CSS_1AL_scaff_110```, where the third field, separated by underscores is used. A simple way to add costum references is to rename the fasta file to follow that convention. Another way is to use the option ```--arm_selection arm_selection_first_two```, where only the first two characters in each contig is used as identifier, useful when pseudomolecules are named after the chromosomes (ie: ">1A" in the fasta file). 

If your contigs follow a different convention, in the file ```polymarker.rb``` it is possible to define new parsers, by adding at the begining, with the rest of the parsers a new lambda like:

```rb
arm_selection_functions[:arm_selection_embl] = lambda do | contig_name|
  arr = contig_name.split('_')
  ret = "U"
  ret = arr[2][0,2] if arr.size >= 3
  ret = "3B" if arr.size == 2 and arr[0] == "v443"
  ret = arr[0][0,2] if arr.size == 1   
  return ret
end
```

The function should return a 2 character string, when the first is the chromosome number and the second the chromosome group. The symbol in the hash is the name to be used in the argument ```--arm_selection```.  If you want your parser to be added to the distribution, feel free to fork and make a pull request. 



##Release Notes

###0.7.1
* BUGFIX: Now the parser for ```arm_selection_embl``` works with the mixture of contigs and pseudomolecules
*  DOC: Added documentation on how to use custom references.  

###0.7.0
* Added flag ```gebines_count``` for number of genomes, to be used on tetraploids, etc. 

###0.6.1


* polymarker.rb now validates that all the files exist. 
* BUGFIX: A reference was required even when it was not used to generate contigs. 

#Notes


* If the SNP is in a gap in the alignment to the chromosomes, it is ignored. 

* BUG: Blocks with NNNs are picked and treated as semi-specific. 
* BUG: If the name of the reference have space, the ID is not chopped. ">gene_1 (G12A)" shouls be treated as ">gene_1". 
* TODO: Add a parameter file to configure the alignments. 
* TODO: Produce primers for products of different sizes. This can probably be done with the primer_3_preferences option, but hasn't been tested. 




