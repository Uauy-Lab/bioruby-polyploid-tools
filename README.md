# bio-polyploid-tools

This project is now maintained by the informatics department at the John Innes Centre and can be found [here](https://github.com/Informatics-John-Innes-Centre/bio-polymarker).

## Introduction

This tools are designed to deal with polyploid wheat. The first tool is to design KASP primers, making them as specific as possible. 


## Installation

```sh
gem install bio-polyploid-tools
```
You need to have in your ```$PATH``` the following programs:

* [MAFFT](http://mafft.cbrc.jp/alignment/software/)
* [primer3](http://primer3.sourceforge.net/releases.php)
* [exonerate](http://www.ebi.ac.uk/~guy/exonerate/) 
* [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE%3DBlastDocs&DOC_TYPE%3DDownload)

The code was originally developed on ruby 2.1, 2.3 and 2.5. It may work on older version. However, it is only actively tested in currently supported ruby versions: 
  
  * 2.1.10
  * 2.2.5
  * 2.3.5
  * 2.4.2
  * 2.5.0

# PolyMarker

To run PolyMarker with the CSS wheat contigs, you need to unzip the reference file from  [ensembl](http://ftp.ensemblgenomes.org/pub/release-25/plants/fasta/triticum_aestivum/dna/Triticum_aestivum.IWGSC2.25.dna.genome.fa.gz).


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
    -i, --min_identity INT           Minimum identity to consider a hit (default 90)
    -a, --arm_selection arm_selection_embl|arm_selection_morex|arm_selection_first_two
                    Function to decide the chromome arm
    -p, --primer_3_preferences FILE  file with preferences to be sent to primer3
    -v, --variation_free_region INT  If present, avoid generating the common primer if there are homoeologous SNPs within the specified distance (not tested)
    -x, --extract_found_contigs      If present, save in a separate file the contigs with matches. Useful to debug.
    -P, --primers_to_order			 If present, saves a file named primers_to_order which contains the KASP tails
```

## Input formats

The following formats are used to define the marker sequences:

### Marker list

If the option ```--marker_list FILE``` is used, the SNP and the flanking sequence is included in the file. The format contains 3 columns (the order is important):

* **snp_name** The ID of the marker. Must be unique. 
* **target chromosome** for the specific primers. Must be in line with the chromosome selection critieria. 
* **sequence** The sequence flanking the SNP with the SNP highligted on square brackets (```[]```) and the two alleles separated by a forward slash (```/```). 

#### Example:

```
BS00068396_51,2A,CGAAGCGATCCTACTACATTGCGTTCCTTTCCCACTCCCAGGTCCCCCTA[T/C]ATGCAGGATCTTGATTAGTCGTGTGAACAACTGAAATTTGAGCGCCACAA
```

### SNP list

If the flanking sequence is unknow, but the position on a reference is available,  the option ```--snp_list``` can be used and the FASTA file with the reference sequence must be provided with the option ```--reference```. This is to allow the use of a different assembly or set of contigs used for the discovery of the SNPs that are different to the reference given in the option ```--contigs```. The format contains the following positional columns: 

* **scaffold** The sacffold where the SNP is. 
* **reference allele** The base in the reference (may or may not be the same as in the reference file.
* **position** Position of the SNP. The first base in the scaffold is base 1. 
* **alternative allele** The base in the alternative allele. 
* **target chromosome** for the specific primers. Must be in line with the chromosome selection critieria. 

#### Example

```
IWGSC_CSS_1AL_scaff_110,C,519,A,2A
```

This file format can be used with ```snp_positions_to_polymarker.rb``` to produce the input for the option```--marker_list```.


### Custom reference sequences. 

By default, the contigs and pseudomolecules from [ensembl](ftp://ftp.ensemblgenomes.org/pub/release-25/plants/fasta/triticum_aestivum/dna/Triticum_aestivum.IWGSC2.25.dna.genome.fa.gz
) are used. However, it is possible to use a custom reference. To define the chromosome where each contig belongs the argument ```arm_selection``` is used.  The defailt uses ids like: ```IWGSC_CSS_1AL_scaff_110```, where the third field, separated by underscores is used. A simple way to add costum references is to rename the fasta file to follow that convention. Another way is to use the option ```--arm_selection arm_selection_first_two```, where only the first two characters in each contig is used as identifier, useful when pseudomolecules are named after the chromosomes (ie: ">1A" in the fasta file). 

If your contigs follow a different convention, in the file ```ChromosomeArm.rb``` it is possible to define new parsers, by adding at the begining, with the rest of the parsers a new lambda like:

```rb
@@arm_selection_functions[:embl] = lambda do | contig_name|
  arr = contig_name.split('_')
  ret = "U"
  ret = arr[2][0,2] if arr.size >= 3
  ret = "3B" if arr.size == 2 and arr[0] == "v443"
  ret = arr[0][0,2] if arr.size == 1   
  return ret
end
```

The function should return a 2 character string, when the first is the chromosome number and the second the chromosome group. The symbol in the hash is the name to be used in the argument ```--arm_selection```.  If you want your parser to be added to the distribution, feel free to fork and make a pull request.  

##Using blast

To use blast instead of exonerate, use the following command:

```
./bin/polymarker.rb --contigs test/data/BS00068396_51_contigs.fa --marker_list test/data/BS00068396_51_for_polymarker.fa  --aligner blast  -a arm_selection_first_two
```


## Release Notes

### 0.9.7 
There was some strange issue with the numbering, so bumped to 0.9.7

* Moved the arm selection function for fields in the chromosome name to the ```ChromosomeArm``` class.

### 0.8.7
* FEATURE: ```polymarker.rb``` now also prints the total number of hits found. 

### 0.8.6

* BUGFIX: ```priemr3.rb``` had a regression when adding the repetitive flag to the ```@values``` array. This lead to the wrong order of the columns in the output and possibly other secondary effects. 

### 0.8.5

* Added the option ```--max_hits``` to ```polyamarker.rb``` to set a maximum number of bast hits to identify repetitive regions. This adds the column ```is_repetitve``` to the output. The mask is not calculated in repetitive regions and the primers are designed as non-specific. 

### 0.8.4 

* Added script ```tag_stats.rb`` That gets the descriptive statistics for a tag in a bam file for each reference. 

```bash
ruby tag_stats.rb -b HI.3206.006.Index_2.CS_125RNA_14d_Leaf8.sorted.bam -r /Users/ramirezr/Dropbox/JIC/expVIPMetadatas/RefSeq1.0/Genes/annotation/IWGSCv1.0_UTR_ALL.cdnas.fasta --tag 'NH'
```

### 0.8.3

* BUGFIX: ```ChromosomeArm.rb``` was fixed to conform the module assumptions for the package. 


### 0.8.2

* FEATURE: The functions to select the chromosome arm are now in ```lib/bio/PolyploidTools/ChromosomeArm.rb``` and the help message is updated automatically with the valid options. 
* FEATURE: Added option ```filter_best``` to replicate the original behaviour of selecting the best hit of each chromosome. Still useful for assemblies which still contain synthetic duplications.

### 0.8.1

* BUGFIX: There was an error which prevented the correct localisation of the SNP in markeres with gaps in the local alignment before the position with the snp.
* FEATURE: PolyMarker now selects the best hit of the target chromosome. This improves the specificity in regions with a recent duplication. The drawback is that if your assembly has artificial repetitions, the primers won't be marked as 'chromosome specific', but as 'chromosome semi-specific '. In a future version this will be addressed. 

### 0.8

* FEATURE: ```polymarker.rb``` added the flag ```--aligner blast|exonerate ``` which lets you pick between ```blast``` or ```exonerate``` as the aligner. For blast the default is to have the database with the same name as the ```--contigs``` file. However, it is possible to use a different name vua the option ```--database```.

### 0.7.3

* FEATURE: ```polymarker.rb``` Added to the flag ```--arm_selection``` the option ```scaffold```, which now supports a scaffold specific primer.
* FEATURE: ```snp_position_to_polymarker``` Added the option ```--mutant_list```  to prepare files for PolyMarker from files with the following columns ```ID,Allele_1,position,Allele_1,target_chromosome```. 

### 0.7.2

* FEATURE: Added a flag ```min_identity``` to set the minimum identity to consider a hit. The default is 90

### 0.7.1
* BUGFIX: Now the parser for ```arm_selection_embl``` works with the mixture of contigs and pseudomolecules
*  DOC: Added documentation on how to use custom references.  

### 0.7.0
* Added flag ```genomes_count``` for number of genomes, to be used on tetraploids, etc. 

### 0.6.1


* polymarker.rb now validates that all the files exist. 
* BUGFIX: A reference was required even when it was not used to generate contigs. 

# Notes

* BUG: Blocks with NNNs are picked and treated as semi-specific. 
* BUG: If the name of the reference have space, the ID is not chopped. ">gene_1 (G12A)" shouls be treated as ">gene_1". 
* TODO: Add a parameter file to configure the alignments. 
* TODO: Produce primers for products of different sizes. This can probably be done with the primer_3_preferences option, but hasn't been tested. 




