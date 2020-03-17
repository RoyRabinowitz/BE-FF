
<!-- PROJECT LOGO -->
<br />
<p align="center">

  <h1 align="center">BE-FF</h1>

  <p align="center">
    Web-based tool which identifies suitable base editors to correct single nucleotide variations.
    <br />
    <a href="http://danioffenlab.pythonanywhere.com/"><strong>Go to website >> </strong></a>
   
  </p>




<!-- TABLE OF CONTENTS -->
## Table of Contents

* [About](#ABOUT)
* [BEsingle.py](#built-with)
  * [Methods to Enter Data](#methods)
* [BEMain.py](#BEMain.py)
* [Contact](#contact)

<!-- ABOUT  -->
## ABOUT

BE-FF is a web-based tool that receives SNV data and matches suitable BEs to correct the
variation. The code for the online tool is available here, as well as batch-mode code that may be used
to generate results for a large number of SNV's.
The 2 main scripts are:


## 1. BEsingle.py: 
This is the main script of the website. It receives a single SNV, via 3 possible methods: <br> 

* <i> Manually entered by user:</i> <br> Here you must enter a 51-nt long DNA sequence, 25 nt upstream to the mutation and 25 nt downstream, as well as the variation and the reading frame. 
![method1](method1.PNG)
 
* <i> Fetched by given rsID: </i> Enters a known rsID. You will then be presented with a table containing all the possible variations. Once you selects one of the options, it will automatically be inserted into the format in (a). 
![method2](method2.PNG)

* <i> Fetche by genomic coordinates: </i> You may select a genome, chromosome number, mutation position, variation, and reading frame, which will be inserted into the format in (a).
![method3](method3.PNG) 

**Note: Only one of these methods is required each time. **
 
You may also use the **'Advanced Options'** button and create you own BE. This base editor will appear as
'User customized BE' in the final results

Press the <u> submit </u> button to show the results. You will be presented with two tables:
1. Table of  base editors that will correct this variance and no other bases around it.
2. Table of base editors that will correctly edit the SNP. However it may also edit other flanking nucleotides.
Such changes are synonymous substitutions and the resulted amino acid sequence will probably be as good as the reference sequence.

## 2. BEMain.py
In the main function of this script, enter a csv table of the following format:
[Template](sample3.csv)

The result will be a new CSV file containing all the SNV's and the possible BE that can corect the variation. 



<!-- GETTING STARTED -->
## Getting Started

To get a local copy up and running follow these simple steps.


### Installation
 
1. Clone the repo
```sh
git clone https://github.com/RoyRabinowitz/BE-FF
```
2. Install biopython packages
```sh
pip install biopython 
```


## Main Functions Details

 ### `matchBE`(snp, BElist)

Recieves an snp object and returns a list of all matching BE that correct the mutation,
as well as a dictionary with all the possible locations (may be a few locations per BE) 

#### Parameters:

*   **snp** – snp object, has attributes such as snp.mutation, snp.wildtype, snp.upstream_sequence, snp.downstream_sequence
*   **BElist** – dictionary of the shape   
{BE: [PAM,activation window start,activation window end, mutation,variation, upstream/downstream]} 

 ### `cleanMatch`(snp,Matches, BElist,rev)
 Finds BE that will result in a prefect correction or a synonymous one (i.e. same smino acid)

#### Parameters:
*   **snp** - snp object
*   **Matches** - list recieved from previous function
*   **BElist** - dictionary of BE
*   **rev** - bool. whether the PAM is found on the given DNA strand or on its reverse complement.  

### `getRevComp`(snp)
Finds reverse complement of the SNP object, including changing its attirbutes respectively 

#### Parameters:
*   **snp** - snp object

### `SpecialCleanMatch`(snp, Matches, BElist, rev, orig_protein)
Finds BE that may correct the mutation by changing a base that is not the mutation. 
when the final amino acid is the same as the original

#### Parameters:
*   **snp** - snp object
*   **Matches** - list received from previous function
*   **BElist** - dictionary of BE
*   **rev** - bool.
*   **orig_protein** - the original amino acid sequence, to be compared with the final one 

## API

### `fetch_dna_coordinates`(genome, chromosome, startpos, endpos, cache_dir)
This function fetches sequence data, by inserting the parameters into the following url:
"http://genome.ucsc.edu/cgi-bin/das/{0}/dna?segment={1}:{2},{3}"
0 - genome
1 - chromosome
2 - start position
3 - end position

#### Parameters:
*   **genome, chromosome,start position,end position** - relevent sequence information
*   **chache_dir** - the local directory to which the data will be saved



`
## Edit the base editors list
To change the base editors used for analysis and their properties, open the [Scripts/baseEditorsTable.py](Scripts/baseEditorsTable.py) file
The base editors are categorized according to their type - CBEs and ABEs. 
Minor editing windows are indicated in a designated dictionary, following the main dictionary. 


<!-- CONTACT -->
## Contact

Roy Rabinowitz - royr2@mail.tau.ac.il

Shiri Almog - shirialmog1@gmail.com

Website: [http://danioffenlab.pythonanywhere.com](http://danioffenlab.pythonanywhere.com)

Project Link: [https://github.com/RoyRabinowitz/BE-FF](https://github.com/RoyRabinowitz/BE-FF)

Please visit our preprint for more information [https://www.biorxiv.org/content/10.1101/2020.01.06.890244v1](https://www.biorxiv.org/content/10.1101/2020.01.06.890244v1)


