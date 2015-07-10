Analyzing & Tagging GoNL Deletions 
==================================

This repository contains the code/scripts for analyzing and tagging deletions for the [Genome of the Netherlands](http://www.nlgenome.nl/). The code serves two purposes: 

* explore the extent in which linkage disequilibria (correlations) between GoNL deletions and GWAS SNPs occur, and 

* finding appropriate tag SNPs for the found deletions. 

The code is written in both C and Python. 

***

## Installation 

### Dependencies 

The compilation of the C code requires the following libraries to be installed:

* The _GNU scientific library_ (GSL - see http://www.gnu.org/software/gsl)

* _CMake_ (see http://www.cmake.org)

The project depends for Python on the following packages: 

* _PyVCF_ (see https://github.com/jamescasbon/PyVCF) for working with VCF files

* _snakemake_ (see https://bitbucket.org/johanneskoester/snakemake/wiki/Home) for using the pipeline, see the file `Snakefile' in the main directory 

### Installation instructions 

In order to compile the C code in the folder `src/`, type in the main directory: 

```
	$ cmake . 
	$ make
	$ make install 
```

The executables `gonl_create_pairs`, `gonl_imputation` and `gonl_tag_deletions` will be placed in the `bin/` folder together with the Python scripts. 

*** 

## Directory structure 

The repository consists of the following directories: 

* `backup/` - contains some older versions of the project. 

* `bin/` - contains the Python scripts and executables used for calling the somatic mutations. 

* `data/` - contains a few of the raw data files used in the project. Some of the data files are ignored (e.g., in the folders `gonl-deletions` and `gonl-snps`) due to their size. Given the original VCF files all the data can be reproduced.

* `include/` - contains the header files for the C-code.

* `matlab/` - contains some Matlab code for plotting the results. 

* `paper/` - some texts and plots for the GoNL-SV paper. 

* `presentations/` - some presentations for the TC meetings on the results. 

* `plots/` - LaTeX files for the figures.

* `results/` - contains various (intermediate) results. The actual results are not in the repository (due to their size) but can be reproduced with the original VCF files and the present scripts/code. 

* `src/` - contains the C-code. 

***

## File formats 

This section contains a description of the two (novel) file formats used in this projects: `.raw-observations` and `.calls`. 

### .pairs

This file format is used for representing SNP-deletion pairs. The first line always represents a SNP; the lines that follow contain the data on deletions in the vicinity of the SNP. The number of deletions that follow the SNP may vary. The data on the SNP is structured as (space-separated):

	<chr> <type> <position> <ref> <alt> <hit-allele> <dist-tss> <region>

where `<chr>` is the chromosome on which the SNP resides. `<type>` can either be `*` in case of a regular SNP and `!` in case of a known GWAS SNP. `<position>` is the position of the SNP as given in the VCF file. `<ref>` and `<alt>` are the reference and the alternative alleles. In case that the SNP is a GWAS SNP, i.e., `<type>` is equal to `!`, then `<hit-allele>` is the allele associated with the disease. In case of a regular SNP, this field simply contains a `.`. The field `<dist-tss>` denotes the distance (in bp) to the closest transcription start site. The field `<region>` can either be `intronic`, `exonic` or `intergenic`, dependent on the location of the SNP. 

The lines that follow the SNP are the deletions. These lines always start with `-`. A line like this is structured as follows: 

	- <pos> <length> <R> <p> <A> <B> <C> <D>

where `<pos>` is the position of the deletion and `<length>` is its length. The last four columns denote the following 2x2 contigency table: 

|               | reference allele | alternative allele | _total_ |
| ------------- |:----------------:|:------------------:|:-------:|
| deletion      | A                | C                  | A + C   |
| no deletion   | B                | D                  | B + D   |
| _total_       | A + B            | C + D              | A+B+C+D | 

`<R>` is the Pearson _R_ and `<p>` is the p-value found when applying Fisher's two-sided exact test on the table. 

### .tagsnps

Every line represents one deletion-tag snp pair and consists of 14 columns in total (space-delimited): 

1. the chromosome on which the deletion resides. 

2. the deletions position as given in the VCF file. 

3. the length of the deletion. 

4. position of the tag SNP with the maximum R-squared. 

5. reference allele of the SNP.

6. alternative allele of the SNP.

7. Pearson R-value.

8. p-value (Fisher's two-sided exact test).

9. Conditional probability of having the deletion given the presence of the reference allele. 

10. Conditional probability of having the deletion given the presence of the alternative allele. 

The columns 11-14 provide the counts for the contingency table A, B, C and D. See the 2x2 contingency table of the previous section.

***

## Contact

Louis Dijkstra

__E-mail__: louisdijkstra (at) gmail.com 


