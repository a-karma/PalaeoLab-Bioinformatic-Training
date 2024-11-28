![bio_logo](../IM/header.png)
### 2. Working with bioinformatic softwares using conda
In Session 1 we have seen three examples of text file that are commonly used in bioinformatics:

- The `fasta` format to store DNA sequence information
- The `GTF` (Gene Transfer Format) developed for the Ensembl genome browser
- The `BED` (Browser Extensible Data) format developed at UCSC for the Genome Browser tool

The `GTF` and the `BED` format are both TAB-separated files used to store genomic regions as coordinates along with their associated annotations (see the previous session for more information on `GTF`). While `GTF` is used to store mostly information about genic structure accross the genome (i.e. the G stands for gene), the `BED` (short for Browser Extensible Data) stores information about any genomic features. It is simpler than `GTF` and widely used in bioinformatics and consist of at least 3 fields:

1. chromsome ID
2. start of feature (0 based coordinate; see link to UCSC below)
3. end of the feature 

Other fields (4 etc.) are free form (as long as they are tab separated) - these can be used to store any information you may want, e.g. the name of the interval. For more information about this widely use format refer to the UCSC website here: https://genome.ucsc.edu/FAQ/FAQformat.html#format1

Although in principles we could manually edit these files using standard text editors this becomes very unpractical when dealing with very large files.
A more efficient option would be combinig command line tools (like `sed`, `awk`, or `grep`) but performing complex tasks using only these tools is not straightforward and often require a good basis understanding of unix. Luckly for us, bioinformaticians have created various software specifically designed to manipulate these file formats. 

Let's have a look at `bedtools` a powerful toolset for genome arithmetic that focuses on `BED` format.

In your terminal, please type: 

```sh
bedtools --help
```

Looks like the program is not installed :( 

To protect the integrity of the file system on a server, normal users do not have permissions to directly install softwares. Moreover, almost any bioinformatic tool will rely on specific libraries or package versions that might create conflicts or even impeed the functionality of other programs. To circumvent these issues, we need to make sure that the software we need are installed in a "confined space" (a.k.a environment) containing all the necessary dependencies which becomes accessible only when we need it. This can be easily implemented using `conda` which is a package, dependency, and environment manager for any programming language. 

You can read more about conda [here](https://docs.conda.io/en/latest/).

We have already installed conda on our server and we have created a different environment for each day of the workshop. In order to have access to conda please run:

```sh
source /home/anaconda3/bin/activate
conda init
```
You need to run this command only once and you should see a change in your prompt: the word `(base)` appears on the left.
This is signalling that you are now in the conda base environment which represent the default space. 
To activate the environment for this session, simply run:

```sh
conda activate Day_1
```

> Question: Have a look at your prompt again, what do you see? 

After activating an environment all software installed in it become immediately accessible. Let's check whether we can use bedtools now:

```sh
bedtools --help
```

Hurray! Now that bedtools is accessible let's see what we can do with it.

The `snp_ch30.bed` file in the folder `/home/DATA/module_1/` is an example of a "customized" bed format. It contains the three mandatory fields (chromosome, start, end) plus an unusual 4th field. In that column I have stored the genotype of 4 individuals at that position. If the 4 th column has a lowercase `m`, that particular site is monomorphic.
You can have a look at it using `less` or inspect just three lines with a combination of head and tail, see for example what you get by running:

```
head -65 ~/session2/raw_data/Day_1/snp_ch30.bed | tail -3
```

Suppose you are interested in analysing neutral evolving sites, therefore, you may want to remove from the analysis all sites that are likely to be under selective pressures. 
As a first approximation, we could take a conservative approach and start to analyse polymorphic sites (SNPs) that are not in coding region (i.e. do not code for protein and so less likely to be under selective constraints; please speak to your instructor if you do not understand this concept). Performing this task manually is obviously tedious and very time consuming but it's super fast using a software like bedtools:

```sh
cd ~/session2/
bedtools intersect -a ./raw_data/Day_1/snp_ch30.bed -b ./raw_data/Day_1/genes_chr30.gtf -v > ./results/snp_filtered.bed
```
The intersect command reports overlapping regions between two BED/GFF/GTF files by comparing the coordinates of the genomic feature listed in them.
The `-v` flag tells `intersect` to report all lines in file A (specified using the `-a` flag) that DO NOT overlap with the genomic intervals listed in file B (-b flag).

> Exercise 2
>
> How many SNPs we have excluded?
>
> hint: remember that each SNP information is recorded on a single line in the bed file format

There is a lot more that you can do with genome arithmetic. Let’s picture a more complex scenario: suppose you are interested in studying the promoter regions of various genes on this chromosome. You have a fasta file with the entire chromosome sequence (see `ptw_ch30.fa` file) and you would like to examine 10 kb upstream the starting codon of each gene. 

How can you get that information?

First of all you need to get the coordinates of the starting codons. This is very easy using grep:
```sh
grep 'start_codon' ./raw_data/Day_1/genes_chr30.gtf > ./raw_data/CDS_start.gtf
```

Now you need to modify this file using the function `flank` implemented in bedtools which will create flanking intervals for each region in a BED/GFF/VCF file.
This function requires also a genome file defining the length of each chromosome, so let’s create this file first.
```sh
echo -e "chr30\t40214260" > ./raw_data/ch30_length.bed
```
The above command `echo` would normally output on screen whatever string you type after it. The value `40214260` is the size of chromosome 30 in the dog genome (which we used for this example).
The `-e` option tells the software to enable the interpretation of backslash escapes and `\t` stands for TAB.
Now we can run:

```sh
bedtools flank -i ./raw_data/CDS_start.gtf -g ./raw_data/ch30_length.bed -l 10000 -r 0 > ./raw_data/ch30_promoters.gtf
```
where the `-l` and `-r` flags stand for left (or upstream) and right (or downstream) and the number following each of these flags represents the length of the flanking region.

Finally we can use the `getfasta` function in bedtools to extract the actual sequences of the promoter regions:

```sh
bedtools getfasta -fi ./raw_data/Day_1/ptw_ch30.fa -bed ./raw_data/ch30_promoters.gtf -fo ./results/ptw_prom_sequences.fa
```

where the `-fi` flag stands for file input, the `-bed` indicates the coordinate file while the `-fo` option stands for file output. 

You can examine the first output line using: 
```sh
cat ./results/ptw_prom_sequences.fa | head -1
```

> Exercise 3
>
> Combine the intersect and flank functions in order to filter the `snp_ch30.bed` file
> by excluding CDS and all regions that are 5 kb from the start and the stop codon of each CDS.


If you finished early and you want to keep playing with bedtools have a look at this tutorial here: https://sandbox.bio/tutorials?id=bedtools-intro
