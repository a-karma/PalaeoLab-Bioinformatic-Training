![bio_logo](../IM/header.png)
## Module 1 - Basic concepts of command line programming - Session 2
In this session we are going to apply what we already know about variables, conditional and loops to explore and manipulate text files focusing on common bioinformatic file formats. Before we start, let's prepare the directory structure for this session. Please connect to our cloud server and run:

```sh
cd ~/module1
mkdir bio_formats
cd bio_formats
mkdir fasta; mkdir gtf; mkdir bed; mkdir raw_data
ln -s /home/DATA/module_1/ ~/module1/bio_formats/raw_data
```
The commands above should be quite straightforward now but if you don't understand them, please revise Session 1.

### 1. Regular expressions & file manipulations
In this first section we will mostly focus on regular expressions, and we will introduce three fundamental tools: `grep`, `sed` and `awk`.

Regular expressions, also known as regex or regexp, are a powerful tool for manipulating text. They allow you to search for, extract, and modify patterns in text data. In Bash scripting, regular expressions are particularly useful for processing text files and automating tasks that involve text manipulation.

#### The FASTA format
Let’s start by looking at the file called `random.fasta`:
```sh
cd raw_data
less random.fasta
Press Q to exit.
```

As the name suggests, this file contains some random DNA sequences of different length stored in a fasta format (a very common format for DNA analysis).

As you can see, each entry consists of two lines: a header (with the sequence identifier) and a second line containing the actual sequence.
Lets calculate how many lines are in this fasta file:
```sh
wc -l random.fasta
```
This gives us the total number of lines. Divide this number by two and you'll get the number of sequences.

Now, let’s look at the first 7 sequences which correspond to the first 14 lines of the file.
We can easily print them on screen using:
```sh
head -14 random.fasta
```
The `head` command is a fundamental tool in Bash scripting used to display the first portion of a file's contents. It's commonly used to quickly preview the beginning of a file or check for specific information at the start.

If instead we were interested in the last 4 entries, we would use:
```sh
tail -8 random.fasta
```
The `tail` command is a versatile tool in Bash scripting used to display the last portion of a file's contents. It's commonly used to quickly review the end of a file, check for recent changes, or monitor log files in real time.

We could use a combination of these two commands to extract a set of sequences in the middle of the file:
```sh
head -20 random.fasta | tail -4
```
With the last command we have selected the 9th and 10th entries corresponding to the lines from 17 to 20 in our fasta file.

The vertical bar (`|`) is called `pipe` and it is used to connect the two commands (`head` and `tail` in this case). 
Specifically, it redirects the standard output of the first command which then serves as input for the second command.

Standard output, also known as stdout, is the default output stream of a process in a Unix-like operating system. It is the channel through which a process (here `head`) sends its output to another process (here `tail`).

Piping the output of head into tail or vice versa is a simple way to extract a block of lines but it becomes very slow if the file you are dealing with is huge. An alternative is using `sed`. 

`sed` is a stream editor, a powerful tool for manipulating text data. It allows you to search, replace, insert, and delete text within files or data streams. sed is commonly used in Bash scripting for tasks like processing log files, cleaning up text files, and performing text-based transformations.

Let’s consider the following commands:

```sh
sed -n '9,12p' random.fasta
sed -n '9,+3p' random.fasta
```
These are alternative ways of printing a range of lines. In this case we are printing lines from 9 to 12 which of course correspond to our 5th and 6th entry.

> `Exercise 1`
>
> Use `sed` to extract the 23rd, 24th, and 25th sequence in the random.fasta file
>
> Redirect the output to a file called `23to25_seq.fasta` inside the `/bio_formats/fasta/` directory 

It is unlikely though that we will know in advance the line numbers of the entries that are relevant to our analysis. Most of the time will have to parse the file and look for patterns. That’s when regular expressions (`regex`) become very useful.

As you may have noticed, all header lines in random.fasta start with `>seq` followed by a number, `Hg`, and a letter, separated by underscores. The string `Hg` stands for haplogroup (A,B,or C) and we might be interested in knowing how many reads we have for each group. We can calculate this using `grep`.

`grep` is a powerful tool for searching for patterns in text files. It is commonly used in Bash scripting to locate specific text strings, analyze log files, and perform text-based searches.

`grep` stands for "global regular expression print" and it has the following syntax:
```sh
grep 'regex' target_file
```

Let's consider the following command as our first example:

```sh
grep '>seq.*_Hg_A' random.fasta | wc -l
```
In this example the regex or the pattern that we are looking for is `>seq.*_Hg_A`. This regular expression is designed not to match a unique sequence of characters but rather a series of sequences that possess similar features.

Let's examine the pattern in more details. The first part of our expression is `>seq.*`. This matches any line containing the string `>seq` followed by any character (represented by the `.` symbol) that appears zero or more times (`*`). The second part is `_Hg_A` which is designed to match any line containing this exact string. Therefore, our regular expression will match any line that:
- starts with the exact series of characters `>seq`
- is followed by any combination of characters repeated zero or more times
- is followed by the string `_Hg_A`.
Thus, `grep` will print all the header lines in the random.fasta corresponding to haplogroup A and this output is the piped (|) into the command `wc -l` which simply counts the number of matching lines.

> `Exercise 2`
>
> Find out how many sequences we have for each haplogroup (A,B,or C) by modifying the pattern of the grep command above.

We can also use a regex inside a sed command. For example, let's extract the 3rd sequence of each haplogroup:
```sh
sed -n '/seq_3_/,+1p' random.fasta > ../fasta/third_seq_all_Hg.fasta
```
As you can see, it looks very similar to the sed command we used before with the exception that instead of providing sed with specific line number, 
here we have specified a pattern (`/seq_3_/`) and asked the program to print each matching line plus and the following one: (+1p). 
Finally, we have redirected the output to store this information into a file called `third_seq_all_Hg.fasta`. Given that our current directory is `raw_data` we need to access the parent directory (`bio_format`) and from there we can then descend into the `fasta` directory where our new file belongs. Therefore, the path for redirection is `../fasta/` followed by the file name.

#### The GTF 
Let’s have a look at a different file format and keep  experimenting with regex. In the`/home/DATA/module_1` folder (which is linked to `~/module1/bio_formats/raw_data` folder) you should see a file called `dog_genes.gtf`.

A GTF (Gene Transfer Format) file is a text-based file format used to store information about the genomic structure of genes and their features. It is commonly used in genomics research to annotate and analyze genome sequences. It is a tab separated file containing annotations for coding sequences in the dog genome.

We can extract the header if of this file (lines starting with #) by running:
```sh
cd ~/module1/bio_formats/raw_data
grep '^#' dog_genes.gtf
```
> `Exercise 2`
>
> Remove the header of this file using the ”select non matching lines” option of grep.
>
> See `grep --help` to identify the correct flag for the job.
>
> Then redirect the output to a file called `dog_genes_no_H.tsv` inside the `bio_formats/gtf/` directory.

The first line of your file without a header should look like:
`X ensembl gene 1575 5716 . + . gene_id "ENSCAFG00000010935"; gene_version "3"; gene_source "ensembl" ... etc`



It contains a lot of information that is not relevant for us at the moment. 

The fields (or columns) that we are interested in are:
- The chromosome (X: 1st field)
- The type of the feature (gene: 3rd field)
- The starting position of the feature (1575: 4th field)
- The ending position of the feature (5716: 5th field)

> `Exercise 4`
> 
> Use `cut` to extract the required fields from `dog_genes_no_H.tsv` file.
>
> Then redirect the output to a file called `dog_genes_table.tsv` inside your `bio_formats/gtf/` directory.
>
> See `cut --help` to identify the option for fields

Now that we have extracted the relevant information, we would like to make a few adjustments to our table. Let’s start with adding the string `chr` at the beginning of each line.

We can do this easily by using the substitution command in `sed` which has the following general syntax:

```sh
sed 's/target/replacement/'
```
In our case we are going to modify the file in-place using the `-i` flag.

```sh
cd ~/module1/bio_formats/gtf/
sed -i 's/^/chr_/' dog_genes_table.tsv
```
In our example, the caret symbol (^) is a regex which denotes the beginning of a line and we replaced this with `chr_`. You can check whether the substitution worked or not by examining the first 10 lines of the table with head.

#### The BED format
The next thing we would like to do is switching the order of the columns in our table:
- 1. Chromosome
- 2. Starting Position
- 3. Ending Position
- 4. Feature Type

This requires a simple awk command:
```sh
awk 'BEGIN {OFS="\t"};{print $1,$3,$4,$2}' dog_genes_table.tsv > dog_genes_tab.bed
```
`awk` is a powerful scripting language designed for processing text files. It is commonly used in Bash scripting to manipulate, analyze, and transform text data. `awk` is particularly well-suited for tasks like parsing log files, extracting information from text files, and performing text-based calculations.

The OFS option before the print command stands for ”Output Filed Separator” and we set it to Tab (`\t`) to ensure our table has the correct delimiter for a ”tab separated file” (TSV). `awk` stores each field in a different variable which is accessible via the `$` symbol.  We made use of this feature to put the columns in the right order.

We are almost done with pre-processing our data but there’s still something that’s not quite right with it. Have a look at the first column:

```sh
cd bed
cat dog_genes_tab.bed | cut -f 1
```
Have you noticed that the chromosomes are not in the right order? Let’s fix it!
```sh
sort -V -o dog_genes_sorted_tab.bed dog_genes_tab.bed
```
The `sort` command is a fundamental tool in Bash scripting used to organize and rearrange data based on specified criteria. It is commonly used to `sort` text files numerically or alphabetically, making it a versatile tool for data manipulation and analysis.

Here we used the -V option because we are dealing with a mixture of numerical and string data. Note also the -o to specify the output file which must precede the input. Our table is now in the corrected format and it's ready to be analysed.

> `Exercise 5`
> 
> How many coding regions (CDS) on the X chromosome are listed in our bed file?
> 
> Use the commands you have learn to find out   

We will have more information about `BED` files in the next section.


### 2. Working with bioinformatic softwares using conda
In the previous section we have seen three examples of text file that are commonly used in bioinformatics:

- The `fasta` format to store DNA sequence information
- The `GTF` (Gene Transfer Format) developed for the Ensembl genome browser
- The `BED` (Browser Extensible Data) format developed at UCSC for the Genome Browser tool

The `GTF` and the `BED` format are both TAB-separated files used to store genomic regions as coordinates along with their associated annotations (see the previous session for more information on `GTF`). While `GTF` is used to store mostly information about genic structure accross the genome (i.e. the G stands for gene), the `BED` (short for Browser Extensible Data) stores information about any genomic features. It is simpler than `GTF` and widely used in bioinformatics and consist of at least 3 fields:

1. chromsome ID
2. start of feature (0 based coordinate; see link to UCSC below)
3. end of the feature 

Other fields (4th, 5th, etc.) are free form (as long as they are tab separated) - these can be used to store any information you may want, e.g. the name of the interval. For more information about this widely use format refer to the UCSC website here: https://genome.ucsc.edu/FAQ/FAQformat.html#format1

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

We have already installed conda on our server and we have created an environment with all software and dependencies need for this course. In order to have access to conda please run:

```sh
source /home/anaconda3/bin/activate
conda init
```
You need to run this command only once and you should see a change in your prompt: the word `(base)` appears on the left.
This is signalling that you are now in the conda base environment which represent the default space. 
To activate the environment for this session, simply run:

```sh
conda activate bio
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

> Exercise 6
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

> Exercise 7
>
> Combine the intersect and flank functions in order to filter the `snp_ch30.bed` file
> by excluding CDS and all regions that are 5 kb from the start and the stop codon of each CDS.


If you want to keep playing with bedtools have a look at this tutorial here: https://sandbox.bio/tutorials?id=bedtools-intro



 
