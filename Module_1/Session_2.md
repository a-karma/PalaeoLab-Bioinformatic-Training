![bio_logo](../IM/header.png)
## Module 1 - Basic concepts of command line programming - Session 2
### 1. Regular expressions & file manipulations
In this section we are going to apply what we already know about variables, conditional and loops to explore text files. 
We will mostly focus on file manipulation using regular expressions, using `grep`, `sed` and `awk`.

Regular expressions, also known as regex or regexp, are a powerful tool for manipulating text. They allow you to search for, extract, and modify patterns in text data. In Bash scripting, regular expressions are particularly useful for processing text files and automating tasks that involve text manipulation.

Let’s start with opening a new terminal, connecting to the server and look at the file called `random.fasta`:
```sh
cd /home/DATA/module_1/
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
With the last command we have selected the 9th and 10th entries corresponding to the lines from 17To assure yourself that all directories have been created you can use the command ls which lists the content of the current directory.

 to 20 in our fasta file.
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

It is unlikely though that we will know in advance the line numbers of the entries that are relevant to our analysis. Most of the time will have to parse the file and look for patterns. That’s when regular expressions (`regex`) become very useful.

As you may have noticed, all header lines in random.fasta start with `>seq` followed by a number, `Hg`, and a letter, separated by underscores. The string `Hg` stands for haplogroup (A,B,or C) and we might be interested in knowing how many reads we have for each group. We can calculate this using `grep`.

`grep` is a powerful tool for searching for patterns in text files. It is commonly used in Bash scripting to locate specific text strings, analyze log files, and perform text-based searches.

`grep` stands for "global regular expression print" and it has the following syntax
```sh
grep 'regex' target_file
```

Let's consider the following command:

```sh
grep '>seq.*_Hg_A' random.fasta | wc -l
```
Let's unpack this command:

In our example the regex or the pattern that we are looking is `>seq.*_Hg_A`. 

This regular expression is designed not to match a unique sequence of character but rather a series of sequences that possess similar features.
The regular expression `>seq.*` matches any line containing the string >seq followed by any character (represented by the `.` symbol) that appears zero or more times (`*`). This means that the regex will match any line that starts with the exact sequence `>seq` and is followed by any combination of characters followed by the string `_Hg_A`. Thus, `grep` will print all the header lines in the random.fasta corresponding to haplogroup A and this output is the piped (|) into the command `wc -l` which simply counts the number of matching lines.

> `Exercise 4`
>
> Find out how many sequences we have for each group by modifying the pattern of the grep command above.

We can also use a regex inside a sed command. For example, let's extract the 3rd sequence of each haplogroup:
```sh
sed -n '/seq_3_/,+1p' random.fasta > ~/project_bash/raw_data/third_seq_all_Hg.fasta
```
As you can see, it looks very similar to the sed command we used before with the exception that instead of providing sed with specific line number, 
here we have specified a pattern (`/seq_3_/`) and asked the program to print each matching line plus and the following one: (+1p). 
Finally, we have redirected the output to store this information into a file called third_seq_all_Hg.fasta.

> `Exercise 5`
> 
> Now navigate to your `raw_data` directory and visualise the content of the file on screen using the command `cat name-of-the-file`.

Let’s have a look at a different file format and keep  experimenting with regex. In the`/home/DATA/module_1` folder you should see a file called `dog_genes.gtf`.

A GTF (Gene Transfer Format) file is a text-based file format used to store information about the genomic structure of genes and their features. It is commonly used in genomics research to annotate and analyze genome sequences. It is a tab separated file containing annotations for coding sequences in the dog genome.

We can extract the header if of this file (lines starting with #) by running:
```sh
grep '^#' /home/DATA/module_1/dog_genes.gtf
```
> `Exercise 6`
>
> Remove the header of this file using the ”select non matching lines” option of grep (-v flag) and redirect the output to a file called `dog_genes_no_H.tsv` inside the stage_1/output directory.

The first line of your file without a header should look like:
`X ensembl gene 1575 5716 . + . gene_id "ENSCAFG00000010935"; gene_version "3"; gene_source "ensembl"; gene_biotype "protein_coding"`
It contains a lot of information that is not relevant for us at the moment. 
The fields that we are interested in are:
- The chromosome (X: 1st field)
- The type of the feature (gene: 3rd field)
- The starting position of the feature (1575: 4th field)
- The ending position of the feature (5716: 5th field)

> `Exercise 7`
> 
> Use `cut` to extract the required fields from dog_genes_no_H.tsv. Then redirect the output to a file called `dog_genes_table.tsv` inside your `stage_2/output/` directory. See cut --help to identify the option for fields

Now that we have extracted the relevant information, we would like to make a few adjustments to our table. Let’s start with adding the string `chr` at the beginning of each line.
We can do this easily by using the substitution command in `sed` which has the following general syntax:

```sh
sed 's/target/replacement/'
```
In our case we are going to modify the file in-place using the `-i` flag.

```sh
cd stage_2/output/
sed -i 's/^/chr_/' dog_genes_table.tsv
```
In our example, the caret symbol (^) is a regex which denotes the beginning of a line and we replaced this with `chr_`. You can check whether the substitution worked or not by examining the first 10 lines of the table with head.

The next thing we would like to do is switching the order of the columns in our table:
- 1. Chromosome
- 2. Starting Position
- 3. Ending Position
- 4. Feature Type

This requires a simple awk command:
```sh
awk 'BEGIN {OFS="\t"};{print $1,$3,$4,$2}' stage_2/output/dog_genes_table.tsv > stage_3/output/d_g_tab_cfp.tsv
```
`awk` is a powerful scripting language designed for processing text files. It is commonly used in Bash scripting to manipulate, analyze, and transform text data. `awk` is particularly well-suited for tasks like parsing log files, extracting information from text files, and performing text-based calculations.

The OFS option before the print command stands for ”Output Filed Separator” and we set it to Tab (`\t`) to ensure our table has the correct delimiter for a ”tab separated file” (TSV). `awk` stores each field in a different variable which is accessible via the `$` symbol.  We made use of this feature to put the columns in the right order.

We are almost done with pre-processing our data but there’s still something that’s not quite right with it. Have a look at the first column:

```sh
cd stage_3/output
cat d_g_tab_cfp.tsv | cut -f 1
```
Have you noticed that the chromosomes are not in the right order? Let’s fix it!
```sh
sort -V -o ../../results/d_g_sorted_table.bed d_g_tab_cfp.tsv
```
The `sort` command is a fundamental tool in Bash scripting used to organize and rearrange data based on specified criteria. It is commonly used to `sort` text files numerically or alphabetically, making it a versatile tool for data manipulation and analysis.

Here we used the -V option because we are dealing with a mixture of numerical and string data. Note also the -o to specify the output file which must precede the input. Our table is now ready to be analysed.

> `Bonus Exercise`
> 
> How many coding regions (CDS) on the X chromosome are listed in our bed file?
> 
> Use the commands you have learn to find out   

We will have more information about `BED` files in the next session.

### 1. Shell Scripting
In Session 1 we have seen how to navigate a Unix-like file system and how to manipulate text files. 
In this section we will revise what we have learnt about variables in Bash, we will then introduce the concept of positional arguments, and we will write our first shell script.

Before we start, let's run some preliminary commads to create the directory structure for this session:

```sh
cd; mkdir session2; cd session2; mkdir raw_data; mkdir scripts; mkdir results 
```
Note the use of `;`, which allows to run multiple commands in short succession.  

> Exercise 1
>
> Create a symbolic link between the `/home/DATA/module_1` folder and your newly created `raw_data` directory 

If you now move into your raw_data directory and run `ls module_1` you should see two files having the `.txt` extension, namely:
```sh
instructors_list.txt
participants_list.txt
```

Let's have a look at `instructors_list.txt` first, you can print the content on screen using `cat`.

Unfortunately, the fields are not well defined because each word is separated by a space. let's try to fix these formatting issues.

First of all we need to separate the last two fields (affiliation and status) from the instructors' names.
We can do this in `awk` and make use of variable `NF` (number of fields variable) which is set to the total number of fields in the input record:
```sh
awk '{print $(NF-1),"\t",$NF}' ~/session2/raw_data/module_1/instructors_list.txt
```
Now let's redirect the output to a file:
```sh
awk '{print $(NF-1),"\t",$NF}' ~/session2/raw_data/module_1/instructors_list.txt > ~/session2/raw_data/aff_status.txt
```
Tabs allow for consistent indentation across different file viewers and analysis tools. Each tab represents a single logical level of indentation, making the alignment structure more apparent and easier to interpret. Tabs are also more compatible with automation tools and scripts used for sequence analysis and manipulation. Many software programs and scripts expect tab-delimited alignment files, and using tabs avoids potential compatibility issues or the need for manual formatting adjustments.

Let's now deal with the names. Given that we don't know how many words each names consit of, we should start by printing all but the last two fields of the original input file:

```sh
awk 'NF-=2 {print $0}'  ~/session2/raw_data/module_1/instructors_list.txt
```

Then we can pipe this into `sed` and replace all white spaces with the character `_`:

```
awk 'NF-=2 {print $0}'  ~/session2/raw_data/module_1/instructors_list.txt | sed -e 's/ /_/g'
```
Note the use of the `g` at the end of the substitution command. The use of `g` means that sed will now change all matching space into a _.

Finally we redirect the output to a file:

```
awk 'NF-=2 {print $0}'  ~/session2/raw_data/module_1/instructors_list.txt | sed -e 's/ /_/g' > ~/session2/raw_data/names.txt
```

Now that we have created these two intermediate files we can stitch them together to reconstruct the initial information correctly formatted:
```sh
cd ~/session2/raw_data/
paste names.txt aff_status.txt > corrected_instructors_list.tsv
```
The `paste` command is a versatile tool that combines text from multiple files. It's particularly useful for merging data, and comparing files.

The series of commands presented above acts as a single unit to accomplish the required task. Therefore we can transform them into a script that we can re-use.
```sh
cd ..
touch ./scripts/formatting.sh
```
Note the `cd ..` command which moves the user up one directory (i.e. from `raw_data` to `session2` in this case).

Now we need to transform into an executable file. to do so run:
```sh
chmod 770 ./scripts/formatting.sh
```
In unix-like systms, `chmod` is the command and used to change the access permissions. In this case the owner of the file (i.e. you) should now be able to read (4), write (2), and execute (1) this file, hence the first 7 which is the sum of the three permission granted. The same is true for other users in the group (the second 7) while external user do not have any permission (0). 

Now let's use nano to edit our script and add this code-block to its content:

```sh
#!/usr/bin/bash
awk '{print $(NF-1),"\t",$NF}' ~/session2/raw_data/module_1/instructors_list.txt > aff_status.txt
awk 'NF-=2 {print $0}' ~/session2/raw_data/module_1/instructors_list.txt | sed -e 's/ /_/g' > names.txt
paste names.txt aff_status.txt >  ~/session2/results/corrected_instructors_list.tsv
```

This version of the script formatting.sh is not very useful because it can work only on the `instructor_list.txt` input file.
If we want to re-use it to format a different input file we would have to open it and edit the file name every time which is not convenient.
Let's modify it to allow for more flexibility in the usage by transforming input and output into variables:

```sh
#!/usr/bin/bash

INPUT_FILE=~/session2/raw_data/module_1/instructors_list.txt
OUTPUT_FILE=~/session2/results/corrected_instructors_list.txt

awk '{print $(NF-1),"\t",$NF}' $INPUT_FILE > aff_status.txt
awk 'NF-=2 {print $0}' $INPUT_FILE | sed -e 's/ /_/g' > names.txt
paste names.txt aff_status.txt > $OUTPUT_FILE
rm names.txt
rm aff_status.txt
```
Note that we have added two lines to delete the intermediate files (names.txt and aff_status.txt) that we don't need anymore.

This version looks slightly better but the input and output are still hard-coded inside the script. 
Ideally, we would like to supply the input and output at the runtime (meaning when we execute the script). To do so we can make use of positional arguments.
The indexing of the arguments starts at one, and the first argument can be accessed inside the script using $1. Similarly, the second argument can be accessed using $2, and so on.
Thus our final version of `formatting.sh` should be:
```sh
#!/usr/bin/bash

INPUT_FILE=$1
OUTPUT_FILE=$2

awk '{print $(NF-1),"\t",$NF}' $INPUT_FILE > aff_status.txt
awk 'NF-=2 {print $0}' $INPUT_FILE | sed -e 's/ /_/g' > names.txt
paste names.txt aff_status.txt > $OUTPUT_FILE
rm names.txt
rm aff_status.txt
```
now we can execute the script from the `session 2` directory and provide the correct input and output at the call:

```sh
./scripts/formatting.sh ./raw_data/module_1/instructors_list.txt ./results/corrected_instructors_list.txt
```

> Exercise 1
>
> Use this script to correctly format the `participants_list.txt`. Make sure you input the right path to this file and use a sensible name for the output file (second argument)
 
