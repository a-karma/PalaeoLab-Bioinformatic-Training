![bio_logo](../IM/header.png)
### 2. Shell Scripting
In this section we will revise what we have learnt about variables in Bash, we will then introduce the concept of positional arguments, and we will write our first shell script.

Before we start, let's run some preliminary commads to create the directory structure for this session:

```sh
cd module1; mkdir shell_scripting; cd shell_scripting; mkdir raw_data; mkdir scripts; mkdir results 
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


If you finished early and you want to keep playing with bedtools have a look at this tutorial here: https://sandbox.bio/tutorials?id=bedtools-intro
