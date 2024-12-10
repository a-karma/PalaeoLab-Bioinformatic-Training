![bio_logo](../IM/header.png)
## Module 1 - Basic concepts of command line programming - Session 3
In this last session of Module 1 we will revise what we have learnt about variables and text manipulation, we will then introduce the concept of positional arguments, and we will write our first shell script.

### 1. Shell Scripting - Theory
A script is just a file that contains a series of commands. A `bash script` is a series of command interpreted by the shell. 
This means that it needs to comply with shell syntax rules or in other words it is written in bash language. 

- The first line is called shebang: `#!/bin/bash`.   
  It means that the script will be executed using the Bash shell (located in the `/bin` folder) as interpreter.
- Lines starting with `#` are interpreted as comments and will be ignored.
- Multiple statements on a single line must be separated by a `;`
- To split a single statement (command) on multiple lines use `\`
- No indentation needed (but might be useful to make your code more readable)
- Call a script by invoking `bash` followed by the script name (with its path)

#### Flexibility
Although each script is designed to find a solution to a specific problem, it needs to posses some flexibility so that it can be used to process different input with a similar data structure.  
For example: imagine that you want to write a script to calculate the sum of 2 plus 3. Even if your script returns the correct value (*i.e.* 5) it is rather useless because of its lack of generalization. On the contrary, a script that calculates the sum of any two integer is way more useful because it solves a whole class of problems.

In order to mantain a certain level of generalization, you will often need to define and use variables in your script:

- A variable is a character string to which we assign a value.readability
- Variables are symbolic references to chunks of memory *i.e* they are pointers to the actual data. 
- Variables’ names are traditionally declared in capital letters.
- Values assigned to a variables are made accessible to commands using the $ sign followed by the variable name.

>[!NOTE]
> Variable names are complitely arbitrary in the sense that they will work no matter what string of characters you use.  
> Choosing a meaningful naming convention though can greatly improve the readability of your code and make dubugging easier.

Another important feature that promotes coding flexibility is the use of `arguments` (a.k.a. *args*). Args are parameters that can be passed to the script at the call as a white-space separated list. Args can be accessed inside the script using the dollar sign followed by their index meaning their position in the list. Therefore, `$1` refers to the first argument, `$2` represents the second and so on. Let's consider the following example to illustrate the concept of positional arguments.

```sh
cat arg_tester.sh
```
The comand above will output on screen the content of the script called `arg_tester.sh` which looks like this:
```sh
#!/bin/bash

echo "Username: $1"
echo "Age: $2"
echo "Status: $3"
```
Now if we call the script from its directory using the following command:
```sh
bash arg_tester.sh Jason 27 "phd student"
```
We obtain the following output on screen:
```sh
Username: Jason
Age: 27
Status: phd student
```
This is because each `echo` command will process a different argument in the order that they have been passed to our script at the call. Note that because our third argument consisted of two words separated by a space we had to use quotes to ensure that both words are interpreted as a single argument.

#### Modularity
A script will make use of the building blocks of programming such as conditionals and loops to make decisions and automate a series of tasks.
A group of commands that act a unit constitutes a module that can be recycled in the future.
These modules are often re-written as `functions` so that they can be easily included in other scripts.  

- Function can be thought as a small script within a script.
- They are particularly useful when a given task need to be performed several times.
- The function definition needs to appear in the script before any function call.

Functions in bash have either of the following syntaxes:

```sh
# Option 1
function function_name {
# function body
<insert code here>
}
```
```sh
# Option 2
function_name(){
# function body
<insert code here>
}
```
Input data can be passed to the function as arguments directly after calling the function in your script. As for scripts, arguments are accessible inside the function as `$1`, `$2`, etc. according to the order in which they are passed to the function at its call. Let's have a look at a simple script in which we have declared and used a function:

```sh
cat funct_tester.sh
``` 
The comand above will output on screen the content of the script called `funct_tester.sh` which looks like this:

```sh
#!/bin/bash

FILE_INPUT1="long_text.txt"
FILE_INPUT2="short_text.txt"

line_counter(){
number_of_lines=$(cat $1 | wc -l)
echo "The file you chose is called $1"
echo "This file has $number_of_lines lines"
}

line_counter $FILE_INPUT1
line_counter $FILE_INPUT2
```
Before examining the function itself, let's have a look at the general structure of the script:

- We first define two variables (`FILE_INPUT1` and `FILE_INPUT2`).
- Then we declare our function called `line_counter`.
- Finally we called the function twice, each time passing to the function a different argument (`$FILE_INPUT1` at the first call and `$FILE_INPUT2` at the second call)

Now let's look at the function body *i.e.* the code enclosed in curly brackets which consists of three lines:
- Line 1 defines the variable `number_of_lines`
- Line 2 consists of a command (`echo`) in which we use the first argument that needs to be provided to the function at its call (`$1`)
- Line 3 consists of a command (`echo`) in which we use the value of the variable defined on the first line ($number_of_lines lines)




If we want to test our script on two previously created filed termed `long_text.txt` and `short_text.txt` we need to run:
```sh
bash funct_tester.sh
```
and its output on screen will look like this:
```sh
The file you chose is called long_text.txt
This file has 54 lines
The file you chose is called short_text.txt
This file has 5 lines
```
### 2. Shell Scripting - Example
Before we start, let's run some preliminary commads to create the directory structure for this session:

```sh
cd ~/module1; mkdir shell_scripting; cd shell_scripting; mkdir raw_data; mkdir scripts; mkdir lists; mkdir results

ln -s  /home/DATA/module_1 ~/module1/shell_scripting/raw_data/
```

If you now run `ls raw_data` you should see two files having the `.txt` extension, namely:
```sh
instructors_list.txt
participants_list.txt
```

Let's have a look at `instructors_list.txt` first, you can print the content on screen using `cat`.

As you can see, the fields are not well defined because each word is separated by a space. Let's try to fix these formatting issues.

First of all we need to separate the last two fields (affiliation and status) from the instructors' names.
We can do this in `awk` and make use of variable `NF` (number of fields variable) which is set to the total number of fields in the input record:
```sh
awk '{print $(NF-1),"\t",$NF}' ./raw_data/instructors_list.txt
```
Now let's redirect the output to a file:
```sh
awk '{print $(NF-1),"\t",$NF}' ./raw_data/instructors_list.txt > ./lists/aff_status.txt
```
Tabs allow for consistent indentation across different file viewers and analysis tools. Each tab represents a single logical level of indentation, making the alignment structure more apparent and easier to interpret. Tabs are also more compatible with automation tools and scripts used for sequence analysis and manipulation. Many software programs and scripts expect tab-delimited alignment files, and using tabs avoids potential compatibility issues or the need for manual formatting adjustments.

Let's now deal with the names. Given that we don't know how many words each names consit of, we should start by printing all but the last two fields of the original input file:

```sh
awk 'NF-=2 {print $0}'  ./raw_data/instructors_list.txt
```

Then we can pipe this into `sed` and replace all white spaces with the character `_`:

```
awk 'NF-=2 {print $0}'  ./raw_data/instructors_list.txt | sed -e 's/ /_/g'
```
Note the use of the `g` at the end of the substitution command. The use of `g` means that sed will now change all matching space into a _.

Finally we redirect the output to a file:

```
awk 'NF-=2 {print $0}'  ./raw_data/instructors_list.txt | sed -e 's/ /_/g' > ./lists/names.txt
```

Now that we have created these two intermediate files we can stitch them together to reconstruct the initial information correctly formatted:
```sh
cd ~/module1/shell_scripting/lists/
paste names.txt aff_status.txt > corrected_instructors_list.tsv
```
The `paste` command is a versatile tool that combines text from multiple files. It's particularly useful for merging data arranged in columns with an equal number of rows.

The series of commands presented above acts as a single unit to accomplish the required task. Therefore we can transform them into a script that we can re-use.
```sh
cd ..
touch ./scripts/formatting.sh
```
Note the `cd ..` command which moves the user up one directory (i.e. from `lists` to `shell_scripting` in this case).

Now we need to transform into an executable file. to do so run:
```sh
chmod 770 ./scripts/formatting.sh
```
In unix-like systms, `chmod` is the command and used to change the access permissions. In this case the owner of the file (i.e. you) should now be able to read (4), write (2), and execute (1) this file, hence the first 7 which is the sum of the three permission granted. The same is true for other users in the group (the second 7) while external user do not have any permission (0). 

Now let's use nano to edit our script and add this code-block to its content:

```sh
#!/usr/bin/bash
awk '{print $(NF-1),"\t",$NF}' ./raw_data/instructors_list.txt > ./lists/aff_status.txt
awk 'NF-=2 {print $0}'  ./raw_data/instructors_list.txt | sed -e 's/ /_/g' > ./lists/names.txt
paste ./lists/names.txt ./lists/aff_status.txt >  ./lists/corrected_instructors_list.tsv
```

This version of the script formatting.sh is not very useful because it can work only on the `instructor_list.txt` input file located in the `raw_data` directory. Moreover it will work only from the directory `shell_scripting` and providing that the directory structure is the same because we have specified relative paths. If we want to re-use it to format a different input file we would have to open it and edit the file name every time which is not convenient. Let's modify it to allow for more flexibility in the usage by transforming input and output into variables:

```sh
#!/usr/bin/bash

INPUT_FILE=~/module1/shell_scripting/raw_data/instructors_list.txt
OUTPUT_FILE=~/module1/shell_scripting/lists/corrected_instructors_list.txt

awk '{print $(NF-1),"\t",$NF}' $INPUT_FILE > aff_status.txt
awk 'NF-=2 {print $0}' $INPUT_FILE | sed -e 's/ /_/g' > names.txt
paste names.txt aff_status.txt > $OUTPUT_FILE
rm names.txt
rm aff_status.txt
```
Note that we have added two lines to delete the intermediate files (names.txt and aff_status.txt) that we don't need anymore.

This version looks slightly better but the input and output are still hard-coded inside the script. 
Ideally, we would like to supply the input and output at the runtime (meaning when we execute the script). To do so we can make use of positional arguments.

As explained in our previous section, indexing of the arguments starts at one, and the first argument can be accessed inside the script using $1. Similarly, the second argument can be accessed using $2, and so on. Thus our final version of `formatting.sh` should be:
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
Now we can execute the script from the `shell_scripting` directory and provide the correct input and output at the call:

```sh
./scripts/formatting.sh ./raw_data/instructors_list.txt ./results/corrected_instructors_list.txt
```

> Exercise 1
>
> Use this script to correctly format the `participants_list.txt`. Make sure you input the right path to this file and use a sensible name for the output file (second argument)

### Shell Scripting - Practice
