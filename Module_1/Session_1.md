![bio_logo](../IM/header.png)

## Module 1 - Basic concepts of command line programming - Session 1

### 1. The filesystem 
This section outlines the basic commands for navigating file systems and creating a structured directory hierarchy for research projects. The naming of files and directories is more crucial than you might realize; a clear and consistent structure will significantly streamline your workflow, especially when dealing with complex projects and pipelines. While you can adapt and customize this template structure to your specific needs, adhering to these principles will greatly enhance project organization and accessibility.

We will start with creating a parent directory for our project. 
```sh
mkdir project_bash
cd project_bash
```
The first command (`mkdir`) creates the required folder while the `cd` command is used to move from your home directory to the newly created directory `project_bash`. To make sure you are now in the right directory, you can use the `pwd` command which prints on screen the name of the current directory *i.e.* your current position in the file system. 

> [!WARNING]
> Never use white spaces when naming files or directories

Now we need to create one directory where we are going to store our scripts, one for the raw data, and one for our results:

```sh
mkdir scripts; mkdir raw_data; mkdir RaSuLts
```
Note the use of the semicolon (;) to separate different commands on the same line.

To assure yourself that all directories have been created you can use the command `ls` which lists the content of the current directory.

Do you notice anything strange? 

The word `results` is misspelt and the directory name is written in a mix of upper and lower case letters. Let's fix this!

One option is to delete the folder and create a new one with the correct name. You can achieve this by running:
```sh
rm -r RaSuLts
mkdir results
```
A second possibility is to move the directory (and its content) to a new directory with the correct name:
```sh
mv RaSuLts results
```
Now that we have corrected our mistake, let's do some more housekeeping to clarify the notion of `path`:
```sh
cd
mkdir module1
cd module1
mv ../project_bash/ ./
```
Now run `pwd` and `ls`: do you see what happened? Let's explain this trick!

The first `cd` command without any dir name is a shortcut to move from any position in the filesystem to your `home` directory.

The `mkdir` command creates a new directory called `module1` and the `cd module1` command allows us to move inside this newly created directory.

Lastly the `mv` command transfers the `project_bash` folder (and all its content) from your `home` directory to the `module1` directory. 

Note how the two paths (source and destination) have been specified: 
- source = `../project_bash/`
- destination = `./`
  
These paths are usually called `relative` in the sense that the location in the filesystem is specified relative to the current directory *i.e.* `module1`. Let's explain this concept in more details:

The `project_bash` directory was located inside your home dir, hence, from the `module1` dir we needed to move upstream one directory in order to find the folder called `project_bash`. This is signalled by the `../` symbols whch allow us to access the parent directory (*i.e.* your home dir). The `./` characters represent instead a general and concise way of specifying the current directory.
 
> `Exercise 1`
>
> Which command should you use to:
> 1. go to your home directory?
> 2. move to the `raw_data` directory from your home dir?
> 3. move to the `scripts` dir from `raw_data`?
> 
> Use relative paths to specify the directories

We can now begin with populating our directories starting with a `log` file. Maintaining a detailed record of executed commands is an essential practice for reproducibility and debugging purposes. By keeping a log of your actions, you can easily retrace your steps, identify potential errors, and effectively reproduce your work, even after a significant time gap.

To facilitate efficient command logging, it's recommended to adopt a consistent naming convention for your log files. In this case, we'll use the file name "what_i_did.txt" to store our command history. To create this empty log file, simply utilize the touch command:

```sh
cd; cd module1/project_bash/
touch what_i_did.txt
```
You can then edit this file manually using any text editor, such as `nano`.

```sh
nano what_i_did.txt
```
After editing the file you can then press `Ctrl X` to close the editor and press `enter` to save the changes.

> `Exercise 2`
>
> log all the command you have run so far in your `what_i_did.txt` file

You can visualize the content of this file using the command `cat path/to/file/file_name`

>[!TIP]
> The tab key triggers autocomplete in bash which is extremely useful to avoid typing complex and long paths or file-names!

### 2. Loops and variables
Let's keep building up our file-system structure by creating a separate directory for each stage of our pipeline/analysis. Because of my lack of immagination, each directory will be termed as `stage` and we are going to number them sequentially but I would encourage you to use a better naming convention when it comes to your research. Given that typing the same command over and over is tedious, we are going to use one of the basics building blocks of any programming language: a loop.  

```sh
for i in $(seq 3)
do
mkdir stage_$i
done
```
Let's unpack this:

The `for` loop is a fundamental programming construct that allows you to execute a block of code repeatedly until a specific condition is met. In this particular instance, the loop iterates over the values generated by the seq command. 

The `seq` command produces a sequence of integers, and in this case, it generates the numbers from 1 to 3.

The `do` command marks the beginning of the block of code that will be executed repeatedly for each value of the iterator variable.

The variable `i` serves as the iterator, keeping track of the current position in the sequence. With each iteration, the value of `i` changes to the next number in the sequence, until it reaches the end of the sequence.

The `mkdir`command: Creates a directory each time (3 directories in total). 

The `done` command marks the end of the for loop, indicating that the loop should continue executing the enclosed block of code until the entire sequence of values has been processed.

This will create 3 directories called: `stage_1`, `stage_2`, `stage_3`.

Now type `ls` to ensure the directories have been indeed created.

Today we will create a very simple directory structure for our research project just to illustrate the principle. Each of our `stage` directories should contain three items:
- a what_i_did.txt file.
- an output sub-directory
- a sub-directory called input

In principles, you could navigate to the each stage directory using the `cd` command and create these objects manually but that involves a lot of typing. You should instead use a loop to avoid this tedious task. An efficient way to do this from a shell terminal requires to create a list of parent directories and then create the child directories only where we need them.
```sh
ls -d stage* > dir_list.txt
```
Here the `-d` option modifies the behaviour of the command `ls` and forces it to list only directories. 

Note the word `stage*` after the `-d` flag: this restricts our list of directories to those that are termed `stage` followed by any other character (`*`). 

Lastly, note the use of the `>`. This symbol in bash has a special meaning: redirect the output of the command that precedes it to a file (`dir_list.txt` in this case).

Now that we have the list, we can easily implement another type of loop using the reserved word `while`:
```sh
while read -r line
do
touch $line/what_i_did.txt
done < dir_list.txt
```
The `$` sign preceding the variable named `line` is crucial in this context because it indicates that we are referencing the actual value stored in the variable, rather than the variable name itself. During each iteration of the loop, the `read` command reads a line from the `dir_list.txt` file and assigns its content to the `line` variable. Finally the command `touch` uses the content of the variable `$line` to create a `what_i_did.txt` file inside the right directory.

> `Exercise 3`
>
> use a similar while loop to create a sub-directory called `output` inside each stage directory

In data processing pipelines, where the output of one step becomes the input for the next, it's often necessary to transfer data between directories. However, blindly copying files can lead to inefficient storage usage and unnecessary data duplication. A more effective approach is to utilize symbolic links, which provide a pointer to the original data rather than creating a new copy.

Consider the scenario where `stage_1/output` and `stage_2/input` contain identical data. Copying the entire contents from `stage_1/output` to `stage_2/input` would result in redundant storage of the same data, consuming unnecessary disk space.

Symbolic links offer a solution to this issue. Instead of replicating the data, a symbolic link is created in `stage_2/input`, pointing to the original data in `stage_1/output`. This link acts as a shortcut, directing the system to the actual data location whenever the input directory is accessed.

Here is an example of symbolic links:

```sh
ln -s ~/module1/project_bash/raw_data/ ~/module1/project_bash/stage_1/input
```
The `ln` command stands for ”link” and it has this general syntax:
```sh
ln full/path/to/source full/path/to/destination/link_name
```
Here we have used the `-s` flag to specify a symbolic link between the `raw_data` directory and a new `input` folder (link name) inside the `stage_1 directory`. This means that the content of `raw_data` is now accessible from `stage_1/input`. Let’s double-check:
```sh
touch raw_data/input_zero.txt
ls stage_1
ls stage_1/input
ls raw_data
```
The first command is just to populate the `raw_data` folder with a file (`input_zero.txt`). By running the second command you should see that an `input` folder has been created via the `ln` command. Now the `stage_1` directory contains all three elements required. The output of the third and fourth commands should be just `input_zero.txt`.

> `Exercise 4`
>
> Use the `ln` command to complete our file system structure
> by linking each `stage_{i}/output` to an input folder inside the `stage_{i+1}` directory.
> Note that in our example stage_4 is actually termed results.

If you have made a mistake with links, do not panic. You can alwayse remove them with `rm` or with the `unlink` command.

>[!NOTE]
> When running the `ln` commands to create our links we have used the abslolute paths to both the source and the destination.
>
>The `~/` characters at the beginning of the paths represent a short cut which stands for `/home/userID/` directory. 

Now that we have a good structure we can start populating our directories. Let's create some files in the `project_bash` directory.

```sh
for i in $(seq 3)
do
touch output_file_$i
done
```
Now we should move each of these file to the corresponding `stage_{i}/output` directory.

> `Exercise 5`
>
> Use a for loop to move each of the output_file_{i} ∀i ∈ {1, 2, 3} to its own directory.
> 
> To do so, you should use the `mv` command which has the following syntax: `mv target_file_name path/to/destination`.


If you now run the `tree` command from the project_bash directory, you should get:

![File-system-structure](../IM/bash_tree.png)


### 3. Conditionals and variables
Another fundamental structure of any programming language are conditional expressions (a.k.a. if statements). They allow us to make decisions and run a piece of code only when it is needed. The general architecture of a conditional construct in bash is:

![Conditionals in bash](../IM/conditionals.png)

There are essentially three syntax rules for conditionals in bash: 
- Always leave a space before and after the conditional expression.
- Always leave a space before and after the comparison operator
- Always terminate line before adding a key word.

Plus there is a good practice rule when dealing with strings: 
- Remember to quote string variables.

This will force the shell to interpret correctly escaped character such as tabs (`\t`), new-lines (`\n`), white-spaces (`\s`), etc.

Let's illustrate these rules and how conditional epressions work with an example:

> Example 1
```sh
a=2; b=3
if [ $a -l $b ]; then
echo "The first number is less than the second one"
fi
```
In the code above we first define two variables (`a` and `b`) and assign to each of them a value (`2` and `3` respectively).
We then build our conditional starting with the reserved word `if` followed by the condition we want to evaluate. In this case we want to check whether the value assigned to the variable `a` (*i.e.* `$a`) is less than (`-l`) the value assigned to the variable `b` (*i.e.* `$b`). Note the space before and after the comparison operator its `$a -l $b` and not `$a-l$b`. We put the condition within square brakets making sure we left a white space between the brackets and our condition. We terminate the line (`;`) before using the other reserved word `then` which introduces our command (`echo`). Finally we close our if statetement with the reserved word `fi`.

The command `echo` followed by a string will be executed only if our condition is `true` and in that case it will result in printing on screen the message `The first number is less than the second one`. If instead the condition is `false` nothing happens.

We can also provide alternatives in our conditionals meaning that we can ask the shell to execute a different command when the condition is false:

> Example 2
```sh
a=2; b=3
if [ $a -g $b ]; then
echo "The first number is bigger"
elif [ $b -g $a ]; then 
echo "The second number is bigger"
else
echo "The two numbers are equal"
fi
```

In this second example we have created a slightly more complex conditional that which covers all possible scenarios that can happen when comparing two numbers. The program will start by evaluating the first condition: " `a` is bigger than `b`". If this statement is true, it will print the first message and terminate its execution. If instead it turns out that the statement is false, the shell will move on and evaluate the second condition: "`b` is bigger than `a`". If this second statement is true it will print the second message and then again terminate. The last case scenario will happen only if both conditions are false which will cause to print on screen the third message and exit.

Complex expression are created by combining multiple condition with logical operators such AND (`&&`) and OR (`||`).

> Example 3
```sh
a=2; b=2
if [ $a -le $b ] && [ $a -ge $b ]; then
echo "The two numbers are equal"
fi
``` 
In this last example we ask the shell to evaluate two condition before making a decision: the command `echo` will be exectued only if both conditions are true. This behaviour is obviously a consequence of the logical operator AND which connects our two conditions.

Now that we have seen ho to built conditionals it's time for practice!

Let's first use what we learn to crate a list of all sub-directories:

> Exercise 6
>
> Go to the `project_bash` dir and list all its content.
> 
> Redirect the stdout to a file called `all_content_list.txt`

We would like to add a file termed `final_results.txt` to the directory `results` and a file called `bash_script.sh` to the directory `scripts`. In order to practice with conditionals, we are now going to parse line by line the file you've generated in `Exercise 6` using a while loop similar to the one you have seen in the previous sections. Each time we will need to evalute two conditions, see pseudo-code below: 

```sh
if line == "results"
create the final_results.txt file inside the directory results
elif line == "scripts"
create the bash_script.sh file inside the directory scripts
```

If the content of the variable `line` is different from both `results` and `scripts`, nothing happens.

> Exercise 7
>
> Modify the following while loop and add the proper conditional using the correct bash syntax:
> make sure you leave a white space before and after the double equal sign (` == `) as well as after `[` and before `]`. 

```sh
while read -r line
do
<insert conditionals here>
done<all_content_list.txt
```
 
If you correctly solved all the exercise in this Session, you should get the following output when running the `tree` command from your `project_bash` directory:

![File-system-structure](../IM/bash_tree_2.png)

This concludes our session, hope you had fun! 


