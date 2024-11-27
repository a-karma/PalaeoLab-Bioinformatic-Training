![bio_logo](../IM/header.png)

## Module 1 - Basic concepts of command line programming - Session 0

The vast majority of genomic analyses require large computational resources way beyond the computing power of a standard laptop. Therefore, we have generated a virtual cloud computer with all the software and the data necessary for the course. You'll need to connect to this server via a `Secure SHell` protocol (ssh), and you’ll run all the commands directly in a terminal. This section provides you with some essential definitions and outlines the procedure for establishing a connection to our cloud server. 

### 1. Preliminary Definitions
- The `OS` (short for Operating System) is the program that manages all other applications in a computer. Windows, Linux and macOS are all examples of Operating Systems. Users can interact directly with the operating system through an interface, such as a command-line interface (CLI) or a graphical user interface (GUI).
- A `shell` is a text-based command-line interpreter or a program that allows the user to execute commands and interact with the OS. Different operating systems have different shell programs, here we will focus on Bash (Bourne Again SHell) which is the most common shell on Unix-like systems.
- A `terminal` is the application that we use to interact with the `shell` i.e. the window were you type commands.
- The `prompt` is the text next to where you type your commands in a terminal
- A `file system` is the structure and logic rules used by the OS to control how data is stored and retrieved.
- A `script` is a series of instructions that automate computer tasks  Here we will focus on shell-scripting meaning executable files interpreted by the shell. 

### 2. Pre-course work
Given that the cloud server is UNIX based and that almost all bioinformatics procedures are conducted via CLI, it is important that you familiarise yourself with it before starting the course. Although we will cover some of the basic concepts of file system structure and shell scripting in this module, if you have never worked from a terminal before we highly recommend you to go through this introductory courses:

- https://edu.sib.swiss/pluginfile.php/2878/mod_resource/content/4/couselab-html/content.html
- https://sandbox.bio/tutorials?id=terminal-basics

> [!NOTE]
> Both courses can be done using only a web browser.

### 3. Setting up your local machine
As we only have a week and a lot of material to cover, unfortunately we will not be able to cover
the basics of R. This is a program we use, mainly for the visualisation of the results. For those of
you who are not familiar with R, it would be very useful to have a go before we all get together.
First everyone needs to download and install the most up to date versions of R and Rstudio on
your laptop:
Rstudio website: https://posit.co/download/rstudio-desktop/
If you need help for the installation you can follow this video:
Rstudio installation: https://www.youtube.com/watch?v=TFGYlKvQEQ4
If you are completely new to R, we would recommend you go through this video before the
workshop: https://www.youtube.com/watch?v=BvKETZ6kr9Q
Extra material: Reference for data visualisation with R
https://r4ds.had.co.nz/data-visualisation.html
Last week we sent you the instructions for installing and setting up R and RStudio. There are
now a few packages that you need to install for the practical sessions. Could you please try to
install the following packages - via the command install.packages(“package_name”):
- tidyverse
- ggplot2
- reshape2
- Rcolourbrewer
- dplyr


### 4. Connecting to the cloud server
##### Windows Users
In the guide we've sent you, you've seen how to install putty and how to connect to our cloud server using this program. 
 - open the putty app
 - under the Saved Session box, click on `PalaeoLab` then on the `Open` button.
 - In the terminal window type your user_ID and press enter
 - You should be now asked for the key-passphrase, please type it and press enter

##### Linux & macOS user
Please open a terminal and then run:

```sh
ssh -i ./.ssh/palaeolab_key -l user_ID 138.246.238.65
```
where user_ID needs to be replaced by the credentials we sent you.
