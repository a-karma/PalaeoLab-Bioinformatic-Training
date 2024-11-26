![bio_logo](../IM/header.png)

## Module 1 - Basic concepts of command line programming - Session 0

The vast majority of genomic analyses require large computational resources way beyond the computing power of a standard laptop. Therefore, we have generated a virtual cloud computer with all the software and the data necessary for course. You'll need to connect to this server via a Secure SHell protocol (ssh), and you’ll run all the commands directly in a terminal. In this section, we will outline the procedure for connecting to our cloud server and provide you with some essential definitions. 

### 1. Preliminary Definitions
- The `OS` (short for Operating System) is the program that manages all other applications in a computer. Windows, Linux and macOS are all examples of Operating Systems. Users can interact directly with the operating system through an interface, such as a command-line interface (CLI) or a graphical user interface (GUI).
- A `shell` is a text-based command-line interpreter or a program that allows the user to execute commands and interact with the OS. Different operating systems have different shell programs, here we will focus on Bash (Bourne Again SHell) which is the most common shell on Unix-like systems.
- A `terminal` is the application that we use to interact with the `shell` i.e. the window were you type commands.
- The `prompt` is the text next to where you type your commands in a terminal
- A `file system` is the structure and logic rules used by the OS to control how data is stored and retrieved.
- A `script` is a series of instructions that automate computer tasks  Here we will focus on shell-scripting meaning executable files interpreted by the shell. 

### 2. Preliminary set up
In order to establish a connection and interact with the server you will need individual credentials: these will be sent to you next week along with detailed instructions to guide you during the process.
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
