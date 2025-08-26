# Introduction to Unix Commands

## Unix System Commands (available on all Unix systems)

### Time, date & user info

```bash
# Screen
$ clear
$ history

# Time and Date
$ date
$ sleep
$ uptime

# User info
$ whoami
$ id
$ groups
$ who
$ last
```

### Navigating the filesystem and simple file operations

```bash
# Navigating filesystems
$ cd
$ cd <dir>             # change directory to <dir>
$ pwd                  # print working (i.e. current) directory
$ mkdir
$ rmdir

# File operations
$ ls <dir>             # list contents of <dir>
$ cp <src> <dest>      # copy <src> to <dest>
$ rm <file>            # remove <file>
$ mv
$ chmod
$ chown
```

### More on navigating the filesystem

```bash
# Listing directories and contents
$ ls -alh
total 120
drwx--x---    7 gruenstaeudl homedir  4096 Feb  5 19:38 .
drwxr-x--x 2273 root         root    57344 Feb  5 10:23 ..
drwxr-x---    2 gruenstaeudl users       6 May 20  2015 .addressbook
-rw-------    1 gruenstaeudl users    6200 Feb  5 18:51 .bash_history
-rw-r--r--    1 gruenstaeudl users     127 Feb  5 19:38 .bash_profile
-rw-r--r--    1 gruenstaeudl users      17 Nov 14  2017 .bashrc
drwxr-x---    3 gruenstaeudl users      21 May 20  2015 .calendars
drwxr-xr-x    2 gruenstaeudl users       6 Feb  5 18:50 .nano
...

# Symbolic directories
.     # working (i.e. current) directory
..    # directory above working directory
-     # last directory (i.e. before last cd)
~     # home directory
/     # root directory

# EXAMPLES:
$ pwd
/home/g/gruenstaeudl
$ cd ..
/home $ pwd
/home/g
```

### Viewing and modifying text files

```bash
# Text file viewing
$ more    # viewing text one screenful at a time
$ less    # same as command 'more' but with more options

$ head    # display only the first n lines of a file
$ tail    # display only the last n lines of a file

$ cat <file1> <file2>  # concatenate <file2> to <file1> and view

# Text editor
$ pico <file1>    # Editing a file
$ nano <file1>    # Editing a file but with more options

# Within nano or pico:
# Ctrl+O    Saving
# Ctrl+X    Exiting
```

### Exercise: Downloading a Genome Record from GenBank using UNIX
```bash
# Download the full genome record of SARS-CoV-2 from GenBank 
# using the following two commands:

$ i=NC_045512.2
$ curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi? \
	db=nucleotide&id=${i}&rettype=gb&retmode=txt" \
	> $i.gbk

# Note: The curl-command must NOT have a line break.
```
