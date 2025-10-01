### Access to Beocat

#### Opening SSH connection to Beocat
```
$ ssh username@headnode.beocat.ksu.edu
Password: 
...
[username@selene ~]$ 
[username@selene ~]$ pwd
/homes/username 
```

#### Available Software
```
module overview   # lists all available modules (summary form)

module avail   # lists all available modules (detailed form)

module list   # lists all currently loaded modules
```

#### Finding Specific Software
```
module keyword a_string   # search modules containing \textit{string}

module spider a_string   # list modules matching \textit{string}

module whatis module_name   # display information about a module
```

#### Installing Miniconda
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh

bash ~/miniconda.sh -b -p $HOME/miniconda3

$HOME/miniconda3/bin/conda init bash
```

#### Closing SSH connection to Beocat
```
[username@selene ~]$ exit
logout
Connection to headnode.beocat.ksu.edu closed.
$
```

#### Connecting to Beocat via sftp
```
$ sftp username@headnode.beocat.ksu.edu
Password:
Connected to headnode.beocat.ksu.edu.
sftp> pwd
Remote working directory: /homes/mgruenstaeudl
sftp>
```

#### Local vs. remote sftp commands
```
# Show working directory on remote system (Beocat)
sftp> pwd
# Show working directory on local computer (laptop)
sftp> lpwd

# Change directories on remote system
sftp> cd ..
sftp> pwd
sftp> cd ~

# Change directories on local computer
sftp> lcd ..
sftp> lpwd
sftp> lcd ~
```

#### Local vs. remote sftp commands -- Cont'd
```
# List files in current remote directory (Beocat)
sftp> ls
# List files in current local directory (laptop)
sftp> lls

# Create and remove a temporary directory on remote system
sftp> mkdir temp99
sftp> ls
sftp> rmdir temp99

# Attempt same on local system
sftp> lmkdir temp99
sftp> lls
sftp> lrmdir temp99   # Note: invalid command, no 'lrmdir' in sftp
```

#### Transferring files between local and remote
```
# Upload file from local → remote
sftp> lls
spike_glycoprotein.fasta
sftp> put spike_glycoprotein.fasta
Uploading spike_glycoprotein.fasta to /homes/mgruenstaeudl/spike_glycoprotein.fasta
spike_glycoprotein.fasta                      100%  710   42.9KB/s   00:00

# Download file from remote → local
sftp> lcd temp99
sftp> get spike_glycoprotein.fasta
Fetching /homes/mgruenstaeudl/spike_glycoprotein.fasta to spike_glycoprotein.fasta
/homes/mgruenstaeudl/spike_glycoprotein.fasta 100%  710    4.7KB/s   00:00

# Close sftp session
sftp> bye
```
