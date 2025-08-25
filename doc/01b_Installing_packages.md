# Installing packages on UNIX

## Installing packages on WSL (Debian), using the apt package manager

```powershell
# Open the Windows PowerShell and type:
$ wsl -u root  # Entering Debain with root rights enabled

# When in Debian, type:
# Principle: apt search myReallyCoolSoftware  # Example of how to search for tools

$ apt install agrep

$ apt install ncbi-entrez-direct

$ apt install sra-toolkit
```

## Installing packages on MacOS, using the brew package manager
```bash
# Installing the package manager 'brew':
$ /bin/bash -c "$(curl -fsSL \
	https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Installing new command-line tools:
$ brew install tre   # includes agrep

$ brew search ncbi-entrez-direct
$ brew install name_of_package_here

$ brew install sratoolkit
```
