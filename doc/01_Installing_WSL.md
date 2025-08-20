# Windows Subsystem for Linux (WSL)

## Installing WSL (Debian) via the Windows PowerShell

```powershell
# Open the PowerShell as an administrator

# 1. Enable WSL feature
$ dism.exe /online /enable-feature `
    /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart

# 2. Enable Virtual Machine Platform (required for WSL2)
$ dism.exe /online /enable-feature `
    /featurename:VirtualMachinePlatform /all /norestart

# 3. Set WSL2 as default version
$ wsl --set-default-version 2

# 4. Install Debian
$ wsl --install -d Debian
```

## Initial WSL Setup
```powershell
# Launch WSL
$ wsl

# First-time setup prompts:

# 1. Enter a new UNIX username:
#    Example: myname

# 2. Enter a password for this user
#    (input is hidden, type carefully)

# 3. Confirm the password

# After setup, you are logged in as your new user:
$ whoami   # Verify your username
$ pwd      # Check your current directory (home directory)
```

## If you have a non-US keyboard: Changing locale in WSL
```powershell
# 1. Check current locale
$ locale

# 2. Reconfigure locales
$ sudo dpkg-reconfigure locales
#  Select desired locales (e.g., en_US.UTF-8)

# 3. Apply new locale
$ export LANG=en_US.UTF-8
$ export LANGUAGE=en_US:en
$ export LC_ALL=en_US.UTF-8
#  Add to ~/.bashrc or ~/.profile for permanence

# 4. Verify
$ locale
```

Installing additional tools on MacOS, using the brew package manager
```bash
# Installing the package manager 'brew':
$ /bin/bash -c "$(curl -fsSL \
	https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Installing new command-line tools:
$ brew install tre   # includes agrep

$ brew search ncbi-entrez-direct
$ brew install ...

$ brew install sratoolkit
```
