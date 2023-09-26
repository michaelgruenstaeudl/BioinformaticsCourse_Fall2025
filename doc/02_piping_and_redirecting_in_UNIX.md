### Piping and redirecting in Unix

#### Combining commands
```
# Perform multiple commands irrespective of whether previous command succeeds: combine commands via ';'

$ echo -n hello; echo -ne '\t'; false; echo -n world
hello	world

# Perform multiple commands only if previous command succeeds: combine commands via '&&'
# Note: The double-ampersand is the boolean AND operator in Unix.

$ echo -n hello && echo -ne '\t' && false && echo -n world
hello	
```

#### Passing and appending
```
# Passing output via '>'

$ echo hello > file1
$ cat file1
hello

# Appending output via '>>'

$ echo hello > file1
$ echo world >> file1
$ cat file1
hello
world
```
#### Passing output to another command (piping) via '|'
```
$ echo 'world hello' | tr h H
world Hello

$ echo 'world hello' | tr h H | tr ' ' '\n'
world
Hello

$ echo 'world hello' | tr h H | tr ' ' '\n' | tac
Hello
world

$ echo 'world hello' | tr h H | tr ' ' '\n' | tac | tr '\n' ' '
Hello world
```
