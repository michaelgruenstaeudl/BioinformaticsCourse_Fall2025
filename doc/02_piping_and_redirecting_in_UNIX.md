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
