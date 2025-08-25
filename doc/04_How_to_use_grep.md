### Introducing grep
#### Text file filtering
```
grep    # Global search for Regular Expression and Print the results
        # i.e., printing lines that match patterns
```
#### General usage
```
grep "PATTERN" <file_to_be_filtered>
```


#### Setting up test file
```
cd ~
echo "This is the first line." > file1.txt
echo "This is the second line." >> file1.txt
echo "This is the third line." >> file1.txt
echo "And here is the fourth line." >> file1.txt
echo "Fifth one and the end." >> file1.txt
```

#### Using grep on test file
```
grep "line" file1.txt    # print all lines with keyword 'line'

grep --color "line" file1.txt    # print all lines with keyword 'line'
                                 # while highlighting the keyword
```

#### Using grep with different options - part 1
```
grep -v "line" file1.txt    # print all lines without keyword 'line'

grep -v "line" file1.txt | grep "here"  # print all lines with 
                                        # keyword 'here' that 
                                        # do not contain keyword 'line'

grep 'line\|end' file1.txt  # print all lines with either keyword 
                            # 'line' or keyword 'end'
                            
grep "^This" file1.txt    # print all lines that start 
                          # with keyword 'This' (case-sensitive)
```

#### Using grep with different options - part 2
```
grep -A1 "second" file1.txt   # print all lines with keyword 'second' and 
                              # any line immediately following it

grep -B2 "fourth" file1.txt  # print all lines with keyword 'fourth' and 
                             # the two lines immediately before it
```

#### Using grep with regular expressions
```
grep "This .* line" file1.txt  # print all lines that start 
                               # with 'This' and end with 'line'

grep "This is the ..... line" file1.txt  # print all lines that 
                                         # start with 'This is the',
                                         # end with 'line', and has
                                         # exactly five characters
                                         # in between.

grep -o "f...t" file1.txt  # from all lines, print only the 5-character
                           # word that starts with 'f' and ends with 't'
                           
grep -o ".* line\.$" file1.txt  # print all lines that end with 'line.'
```
