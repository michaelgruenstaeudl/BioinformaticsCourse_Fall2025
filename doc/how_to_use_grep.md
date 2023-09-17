### Introducing grep
```
# Text file filtering
grep    # Global search for Regular Expression and Print the results
        # i.e., printing lines that match patterns

# General usage
grep "PATTERN" <file_to_be_filtered>
```


```
# Setting up test file
cd ~
echo "This is the first line." > file1.txt
echo "This is the second line." >> file1.txt
echo "This is the third line." >> file1.txt
echo "And here is the fourth line." >> file1.txt
echo "Fifth one and the end." >> file1.txt

# Using grep on test file
grep "line" file1.txt    # print all lines with keyword 'line'

grep --color "line" file1.txt    # print all lines with keyword 'line'
                                 # while highlighting the keyword
```

```
# Using grep on test file
grep -v "line" file1.txt    # print all lines without keyword 'line'

grep -v "line" file1.txt | grep "here"  # print all lines with 
                                        # keyword 'here' that 
                                        # do not contain keyword 'line'

grep 'line\|end' file1.txt  # print all lines with either keyword 
                            # 'line' or keyword 'end'
                            
grep "^This" file1.txt    # print all lines that start 
                          # with keyword 'This' (case-sensitive)
```
