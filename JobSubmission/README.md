## How to submit analysis jobs in virgo cluster  
1. `ssh -x <username>@lx-pool.gsi.de`    
2. `ssh -x virgo.hpc.gsi.de`  
3. Say you have a list of root files named `root_list.txt`.  
You can split this list into a number `n` of smaller lists using `split` shell command. For example  
`split -l 50 --numeric-suffixes root_list.txt small_root_list`. THis will create smaller lists named `small_root_list_xx`  
where `xx` is numueric prefix.  
4. Create a list of lists using wildcards `ls small_root_list_* > list_of_lists.txt`  
5. Finally you can submit your jobs by typing `. sendScript.sh <n> list_of_lists.txt <out_prefix>`, where `<out_prefix>` is the name of the output root file (in case you save some plots in an output root file) and `<n>` is the number of jobs, this MUST be equal to the number of lines in `list_of_lists.txt`. You can get the number of lines (corresponding to the number of inner lists) by simply typing in your shell `wc -l list_of_lists.txt`.  
