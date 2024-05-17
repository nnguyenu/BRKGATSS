This folder contains experiments result of 4 algorithms: **BRKGA**, **BRKGAreverse**, **fastBRKGA**, **fastBRKGAreverse** on 24 graphs.

For 'BRKGA.txt' and 'fastBRKGA.txt' file, each line of a file denotes a new best solution found by the algorithm and has the format:
```
current_time generation result list_of_selected_nodes
```
For example, the first line of 'jazz/BRKGA.txt' is:
```
0.00467999 0 29 6 31 34 39 48 59 64 68 76 80 85 97 99 100 105 108 110 127 130 131 134 135 148 166 167 169 170 191 194 
```
This means that the algorithm found a new best solution at 0.00467999 seconds, generation 0, with 29 nodes selected. The next 29 numbers are the list of selected nodes.

There are 10 runs for each algorithm & graph (though some graphs may have more runs). The best result of each run has the format:
```
Finish best_time best_result generation list_of_selected_nodes
```
For example, the result line of the first run of 'jazz/BRKGA.txt' is:
```
Finish 1.39932 20 582 6 53 59 68 69 82 85 95 113 121 129 131 135 148 149 157 166 167 171 191 
```
This means that the algorithm got its best result at 1.39932 seconds and generation 582, with 20 nodes selected. The next 20 numbers are the list of selected nodes.

For 'BRKGAreverse.txt' and 'fastBRKGAreverse.txt' file, each line of a file denotes the new solution after applied the reverse greedy heuristic. It has the format:
```
result gereration list_of_selected_nodes
```
For example, the first line of 'jazz/BRKGAreverse.txt' is:
```
20 6 53 59 68 69 82 85 95 113 121 129 131 135 148 149 157 166 167 171 191
```
This means that the reverse greedy heuristic found a solution with 20 nodes selected. The next 20 numbers are the list of selected nodes. The original solution is the best solution found by **BRKGA**.