# First Bioinformatics assignment:
The swalign.py file and nwalign.py file contain the code for implementing respectively Smith-Waterman and Needleman-Wunsch alignment algorithm.
Smith-Waterman only comes in a recursive implementation, with linear gap penalty, while Needleman-Wunsch also has linear gap penalty, but
provides the user with two alternative backtracking strategies: iterative and recursive, which can be easily switched when writing code, and
even more easily in test.py, which supports specific flags to change the method. The double implementation covers some issues emerged with
the recursive version testing: the number of recursive calls Python can make is limited by its high memory usage, but performs well with shorter sequences
if compared to an iterative implementation. When working with longer sequences, the iterative version shows better performances, and does not
technically have a limit as that provided by recursion.

The expected way of usage is already shown in the test.py file, which uses Python "sys" library to allow the user to provide the sequences to align
next to the file call in the terminal. test.py aligns the sequence pair both glabally  and locally, but the NWAligner class and the SWAligner class
can be used with different interfaces. Furthermore, test.py allows the user to change score settings for both algorithms from the call, by inserting
specific flags which are illustrated directly inside test.py.

No external libraries are required, as everything has been implemented in core Python, although this means renouncing to some efficient
data structures as those provided by scientific Python libraries.

The outputs from my personal testing were compared to the results produced by the alignment tools provided at https://rna.informatik.uni-freiburg.de/Teaching/,
and all the algorithms only compute the first 20 equivalent alignments, for computational efficiency's sake.
