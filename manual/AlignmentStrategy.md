## Smith Waterman Alignment
Adapter sequences are aligned to reads using the Smith-Waterman
Algorithm. For example, say we want to align the sequence *TACGTTGACACACGTC* against some
reference region: *AATATACGATGAACACCGTCATA*

The smith waterman alignment takes place in a few successive steps

1) score Matrix
The score matrix is computed using a variety of alignment parameters
* match = 1
* mismatch = -2
* insertion_open = -2
* insertion_extend = -1
* deletion_open = -1
* deletion_extend = -1
* minimum_alignment_score = 9

```
    -  T  A  C  G  T  T  G  A  C  A  C  A  C  G  T  C
-   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
A   0  0  1  0  0  0  0  0  1  0  1  0  1  0  0  0  0 
A   0  0  1  0  0  0  0  0  1  0  1  0  1  0  0  0  0 
T   0  1  0  0  0  1  1  0  0  0  0  0  0  0  0  1  0 
A   0  0  2  0  0  0  0  0  1  0  1  0  1  0  0  0  0 
T   0  1  1  0  0  1  1  0  0  0  0  0  0  0  0  1  0 
A   0  0  2  0  0  0  0  0  1  0  1  0  1  0  0  0  0 
C   0  0  1  3  1  0  0  0  0  2  0  2  0  2  0  0  1 
G   0  0  0  2  4  2  1  1  0  1  0  1  0  1  3  1  0 
A   0  0  1  1  3  2  1  0  2  0  2  0  2  0  2  1  0 
T   0  1  0  0  2  4  3  1  1  0  1  0  1  0  1  3  1 
G   0  0  0  0  1  3  2  4  2  1  0  0  0  0  1  2  1 
A   0  0  1  0  0  2  1  3  5  3  2  0  1  0  0  1  0 
A   0  0  1  0  0  1  0  2  4  3  4  2  1  0  0  0  0 
C   0  0  0  2  0  0  0  1  3  5  3  5  3  2  0  0  1 
A   0  0  1  1  0  0  0  0  2  4  6  4  6  4  3  2  1 
C   0  0  0  2  0  0  0  0  1  3  5  7  5  7  5  4  3 
C   0  0  0  1  0  0  0  0  0  2  4  6  5  6  5  4  5 
G   0  0  0  0  2  0  0  1  0  1  3  5  4  5  7  5  4 
T   0  1  0  0  1  3  1  0  0  0  2  4  3  4  6  8  6 
C   0  0  0  1  0  2  1  0  0  1  1  3  2  4  5  7  9 
A   0  0  1  0  0  1  0  0  1  0  2  2  4  3  4  6  8 
T   0  1  0  0  0  1  2  0  0  0  1  1  3  2  3  5  7 
A   0  0  2  0  0  0  1  0  1  0  1  0  2  1  2  4  6 

```

2) Traceback matrix

The traceback matrix is also computed at the same time the score matrix is being filled.
This matrix specifies the path taken for an alignment.
* -1 = diagonal (mismatch)
*  0  = diagonal (match)
*  1  = insertion
*  2  = deletion


```
    -  T  A  C  G  T  T  G  A  C  A  C  A  C  G  T  C
-   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
A   0  2  0  1  1  1  1  1  0  1  0  1  0  1  1  1  1 
A   0  2  0 -1  1  1  1  1  0 -1  0 -1  0 -1  1  1  1 
T   0  0  2 -1  1  0  0  1  2 -1  2 -1  2 -1  1  0  1 
A   0  2  0  1  1  2  2 -1  0  1  0  1  0  1  1  2 -1 
T   0  0  2 -1  1  0  0  1  2 -1  2 -1  2 -1  1  0  1 
A   0  2  0  1  1  2  2 -1  0  1  0  1  0  1  1  2 -1 
C   0  2  2  0  1  1  1  1  2  0  1  0  1  0  1  1  0 
G   0  2  2  2  0  1  1  0  1  2 -1  2 -1  2  0  1  1 
A   0  2  0  2  2 -1  1  1  0  1  0  1  0  1  2 -1  1 
T   0  0  2  2  2  0  0  1  2 -1  2 -1  2 -1  2  0  1 
G   0  2 -1  1  0  2 -1  0  1  1  1 -1  2 -1  0  2 -1 
A   0  2  0  1  2  2 -1  2  0  1  0  1  0  1  2  2 -1 
A   0  2  0 -1  1  2 -1  2  0 -1  0  1  0 -1  1  2 -1 
C   0  2  2  0  1  2 -1  2  2  0  1  0  1  0  1  1  0 
A   0  2  0  2 -1  1  1  2  0  2  0  1  0  1  1  1  1 
C   0  2  2  0  1  1  1  1  2  0  2  0  1  0  1  1  0 
C   0  2  2  0 -1  1  1  1  2  0  2  0 -1  0 -1  1  0 
G   0  2  2  2  0  1  1  0  1  2  2  2 -1  2  0  1  1 
T   0  0  1  1  2  0  0  2 -1  2  2  2 -1  2  2  0  1 
C   0  2 -1  0  2  2 -1  1  1  0  2  0 -1  0  2  2  0 
A   0  2  0  2 -1  2 -1 -1  0  2  0  2  0  2  2  2  2 
T   0  0  2 -1  1  0  0  1  2 -1  2  2  2 -1  2  0  2 
A   0  2  0  1  1  2  2 -1  0  1  0  2  0 -1  2  2  2 

```

3) Maximum Alignment Score
The maximum alignment score is the best possible alignment between the two strings given the scoring parameters. Using the score matrix above we can find the reference position (by row) and query position (by column) yielding the best alignment.  The last parameter

**a) Score array**
maximum alignment score by reference position
```
                                                          
A  A  T  A  T  A  C  G  A  T  G  A  A  C  A  C  C  G  T  C  A  T  A
1  1  1  2  1  2  3  4  3  4  4  5  4  5  6  7  6  7  8  9  8  7  6 
                                                         ^
```

position 20 in the reference (see ^), the query sequenced yields a score of 9.

**b) Score index**
query position yielding maximum alignment score for each reference position
```
A  A  T  A  T  A  C  G  A  T  G  A  A  C  A  C  C  G  T  C  A  T  A
2  2  1  2  1  2  3  4  4  5  7  8  8  9  10 11 11 14 15 16 16 16 16
                                                         ^
```
At position 20 in the reference, we find the maximum score to be position 16 of the query.

4) Trace the alignment 
The maximum alignment score told us the end of our alignment (query_end = 16, reference_end = 20). Now, using the traceback matrix, we can trace the alignment between the two strings
Following the steps in the traceback matrix, continuing until the score is zero or we reach the end of either sequence.

```
TACGTTG-ACACACGTC
||||*|| |||| ||||
TACGATGAACAC-CGTC
```

