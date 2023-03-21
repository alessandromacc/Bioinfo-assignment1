from nwalign import NWAligner
from swalign import SWAligner
import sys

'''Create an NWAligner object, parsing mandatorily two strings as the two sequence to be aligned
and optionally the values for the alignment scores, as the object definition shows.
The alignment is carried out upon class instancing.
Using Python sys library the test file requires the user to pass the sequences as inputs at file call as two distinct strings;
Moreover, the user can use some available flags to change at runtime the matrix scores for match/mismatch/gap and the treshold:
    "mNW\int" for changing the match score in global alignment;
    "msNW\int" for changing the mismatch score in global alignment;
    "gapNW\int" for changing the gap score in global alignment;
    "mSW\int" for changing the match score in local alignment;
    "msSW\int" for changing the mismatch score in local alignment;
    "gapSW\int" for changing the gap score in local alignment;
    "trSW\int" for changing the treshold for local alignment;
    "metNW\str" for swithcing between recursive and iterative execution in global alignment: it should specify "it" for iterative, and "rec" for recursive. 
If not differently specified, the default values will be used, which include iterative procedure for global alignment.'''

#setting the standard parameters
NW_match_score = 1
NW_mismatch_score = -1
NW_gap_score = -2
SW_match_score = 1
SW_mismatch_score = -1
SW_gap_score = -2
threshold = 1
method = 'iterative'

#updating the standard parameters with user indications
seqs = []
for i in sys.argv[1:]:
    if 'mNW' in i:
        NW_match_score = int(i[3:])
    if 'mSW' in i:
        SW_match_score = int(i[3:])
    if 'msNW' in i:
        NW_mismatch_score = int(i[4:])
    if 'msSW' in i:
        SW_mismatch_score = int(i[4:])
    if 'gapNW' in i:
        NW_gap_score = int(i[5:])
    if 'gapSW' in i:
        SW_gap_score = int(i[5:])
    if 'thrSW' in i:
        threshold = int(i[5:])
    if 'metNW' in i:
        method = i[5:]
        if method == 'it':
            method = 'iterative'
        elif method == 'rec':
            method = 'recursive'
        else:
            print('Bad flag inserted for global alignment method: defaulting to iterative algorithm')
            method = 'iterative'
        
    if 'mNW' not in i and 'msNW' not in i and 'gapNW' not in i and 'thrSW' not in i and 'metNW' not in i and 'mSW' not in i and 'msSW' not in i and 'gapSW' not in i:
        seqs.append(i)

if len(seqs) == 2:
    z = NWAligner(seqs[0], seqs[1], NW_match_score, NW_mismatch_score, NW_gap_score, method=method)
    #The method NWAligner.getAlignment() outputs in the terminal the result of the performed global alignment
    z.getAlignment()

    x = SWAligner(seqs[0], seqs[1], SW_match_score, SW_mismatch_score, SW_gap_score,threshold=threshold)
    #The method SWAligner.getAlignment() outputs in the terminal the result of the performed global alignment
    x.getAlignment()
else:
    raise(ValueError('ArgumentError: you must pass exactly two strings as input'))