from nwalign import NWAligner
from swalign import SWAligner
import sys

'''Create an NWAligner object, parsing mandatorily two strings as the two sequence to be aligned
and optionally the values for the alignment scores, as the object definition shows.
The alignment is carried out upon class instancing.
Using Python sys library the test file requires the user to pass the sequences as inputs at file call as two distinct strings;
Moreover, the user can use some available flags to change at runtime the matrix scores for match/mismatch/gap and the treshold:
    "m\int" for changing the match score;
    "ms\int" for changing the mismatch score;
    "gap\int" for changing the gap score;
    "tr\int" for changing the treshold for local alignment.
If not differently specified, the default values will be used.'''

#setting the standard parameters
match_score = 1
mismatch_score = -1
gap_score = -2
treshold = 1

#updating the standard parameters with user indications
seqs = []
for i in sys.argv[1:]:
    if 'm' in i:
        match_score = int(i[1:])
    if 'ms' in i:
        mismatch_score = int(i[2:])
    if 'gap' in i:
        gap_score = int(i[3:])
    if 'tr' in i:
        treshold = int(i[2:])
    if 'm' not in i and 'ms' not in i and 'gap' not in i and 'tr' not in i:
        seqs.append(i)

if len(seqs) == 2:
    z = NWAligner(seqs[0], seqs[1], match_score, mismatch_score, gap_score)
    #The method NWAligner.getAlignment() outputs in the terminal the result of the performed global alignment
    z.getAlignment()

    x = SWAligner(seqs[0], seqs[1], match_score, mismatch_score, gap_score,treshold=1)
    #The method SWAligner.getAlignment() outputs in the terminal the result of the performed global alignment
    x.getAlignment()
else:
    raise(ValueError('ArgumentError: you must pass exactly two strings as input'))