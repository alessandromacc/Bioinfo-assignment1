
class SWAligner:
    '''The SWAligner objects performs local alingment with the Smith-Waterman algorithm and allows the user to easily change
    the alignment scores and the treshold for selecting minimum overlap length. The alignment is performed within the initiator,
    therefore upon instancing of the class, and the user only needs to call the specific method for seeing the output of such alignment, self.getAlignment().
    The user can also use self.getMatrix() to see the alignment matrix, but no other action is allowe, since the object takes care of everithing else automatically.'''
    def __init__(self, s1: str, s2: str, match_score: int = 1, mismatch_score: int = -1, gap_score: int = -2, treshold: int = 1):
        '''Initializing the needed variables for the alignment'''
        self.s1 = s1.upper()
        self.s2 = s2.upper()
        self.dims = (len(s1), len(s2))
        self.rows = self.dims[0] + 1
        self.cols = self.dims[1] + 1
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_score = gap_score
        self.treshold = treshold
        self.path = []
        self.paths = []
        self.alignment = []

        '''Setting up the alignment matrix, still empty'''
        self.alMat = []
        self.alMat.append([0 for i in range(self.cols)])
        for i in range(1, self.rows):
            row = []
            row.append(0)
            for j in range(self.cols-1):
                row.append(None)
            self.alMat.append(row)
        
        '''Setting up the boolena matrix for matches/mismatches'''
        self.match_matrix = [[False for i in range(self.cols)]]
        for i in range(self.rows-1):
            row = []
            row.append(False)
            for j in range(self.cols-1):
                row.append(self.s1[i] == self.s2[j])
            self.match_matrix.append(row)

        '''Filling the alignment matrix with the proper values and determining the possible local alignment starting points as a list of matrix coordinates'''
        self.__fill()
        
        '''Necessary setup for the backtracking and storing paths in self.paths'''
        for i in self.starts:
            self.path = []
            self.am = [[k for k in j[:i[1]+1]] for j in self.alMat[:i[0]+1]]
            self.mm = [[k for k in j[:i[1]+1]] for j in self.match_matrix[:i[0]+1]]
            self.__backtrack()
            if len(self.path) >= self.treshold:
                self.paths.append(self.path[::-1])
        
        '''Final alignment, one for each path produced'''
        self.__align()


    def __fill(self):
        '''Along with filling the matrix with the proper values from the __best_move() method,
        takes care of finding the starting points with the highest scores and storing them properly for later use'''
        max = 0
        max_coord = []
        for i in range(1, self.rows):
            for j in range(1, self.cols):
                self.alMat[i][j] = self.__best_move((i,j))
                if self.alMat[i][j] > max:
                    max_coord = []
                    max = self.alMat[i][j]
                    max_coord.append((i,j))
                elif self.alMat[i][j] == max:
                    max_coord.append((i,j))
                else:
                    continue
        self.score = max
        self.starts = max_coord

    
    def __best_move(self, coord):
        '''Determines the highest possible score for a matrix cell with given coordinates coord, according to the Smith-Waterman algorithm'''
        diag = self.alMat[coord[0]-1][coord[1]-1]
        if self.s1[coord[0]-1] == self.s2[coord[1]-1]:
            val = diag + self.match_score
        else:
            val = diag + self.mismatch_score
        vv = self.alMat[coord[0]-1][coord[1]] + self.gap_score
        hv = self.alMat[coord[0]][coord[1]-1] + self.gap_score
        return max(val, vv, hv, 0)
    
    def __backtrack(self):
        '''Backtracks recursively along the matrix to identify the alignment path associated with a start;
        quite complex in its case handling, since the aim was that of covering the most cases where a precise decision has to be taken,
        for the sake of the algorithmic correctness. Since it is used iteratively by the object over many paths, employs copies of the main
        alignment matrix and match matrix to cut recursively and carry along in the backtracking'''
        if len(self.am) == 1:
            if self.am[-1][-2] == 0:
                self.path.append((len(self.am)-1, len(self.am[-1])-1))
                return
            else:
                self.path.append((len(self.am)-1, len(self.am[-1])-1))
                self.am = [[i for i in self.am[-1][:-1]]]
                self.mm = [[i for i in self.mm[-1][:-1]]]
                self.__backtrack()
        elif len(self.am[-1]) == 1:
            if self.am[-2][-1] == 0:
                self.path.append((len(self.am)-1, len(self.am[-1])-1))
                return
            else:
                self.path.append((len(self.am)-1, len(self.am[-1])-1))
                self.am = self.am[:-1]
                self.mm = self.mm[:-1]
                self.__backtrack()
        elif self.am[-1][-2] == self.am[-2][-1] == self.am[-2][-2] == 0:
            self.path.append((len(self.am)-1, len(self.am[-1])-1))
            return
        elif self.am[-1][-2] == self.am[-2][-1] == 0:
            if (self.am[-1][-1] - self.am[-2][-2] == self.match_score and self.mm[-1][-1] == True) or (self.am[-1][-1] - self.am[-2][-2] == self.mismatch_score and self.mm[-1][-1] == False):
                self.path.append((len(self.am)-1, len(self.am[-1])-1))
                self.am = [[k for k in j[:-1]] for j in self.am[:-1]]
                self.mm = [[k for k in j[:-1]] for j in self.mm[:-1]]
                self.__backtrack()
            else:
                raise(ValueError('BacktrackingError: attempted diagonal movement'))
        elif self.am[-1][-2] == self.am[-2][-2] == 0:
            if self.am[-1][-1] - self.am[-2][-1] == self.gap_score:
                self.path.append((len(self.am)-1, len(self.am[-1])-1))
                self.am = self.am[:-1]
                self.mm = self.mm[:-1]
                self.__backtrack()
            else:
                self.path.append((len(self.am)-1, len(self.am[-1])-1))
                return
        elif self.am[-2][-1] == self.am[-2][-2] == 0:
            if self.am[-1][-1] - self.am[-1][-2] == self.gap_score:
                self.path.append((len(self.am)-1, len(self.am[-1])-1))
                self.am = [[k for k in j[:-1]] for j in self.am]
                self.mm = [[k for k in j[:-1]] for j in self.mm]
                self.__backtrack()
            else:
                self.path.append((len(self.am)-1, len(self.am[-1])-1))
                return
        elif self.am[-1][-2] == 0:
            if (self.am[-1][-1] - self.am[-2][-2] == self.match_score and self.mm[-1][-1] == True) or (self.am[-1][-1] - self.am[-2][-2] == self.mismatch_score and self.mm[-1][-1] == False) and (self.am[-2][-2] >= self.am[-2][-1]):
                self.path.append((len(self.am)-1, len(self.am[-1])-1))
                self.am = [[k for k in j[:-1]] for j in self.am[:-1]]
                self.mm = [[k for k in j[:-1]] for j in self.mm[:-1]]
                self.__backtrack()
        elif self.am[-2][-1] == 0:
            if (self.am[-1][-1] - self.am[-2][-2] == self.match_score and self.mm[-1][-1] == True) or (self.am[-1][-1] - self.am[-2][-2] == self.mismatch_score and self.mm[-1][-1] == False) and (self.am[-2][-2] >= self.am[-1][-2]):
                self.path.append((len(self.am)-1, len(self.am[-1])-1))
                self.am = [[k for k in j[:-1]] for j in self.am[:-1]]
                self.mm = [[k for k in j[:-1]] for j in self.mm[:-1]]
                self.__backtrack()
        elif self.am[-2][-2] == 0:
            if self.am[-1][-1] - self.am[-1][-2] != self.gap_score and self.am[-1][-1] - self.am[-2][-1] != self.gap_score:
                raise(ValueError('BacktrackingError: ambiguous NO-DIAGONAL situation - 1'))
            elif self.am[-1][-1] - self.am[-1][-2] == self.gap_score and self.am[-1][-1] - self.am[-2][-1] != self.gap_score:
                self.path.append((len(self.am)-1, len(self.am[-1])-1))
                self.am = [[k for k in j[:-1]] for j in self.am]
                self.mm = [[k for k in j[:-1]] for j in self.mm]
                self.__backtrack()
            elif self.am[-1][-1] - self.am[-1][-2] != self.gap_score and self.am[-1][-1] - self.am[-2][-1] == self.gap_score:
                self.path.append((len(self.am)-1, len(self.am[-1])-1))
                self.am = self.am[:-1]
                self.mm = self.mm[:-1]
                self.__backtrack()
            else:
                raise(ValueError('BacktrackingError: ambiguous NO-DIAGONAL situation - 2'))
        else:
            if (self.am[-1][-1] - self.am[-2][-2] == self.match_score and self.mm[-1][-1] == True) or (self.am[-1][-1] - self.am[-2][-2] == self.mismatch_score and self.mm[-1][-1] == False):
                self.path.append((len(self.am)-1, len(self.am[-1])-1))
                self.am = [[k for k in j[:-1]] for j in self.am[:-1]]
                self.mm = [[k for k in j[:-1]] for j in self.mm[:-1]]
                self.__backtrack()
            else:
                raise(ValueError('BacktrackingError: attempted diagonal movement'))
    
    def __align(self):
        '''Final alingment executed automatically in an iterative way over the paths in self.paths; 
        appends to self.alignment the graphical version of the aligned sequences to be output. As in common
        depictions of Smith-Waterman aligned sequences, non-ovelaying parts are included.'''
        for i in self.paths:
            al1 = []
            al2 = []
            matches = []
            al1.append(self.s1[i[0][0]-1])
            al2.append(self.s2[i[0][1]-1])
            for j in range(1, len(i)):
                if i[j][0] == i[j-1][0] and i[j][1] == i[j-1][1]:
                    raise(ValueError('AlignemntError - 1'))
                elif i[j][0] == i[j-1][0]:
                    al1.append('-')
                    al2.append(self.s2[i[j-1][1]])
                elif i[j][1] == i[j-1][1]:
                    al1.append(self.s1[i[j-1][0]])
                    al2.append('-')
                else:
                    al1.append(self.s1[i[j-1][0]])
                    al2.append(self.s2[i[j-1][1]])
            
            ls = len(max(self.s1[:i[0][0]-1], self.s2[:i[0][1]-1], key=len))
            rs = len(max(self.s1[i[-1][0]:], self.s2[i[-1][1]:], key=len))
            
            al1 = [' ' for i in range(ls-len(self.s1[:i[0][0]-1]))] + list(self.s1[:i[0][0]-1]) + al1 + list(self.s1[i[-1][0]:]) + [' ' for i in range(rs - len(self.s1[i[-1][0]:]))]
            al2 = [' ' for i in range(ls-len(self.s2[:i[0][1]-1]))] + list(self.s2[:i[0][1]-1]) + al2 + list(self.s2[i[-1][1]:]) + [' ' for i in range(rs - len(self.s2[i[-1][1]:]))]
            
            for i in range(len(al1)):
                if al1[i] == ' ' or al2[i] == ' ':
                    matches.append(' ')
                elif al1[i] != al2[i]:
                    matches.append('|')
                else:
                    matches.append('*')
            al = (' '.join(al1) + '\n' + ' '.join(matches) + '\n' + ' '.join(al2))
            if al not in self.alignment:
                self.alignment.append(al)
            
    def getMatrix(self):
        '''Outputs nicely the alignment matrix, useful in testing and checking.'''
        print('', end='\t')
        print('§', end = '\t')
        for i in range(len(self.s2)):
            print(self.s2[i], end = '\t')
        print('\n')
        for i in range(len(self.alMat)):
            if i != 0:
                print(self.s1[i-1], end='\t')
            else:
                print('§', end='\t')
            for j in range(len(self.alMat[i])):
                print(self.alMat[i][j], end = '\t')
            print('\n')
    
    def getAlignment(self):
        '''Outputs nicely all the produced alignment; this was implemented to be the only legal way for the use to access the alignment'''
        if len(self.alignment) > 0:
            print('='*(len(max(self.alignment, key=len))//3 + 10), end='\n\n')
            print('\t'+f'Local Score: {self.score}', end='\n\n')
            for i in range(len(self.alignment)):
                print(f'[#{i+1}]', '-'*(len(max(self.alignment, key=len))//3), f'[Local score:{self.score}]', end='\n\n')
                print(self.alignment[i], end='\n\n')
            print('='*(len(max(self.alignment, key=len))//3 + 10), end='\n\n')
        else:
            print('='*50, end='\n\n')
            print('No local alignment found, have you tried lowering the selection treshold yet?', end='\n\n')
            print('='*50, end='\n\n')