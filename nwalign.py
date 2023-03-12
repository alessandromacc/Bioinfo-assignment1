
class NWAligner:
    '''The NWAligner objects performs global alingment with the Needleman-Wunsch algorithm and allows the user to easily change
    the alignment scores. The alignment is performed within the initiator, therefore upon instancing of the class, and the user
    only needs to call the specific method for seeing the output of such alignment, self.getAlignment(). The user can also
    use self.getMatrix() to see the alignment matrix, but no other action is allowe, since the object takes care of everithing else automatically.'''
    def __init__(self, s1: str, s2: str, match_score: int = 1, mismatch_score: int = -1, gap_score: int = -2):
        '''Initializing the needed variables for the alignment'''
        self.s1 = s1.upper()
        self.s2 = s2.upper()
        self.dims = (len(s1), len(s2))
        self.rows = self.dims[0] + 1
        self.cols = self.dims[1] + 1
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_score = gap_score
        self.path = []
        self.paths = []
        self.__alignment = []

        '''Setting up the actual Matrix for the alignment, still empty'''
        self.alMat = []
        self.alMat.append([self.gap_score*i for i in range(self.cols)])
        for i in range(1, self.rows):
            row = []
            row.append(self.gap_score*i)
            for j in range(self.cols-1):
                row.append(None)
            self.alMat.append(row)
        
        '''Storing the value of the global alignment score and filling the matrix with fill() method'''
        self.score = self.__fill()

        self.m = self.alMat[::][::]
        
        '''Filling the boolean matrix for matches/mismatches'''
        self.match_matrix = [[False for i in range(self.cols)]]
        for i in range(self.rows-1):
            row = []
            row.append(False)
            for j in range(self.cols-1):
                row.append(self.s1[i] == self.s2[j])
            self.match_matrix.append(row)
        
        '''Once the previous steps are completed, pass to backtracking of the different eqiuvalent paths and graphical alignment for the terminal'''
        self.__backtrack()
        self.__recompose_paths()
        self.__align()

    
    def __fill(self) -> int:
        '''Fill the matrix automatically with the proper values calculated by __best_move();
        returns the value of the score of the global alignment.'''
        for i in range(1, self.rows):
            for j in range(1, self.cols):
                self.alMat[i][j] = self.__best_move((i,j))
        return self.alMat[i][j]
        
    def __best_move(self, coord):
        '''Determines the highest possible score for a matrix cell with given coordinates coord, according to the Needleman-Wunsch algorithm'''
        diag = self.alMat[coord[0]-1][coord[1]-1]
        if self.s1[coord[0]-1] == self.s2[coord[1]-1]:
            val = diag + self.match_score
        else:
            val = diag + self.mismatch_score
        vv = self.alMat[coord[0]-1][coord[1]] + self.gap_score
        hv = self.alMat[coord[0]][coord[1]-1] + self.gap_score
        return max(val, vv, hv)
    
    def __backtrack(self):
        '''Backtracks recursively along the matrix storing possible branchings of the main path as lists of coordinates'''
        #base case
        if (len(self.m), len(self.m[-1])) == (1,1):
            self.path.append((0,0))
        #horizontal only
        elif len(self.m) == 1:
            c = (len(self.m)-1, len(self.m[-1])-1)
            self.m = [[i for i in self.m[-1][:-1]]]
            self.match_matrix = [[i for i in self.match_matrix[-1][:-1]]]
            self.__backtrack()
            self.path.append(c)
        #vertical only
        elif len(self.m[-1]) == 1:
            c = (len(self.m)-1, len(self.m[-1])-1)
            self.m = self.m[:-1]
            self.match_matrix = self.match_matrix[:-1]
            self.__backtrack()
            self.path.append(c)
        else:
            #handling of general cases and more complex situations with branchings
            if (self.m[-1][-1] - self.m[-2][-2] == self.match_score and self.match_matrix[-1][-1] == True) or (self.m[-1][-1] - self.m[-2][-2] == self.mismatch_score and self.match_matrix[-1][-1] == False):
                current = self.m[::][::]
                current_match = self.match_matrix[::][::]
                if self.m[-1][-1] - self.m[-1][-2] != self.gap_score and self.m[-1][-1] - self.m[-2][-1] == self.gap_score:
                    c = (len(self.m)-1, len(self.m[-1])-1)
                    self.m = self.m[:-1]
                    self.match_matrix = self.match_matrix[:-1]
                    self.__backtrack()
                    self.path.append(c)
                    self.paths.append(self.path[::])
                    self.path = []
                elif self.m[-1][-1] - self.m[-1][-2] == self.gap_score and self.m[-1][-1] - self.m[-2][-1] != self.gap_score:
                    c = (len(self.m)-1, len(self.m[-1])-1)
                    self.m = [[i for i in j[:-1]] for j in self.m]
                    self.match_matrix = [[i for i in j[:-1]] for j in self.match_matrix]
                    self.__backtrack()
                    self.path.append(c)
                    self.paths.append(self.path[::])
                    self.path = []
                
                self.m = current
                self.match_matrix = current_match
                c = (len(self.m)-1, len(self.m[-1])-1)
                self.m = [[i for i in j[:-1]] for j in self.m[:-1]]
                self.match_matrix = [[i for i in j[:-1]] for j in self.match_matrix[:-1]]
                self.__backtrack()
                self.path.append(c)
            elif self.m[-1][-1] - self.m[-1][-2] != self.gap_score and self.m[-1][-1] - self.m[-2][-1] == self.gap_score:
                c = (len(self.m)-1, len(self.m[-1])-1)
                self.m = self.m[:-1]
                self.match_matrix = self.match_matrix[:-1]
                self.__backtrack()
                self.path.append(c)
            elif self.m[-1][-1] - self.m[-1][-2] == self.gap_score and self.m[-1][-1] - self.m[-2][-1] != self.gap_score:
                c = (len(self.m)-1, len(self.m[-1])-1)
                self.m = [[i for i in j[:-1]] for j in self.m]
                self.match_matrix = [[i for i in j[:-1]] for j in self.match_matrix]
                self.__backtrack()
                self.path.append(c)
            else:
                raise ValueError('Backtracking error')
    
    def __recompose_paths(self):
        '''Handles recompasition of the complete list of possible paths starting from the
        branches produced by the backtracking. Combinatorial procedure, it tests all possible combinations and
        stores as final paths only not-already-stored ones'''
        self.final_paths = [self.path[::]]
        self.paths.sort(reverse=True, key=len)
        for i in range(len(self.paths)):
            for j in range(len(self.paths[i])-1, -1, -1):
                for z in self.final_paths:
                    found = False
                    if self.paths[i][j] not in z:
                        found = True
                        self.paths[i] = self.paths[i][:j+1] + z[j+1:]
                        if self.paths[i] not in self.final_paths:
                            for k in range(1, len(self.paths[i])):
                                if self.paths[i][k][0] - self.paths[i][k-1][0] >= 0 and self.paths[i][k][1] - self.paths[i][k-1][1] >= 0:
                                    if k == len(self.paths[i])-1:
                                        self.final_paths.append(self.paths[i])
                                    else:
                                        continue
                                else:
                                    break
                if found:
                    break
    
    
    def __align(self):
        '''Last step of the alignment, moslty graphical, but contains some checking of the alignment score, since during benchmarking
        some cases produced final paths with different scores from the expected one. ierates over the list of paths and appends to
        self.__alignment the command line representation of an alignment, if it has passed the final check.'''
        for i in self.final_paths:
            al1 = []
            al2 = []
            for j in range(1, len(i)):
                if i[j][0] == i[j-1][0] and i[j][1] == i[j-1][1]:
                    raise ValueError('Something weird in alignment')
                elif i[j][0] == i[j-1][0]:
                    al1.append('-')
                    al2.append(self.s2[j-1-al2.count('-')])
                elif i[j][1] == i[j-1][1]:
                    al1.append(self.s1[j-1-al1.count('-')])
                    al2.append('-')
                else:
                    al1.append(self.s1[j-1-al1.count('-')])
                    al2.append(self.s2[j-1-al2.count('-')])
            matches = []
            for i in range(len(al1)):
                if al1[i] == al2[i]:
                    matches.append('*')
                elif al1[i] != al2[i] and (al1[i] != '-' and al2[i] != '-'):
                    matches.append('|')
                else:
                    matches.append(' ')
            validate_score = matches.count('*')*self.match_score + matches.count('|')*self.mismatch_score + matches.count(' ')*self.gap_score
            if self.score == validate_score:
                self.__alignment.append(' '.join(al1) + '\n' + ' '.join(matches) + '\n' + ' '.join(al2))
    
    
    def getAlignment(self):
        '''Outputs nicely all the produced alignment; this was implemented to be the only legal way for the use to access the alignment'''
        if len(self.__alignment) > 0:
            print('='*(len(max(self.__alignment, key=len))//3 + 10), end='\n\n')
            print('\t'+f'Global Score: {self.score}', end='\n\n')
            for i in range(len(self.__alignment)):
                print(f'[#{i+1}]', '-'*(len(max(self.__alignment, key=len))//3), f'[Global score: {self.score}]', end='\n\n')
                print(self.__alignment[i], end='\n\n')
            print('='*(len(max(self.__alignment, key=len))//3 + 10), end='\n\n')
        else:
            print('='*50, end='\n\n')
            print('No global alignment found, probably due to a mistake in the processing or bad parsing', end = '\n\n')
            print('='*50, end='\n\n')


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
