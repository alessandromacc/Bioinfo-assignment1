
class NWAligner:
    '''The NWAligner objects performs global alingment with the Needleman-Wunsch algorithm and allows the user to easily change
    the alignment scores. The alignment is performed within the initiator, therefore upon instancing of the class, and the user
    only needs to call the specific method for seeing the output of such alignment, self.getAlignment(), with default word-wrap of 75 characters. The user can also
    use self.getMatrix() to see the alignment matrix, but no other action is allowe, since the object takes care of everithing else automatically.
    The user can also select which method is to be used in backtracking: recursive or iterative, by specifying the "method" argument;
    whichever is specified will be used, and the user cannot furthermore interfere with the internal processes.'''
    def __init__(self, s1: str, s2: str, match_score: int = 1, mismatch_score: int = -1, gap_score: int = -2, method: str = 'iterative'):
        '''Initializing the needed variables for the alignment'''
        self.s1 = s1.upper()
        self.s2 = s2.upper()
        self.dims = (len(s1), len(s2))
        self.rows = self.dims[0] + 1
        self.cols = self.dims[1] + 1
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_score = gap_score
        self.method = method
        self.path = []
        self.paths = []
        self.__alignment = []
        self.lenError = False
        if len(self.s1) >1 and len(s2) > 1:
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

            '''Once the previous steps are completed, pass to backtracking of the different eqiuvalent paths and graphical alignment for the terminal,
            according to the requested method'''
            if self.method == 'iterative':
                self.__backtrack_it()
            elif self.method == 'recursive':
                self.__backtrack()
            else:
                raise ValueError(f'method argument must be either the string "iterative" or "recursive"; Received: {self.method}')
            self.__align()
        else:
            self.lenError = True

    
    def __fill(self) -> int:
        '''Fill the matrix automatically with the proper values calculated by __best_move();
        returns the value of the score of the global alignment.'''
        for i in range(1, self.rows):
            for j in range(1, self.cols):
                self.alMat[i][j] = self.__best_move((i,j))
        return self.alMat[i][j]
        
    def __best_move(self, coord: tuple) -> int:
        '''Determines the highest possible score for a matrix cell with given coordinates coord, according to the Needleman-Wunsch algorithm'''
        diag = self.alMat[coord[0]-1][coord[1]-1]
        if self.s1[coord[0]-1] == self.s2[coord[1]-1]:
            val = diag + self.match_score
        else:
            val = diag + self.mismatch_score
        vv = self.alMat[coord[0]-1][coord[1]] + self.gap_score
        hv = self.alMat[coord[0]][coord[1]-1] + self.gap_score
        return max(val, vv, hv)
    
    def __backtrack(self) -> None:
        '''Backtracks recursively along the matrix storing the first 20 alignment paths as lists of coordinates'''
        #base case
        if (len(self.m), len(self.m[-1])) == (1,1):
            self.path.append((0,0))
            if len(self.paths) < 20:
                self.paths.append(self.path[::-1])
                self.path = []
            return
        #horizontal only
        elif len(self.m) == 1:
            c = (len(self.m)-1, len(self.m[-1])-1)
            self.m = [[i for i in self.m[-1][:-1]]]
            self.match_matrix = [[i for i in self.match_matrix[-1][:-1]]]
            if len(self.paths) < 20:
                self.path.append(c)
                self.__backtrack()
        #vertical only
        elif len(self.m[-1]) == 1:
            c = (len(self.m)-1, len(self.m[-1])-1)
            self.m = self.m[:-1]
            self.match_matrix = self.match_matrix[:-1]
            if len(self.paths) < 20:
                self.path.append(c)
                self.__backtrack()
        else:
            #handling of general cases and more complex situations with branchings
            if (self.m[-1][-1] - self.m[-2][-2] == self.match_score and self.match_matrix[-1][-1] == True) or (self.m[-1][-1] - self.m[-2][-2] == self.mismatch_score and self.match_matrix[-1][-1] == False):
                current = self.m[::][::]
                current_match = self.match_matrix[::][::]
                current_path = self.path[::]
                if self.m[-1][-1] - self.m[-2][-1] == self.gap_score:
                    c = (len(self.m)-1, len(self.m[-1])-1)
                    self.m = self.m[:-1]
                    self.match_matrix = self.match_matrix[:-1]
                    if len(self.paths) < 20:
                        self.path.append(c)
                        self.__backtrack()
                    self.m = current
                    self.match_matrix = current_match
                    self.path = current_path[::]

                if self.m[-1][-1] - self.m[-1][-2] == self.gap_score:
                    c = (len(self.m)-1, len(self.m[-1])-1)
                    self.m = [[i for i in j[:-1]] for j in self.m]
                    self.match_matrix = [[i for i in j[:-1]] for j in self.match_matrix]
                    if len(self.paths) < 20:
                        self.path.append(c)
                        self.__backtrack()
                
                self.path = current_path[::]
                self.m = current
                self.match_matrix = current_match
                c = (len(self.m)-1, len(self.m[-1])-1)
                self.m = [[i for i in j[:-1]] for j in self.m[:-1]]
                self.match_matrix = [[i for i in j[:-1]] for j in self.match_matrix[:-1]]
                if len(self.paths) < 20:
                    self.path.append(c)
                    self.__backtrack()
            elif self.m[-1][-1] - self.m[-2][-1] == self.gap_score:
                current = self.m[::][::]
                current_match = self.match_matrix[::][::]
                current_path = self.path[::]
                if self.m[-1][-1] - self.m[-1][-2] == self.gap_score:
                    c = (len(self.m)-1, len(self.m[-1])-1)
                    self.m = [[i for i in j[:-1]] for j in self.m]
                    self.match_matrix = [[i for i in j[:-1]] for j in self.match_matrix]
                    if len(self.paths) < 20:
                        self.path.append(c)
                        self.__backtrack()
                    
                self.path = current_path
                self.m = current
                self.match_matrix = current_match
                c = (len(self.m)-1, len(self.m[-1])-1)
                self.m = self.m[:-1]
                self.match_matrix = self.match_matrix[:-1]
                if len(self.paths) < 20:
                    self.path.append(c)
                    self.__backtrack()
            elif self.m[-1][-1] - self.m[-1][-2] == self.gap_score:
                current = self.m[::][::]
                current_match = self.match_matrix[::][::]
                current_path = self.path[::]
                if self.m[-1][-1] - self.m[-2][-1] == self.gap_score:
                    c = (len(self.m)-1, len(self.m[-1])-1)
                    self.m = self.m[:-1]
                    self.match_matrix = self.match_matrix[:-1]
                    if len(self.paths) < 20:
                        self.path.append(c)
                        self.__backtrack()
                
                self.path = current_path
                self.m = current
                self.match_matrix = current_match
                c = (len(self.m)-1, len(self.m[-1])-1)
                self.m = [[i for i in j[:-1]] for j in self.m]
                self.match_matrix = [[i for i in j[:-1]] for j in self.match_matrix]
                if len(self.paths) < 20:
                    self.path.append(c)
                    self.__backtrack()
            else:
                raise ValueError('Backtracking error')
    
    def __backtrack_it(self) -> None:
        '''Backtracks iteratively to find the first 20 equivalent alignment paths and stores them as lists of coordinates'''
        self.pointer = (len(self.m)-1, len(self.m[-1])-1)
        self.branches = []
        self.solved_branches = []
        self.current_branch = (None, None, None)
        self.directive = None
        done = False

        while not done:
            if len(self.paths) > 19:
                break
            branch_index = [k[0] for k in self.branches]
            solved_index = [k[0] for k in self.solved_branches]
            if self.pointer == (0,0):
                #termination of path detection
                self.path.append((0,0))
                if self.path[::-1] not in self.paths:
                    self.paths.append(self.path[::-1])
                #print(self.paths)
                self.directive = None

                #redirecting block
                if self.branches != []:
                    if self.current_branch in self.branches:
                        self.branches.remove(self.current_branch)
                        self.solved_branches.append(((self.current_branch[0][0],self.current_branch[0][1]),self.current_branch[1],[(i[0],i[1]) for i in self.current_branch[2]]))
                        #recompute easy access branch infos
                        branch_index = [k[0] for k in self.branches]
                        solved_index = [k[0] for k in self.solved_branches]
                    if self.branches != []:
                        minimum = max([i[0]+i[1] for i in branch_index])
                        summed_indexes = [i[0]+i[1] for i in branch_index]
                        min_index = summed_indexes.index(minimum)
                        closest_branch = (branch_index[min_index], self.branches[min_index][1], [(i[0], i[1]) for i in self.branches[min_index][2]])
                        self.pointer = closest_branch[0]
                        self.path = closest_branch[2]
                        self.directive = closest_branch[1]
                        #print(self.branches)
                    else:
                        done = True
                else:
                    done = True
                
            # horizontal only
            elif self.pointer[0] == 0:
                c = (self.pointer[0], self.pointer[1])
                self.path.append(c)
                self.pointer = (self.pointer[0], self.pointer[1]-1)
            # vertical only
            elif self.pointer[1] == 0:
                c = (self.pointer[0], self.pointer[1])
                self.path.append(c)
                self.pointer = (self.pointer[0]-1, self.pointer[1])
            else:
                #handling of general cases and more complex situations with branchings
                #can go diagonal
                if (self.m[self.pointer[0]][self.pointer[1]] - self.m[self.pointer[0]-1][self.pointer[1]-1] == self.match_score and self.match_matrix[self.pointer[0]][self.pointer[1]] == True) or (self.m[self.pointer[0]][self.pointer[1]] - self.m[self.pointer[0]-1][self.pointer[1]-1] == self.mismatch_score and self.match_matrix[self.pointer[0]][self.pointer[1]] == False):
                    current_path = [(i[0],i[1]) for i in self.path]
                    c = (self.pointer[0], self.pointer[1])
                    #Also vertical
                    if self.m[self.pointer[0]][self.pointer[1]] - self.m[self.pointer[0]-1][self.pointer[1]] == self.gap_score:
                        #register branch only if unknown
                        if (c,1,current_path) not in self.branches and (c,1,current_path) not in self.solved_branches:
                            self.branches.append((c,1,[(i[0],i[1]) for i in current_path]))
                        if (c,0,current_path) not in self.branches and (c,0,current_path) not in self.solved_branches:
                            self.branches.append((c,0,[(i[0],i[1]) for i in current_path]))
                            self.current_branch = (c,0,[(i[0],i[1]) for i in current_path])
                    #Also horizontal
                    if self.m[self.pointer[0]][self.pointer[1]] - self.m[self.pointer[0]][self.pointer[1]-1] == self.gap_score:
                        #register branch only if unknown
                        if (c,2,current_path) not in self.branches and (c,2,current_path) not in self.solved_branches:
                            self.branches.append((c,2,[(i[0],i[1]) for i in current_path]))
                        if (c,0,current_path) not in self.branches and (c,0,current_path) not in self.solved_branches:
                            self.branches.append((c,0,[(i[0],i[1]) for i in current_path]))
                            self.current_branch = (c,0,[(i[0],i[1]) for i in current_path])
                    #pointer block according to directive
                    self.path.append(c)
                    if self.directive != None:
                        if self.directive == 0:
                            self.current_branch = (c,0,[(i[0],i[1]) for i in current_path])
                            self.pointer = (self.pointer[0]-1, self.pointer[1]-1)
                            self.directive = None
                        elif self.directive == 1:
                            self.current_branch = (c,1,[(i[0],i[1]) for i in current_path])
                            self.pointer = (self.pointer[0]-1, self.pointer[1])
                            self.directive = None
                        elif self.directive == 2:
                            self.current_branch = (c,2,[(i[0],i[1]) for i in current_path])
                            self.pointer = (self.pointer[0], self.pointer[1]-1)
                            self.directive = None
                    else:
                        self.pointer = (self.pointer[0]-1, self.pointer[1]-1)
                        

                # NO diagonal, surely vertical
                elif self.m[self.pointer[0]][self.pointer[1]] - self.m[self.pointer[0]-1][self.pointer[1]] == self.gap_score:
                    current_path = [(i[0],i[1]) for i in self.path]
                    c = (self.pointer[0], self.pointer[1])
                    #Also horizontal
                    if self.m[self.pointer[0]][self.pointer[1]] - self.m[self.pointer[0]][self.pointer[1]-1] == self.gap_score:
                        #register branch only if unknown
                        if (c,2,current_path) not in self.branches and (c,2,current_path) not in self.solved_branches:
                            self.branches.append((c,2,[(i[0],i[1]) for i in current_path]))
                        if (c,1,current_path) not in self.branches and (c,1,current_path) not in self.solved_branches:
                            self.branches.append((c,1,[(i[0],i[1]) for i in current_path]))
                            self.current_branch = (c,1,[(i[0],i[1]) for i in current_path])
                    #pointer block according to directive
                    self.path.append(c)
                    if self.directive != None:
                        if self.directive == 0:
                            self.current_branch = (c,0,[(i[0],i[1]) for i in current_path])
                            self.pointer = (self.pointer[0]-1, self.pointer[1]-1)
                            self.directive = None
                        elif self.directive == 1:
                            self.current_branch = (c,1,[(i[0],i[1]) for i in current_path])
                            self.pointer = (self.pointer[0]-1, self.pointer[1])
                            self.directive = None
                        elif self.directive == 2:
                            self.current_branch = (c,2,[(i[0],i[1]) for i in current_path])
                            self.pointer = (self.pointer[0], self.pointer[1]-1)
                            self.directive = None
                    else:
                        self.pointer = (self.pointer[0]-1, self.pointer[1])
                
                # NO diagonal, surely horizontal
                elif self.m[self.pointer[0]][self.pointer[1]] - self.m[self.pointer[0]][self.pointer[1]-1] == self.gap_score:
                    current_path = [(i[0],i[1]) for i in self.path]
                    c = (self.pointer[0], self.pointer[1])
                #Also vertical
                    if self.m[self.pointer[0]][self.pointer[1]] - self.m[self.pointer[0]-1][self.pointer[1]] == self.gap_score:
                        #register branch only if unknown
                        if (c,1,current_path) not in self.branches and (c,1,current_path) not in self.solved_branches:
                            self.branches.append((c,1,[(i[0],i[1]) for i in current_path]))
                        if (c,2,current_path) not in branch_index and (c,2,current_path) not in solved_index:
                            self.branches.append((c,2,[(i[0],i[1]) for i in current_path]))
                            self.current_branch = (c,2,[(i[0],i[1]) for i in current_path])
                    #pointer block according to directive
                    self.path.append(c)
                    if self.directive != None:
                        if self.directive == 0:
                            self.current_branch = (c,0,[(i[0],i[1]) for i in current_path])
                            self.pointer = (self.pointer[0]-1, self.pointer[1]-1)
                            self.directive = None
                        elif self.directive == 1:
                            self.current_branch = (c,1,[(i[0],i[1]) for i in current_path])
                            self.pointer = (self.pointer[0]-1, self.pointer[1])
                            self.directive = None
                        elif self.directive == 2:
                            self.current_branch = (c,2,[(i[0],i[1]) for i in current_path])
                            self.pointer = (self.pointer[0], self.pointer[1]-1)
                            self.directive = None
                    else:
                        self.pointer = (self.pointer[0], self.pointer[1]-1)
                # Problems
                else:
                    raise ValueError('Backtracking error')


    def __align(self) -> None:
        '''Last step of the alignment, moslty graphical, but contains some checking of the alignment score, since during benchmarking
        some cases produced final paths with different scores from the expected one. ierates over the list of paths and appends to
        self.__alignment the command line representation of an alignment, if it has passed the final check.'''
        
        for i in self.paths:
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
                    matches.append('|')
                elif al1[i] != al2[i] and (al1[i] != '-' and al2[i] != '-'):
                    matches.append('*')
                else:
                    matches.append(' ')
            validate_score = matches.count('|')*self.match_score + matches.count('*')*self.mismatch_score + matches.count(' ')*self.gap_score
            if self.score == validate_score and ' '.join(al1) + '\n' + ' '.join(matches) + '\n' + ' '.join(al2) not in self.__alignment:
                if len(al1) <= 75:
                    self.__alignment.append([' '.join(al1) + '\n' + ' '.join(matches)+ f'  ({len(al1)})' + '\n' + ' '.join(al2)])
                else:
                    als = []
                    end = 0
                    for z in range(75,len(al1)+75,75):
                        if len(al1[z-75:z]) < 75:
                            al = ' '.join(al1[z-75:]) + '\n' + ' '.join(matches[z-75:])+ f'  ({end + len(al1[z-75:])})' + '\n' + ' '.join(al2[z-75:])
                            als.append(al)
                        else:
                            end += 75
                            al = ' '.join(al1[z-75:z]) + '\n' + ' '.join(matches[z-75:z])+ f'  ({end})' + '\n' + ' '.join(al2[z-75:z])
                            als.append(al)
                    self.__alignment.append(als)
    
    
    def getAlignment(self) -> None:
        '''Outputs nicely all the produced alignment; this was implemented to be the only legal way for the use to access the alignment'''
        if len(self.__alignment) > 0:
            print('='*150, end='\n\n')
            print('\t'+f'Global Score: {self.score}', '||', f'[{self.method.capitalize()} bactracking]', '||', f'Scoring: {self.match_score},{self.mismatch_score},{self.gap_score} [match,mismatch,gap]',end='\n\n')
            for i in range(len(self.__alignment)):
                print(f'[#{i+1}]', '-'*125, f'[Global score: {self.score}]', end='\n\n')
                for j in self.__alignment[i]:
                    print(j, end='\n\n')
            print('='*150, end='\n\n')
        else:
            if self.lenError:
                raise(ValueError('InputSequenceError: sequences to align must be at least 2 characters long'))
            else:
                print('='*50, end='\n\n')
                print('No global alignment found, probably due to a mistake in the processing or bad parsing', end = '\n\n')
                print('='*50, end='\n\n')
            

    def getMatrix(self) -> None:
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
