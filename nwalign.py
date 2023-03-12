
class NWAligner:
    def __init__(self, s1, s2, match_score = 1, mismatch_score = -1, gap_score = -2):
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
        self.alMat = []
        self.alMat.append([self.gap_score*i for i in range(self.cols)])
        for i in range(1, self.rows):
            row = []
            row.append(self.gap_score*i)
            for j in range(self.cols-1):
                row.append(None)
            self.alMat.append(row)
        self.score = self.fill()
        self.m = self.alMat[::][::]
        
        
        self.match_matrix = [[False for i in range(self.cols)]]
        for i in range(self.rows-1):
            row = []
            row.append(False)
            for j in range(self.cols-1):
                row.append(self.s1[i] == self.s2[j])
            self.match_matrix.append(row)
        
        self.backtrack()
        self.__recompose_paths()
        self.align()

    
    def fill(self):
        for i in range(1, self.rows):
            for j in range(1, self.cols):
                self.alMat[i][j] = self.__best_move((i,j))
        return self.alMat[i][j]
        
    def __best_move(self, coord):
        diag = self.alMat[coord[0]-1][coord[1]-1]
        if self.s1[coord[0]-1] == self.s2[coord[1]-1]:
            val = diag + self.match_score
        else:
            val = diag + self.mismatch_score
        vv = self.alMat[coord[0]-1][coord[1]] + self.gap_score
        hv = self.alMat[coord[0]][coord[1]-1] + self.gap_score
        return max(val, vv, hv)
    
    def backtrack(self):

        if (len(self.m), len(self.m[-1])) == (1,1):
            
            self.path.append((0,0))
        elif len(self.m) == 1:
            
            c = (len(self.m)-1, len(self.m[-1])-1)
            self.m = [[i for i in self.m[-1][:-1]]]
            self.match_matrix = [[i for i in self.match_matrix[-1][:-1]]]
            self.backtrack()
            self.path.append(c)
            
        elif len(self.m[-1]) == 1:
            
            c = (len(self.m)-1, len(self.m[-1])-1)
            self.m = self.m[:-1]
            self.match_matrix = self.match_matrix[:-1]
            self.backtrack()
            self.path.append(c)
            
        else:
            if (self.m[-1][-1] - self.m[-2][-2] == self.match_score and self.match_matrix[-1][-1] == True) or (self.m[-1][-1] - self.m[-2][-2] == self.mismatch_score and self.match_matrix[-1][-1] == False):
                
                current = self.m[::][::]
                current_match = self.match_matrix[::][::]
                if self.m[-1][-1] - self.m[-1][-2] != self.gap_score and self.m[-1][-1] - self.m[-2][-1] == self.gap_score:
                    
                    c = (len(self.m)-1, len(self.m[-1])-1)
                    self.m = self.m[:-1]
                    self.match_matrix = self.match_matrix[:-1]
                    self.backtrack()
                    self.path.append(c)
                    
                    self.paths.append(self.path[::])
                    self.path = []
                elif self.m[-1][-1] - self.m[-1][-2] == self.gap_score and self.m[-1][-1] - self.m[-2][-1] != self.gap_score:
                    c = (len(self.m)-1, len(self.m[-1])-1)
                    self.m = [[i for i in j[:-1]] for j in self.m]
                    self.match_matrix = [[i for i in j[:-1]] for j in self.match_matrix]
                    self.backtrack()
                    self.path.append(c)

                    self.paths.append(self.path[::])
                    self.path = []
                
                
                self.m = current
                self.match_matrix = current_match
                c = (len(self.m)-1, len(self.m[-1])-1)
                self.m = [[i for i in j[:-1]] for j in self.m[:-1]]
                self.match_matrix = [[i for i in j[:-1]] for j in self.match_matrix[:-1]]
                self.backtrack()
                self.path.append(c)
                
                
                
            elif self.m[-1][-1] - self.m[-1][-2] != self.gap_score and self.m[-1][-1] - self.m[-2][-1] == self.gap_score:
                
                c = (len(self.m)-1, len(self.m[-1])-1)
                self.m = self.m[:-1]
                self.match_matrix = self.match_matrix[:-1]
                self.backtrack()
                self.path.append(c)
                
            elif self.m[-1][-1] - self.m[-1][-2] == self.gap_score and self.m[-1][-1] - self.m[-2][-1] != self.gap_score:
                
                c = (len(self.m)-1, len(self.m[-1])-1)
                self.m = [[i for i in j[:-1]] for j in self.m]
                self.match_matrix = [[i for i in j[:-1]] for j in self.match_matrix]
                
                self.backtrack()
                self.path.append(c)
                
            else:
                raise ValueError('Backtracking error')
    
    def __recompose_paths(self):
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
        
    
    def align(self):
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
        print('='*50, end='\n\n')
        print('\t'+f'Score: {self.score}', end='\n\n')
        for i in self.alignment:
            print(i, end='\n\n')
        print('='*50, end='\n\n')
    
    @property
    def alignment(self):
        return self.__alignment


    def getMatrix(self):
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
