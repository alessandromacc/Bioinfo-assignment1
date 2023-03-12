
class SWAligner:
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
        self.alignment = []
        self.alMat = []
        self.alMat.append([0 for i in range(self.cols)])
        for i in range(1, self.rows):
            row = []
            row.append(0)
            for j in range(self.cols-1):
                row.append(None)
            self.alMat.append(row)
        
        self.match_matrix = [[False for i in range(self.cols)]]
        for i in range(self.rows-1):
            row = []
            row.append(False)
            for j in range(self.cols-1):
                row.append(self.s1[i] == self.s2[j])
            self.match_matrix.append(row)

        self.fill()
        
        for i in self.starts:
            self.path = []
            self.am = [[k for k in j[:i[1]+1]] for j in self.alMat[:i[0]+1]]
            self.mm = [[k for k in j[:i[1]+1]] for j in self.match_matrix[:i[0]+1]]
            self.backtrack()
            self.paths.append(self.path[::-1])
        
        self.align()


        
    
    def fill(self):
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
        diag = self.alMat[coord[0]-1][coord[1]-1]
        if self.s1[coord[0]-1] == self.s2[coord[1]-1]:
            val = diag + self.match_score
        else:
            val = diag + self.mismatch_score
        vv = self.alMat[coord[0]-1][coord[1]] + self.gap_score
        hv = self.alMat[coord[0]][coord[1]-1] + self.gap_score
        return max(val, vv, hv, 0)
    
    def backtrack(self):
        
        if len(self.am) == 1:
            if self.am[-1][-2] == 0:
                self.path.append((len(self.am)-1, len(self.am[-1])-1))
                return
            else:
                self.path.append((len(self.am)-1, len(self.am[-1])-1))
                self.am = [[i for i in self.am[-1][:-1]]]
                self.mm = [[i for i in self.mm[-1][:-1]]]
                self.backtrack()
        elif len(self.am[-1]) == 1:
            if self.am[-2][-1] == 0:
                self.path.append((len(self.am)-1, len(self.am[-1])-1))
                return
            else:
                self.path.append((len(self.am)-1, len(self.am[-1])-1))
                self.am = self.am[:-1]
                self.mm = self.mm[:-1]
                self.backtrack()
        elif self.am[-1][-2] == self.am[-2][-1] == self.am[-2][-2] == 0:
            self.path.append((len(self.am)-1, len(self.am[-1])-1))
            return
        elif self.am[-1][-2] == self.am[-2][-1] == 0:
            if (self.am[-1][-1] - self.am[-2][-2] == self.match_score and self.mm[-1][-1] == True) or (self.am[-1][-1] - self.am[-2][-2] == self.mismatch_score and self.mm[-1][-1] == False):
                self.path.append((len(self.am)-1, len(self.am[-1])-1))
                self.am = [[k for k in j[:-1]] for j in self.am[:-1]]
                self.mm = [[k for k in j[:-1]] for j in self.mm[:-1]]
                self.backtrack()
            else:
                raise(ValueError('BacktrackingError: attempted diagonal movement'))
        elif self.am[-1][-2] == self.am[-2][-2] == 0:
            if self.am[-1][-1] - self.am[-2][-1] == self.gap_score:
                self.path.append((len(self.am)-1, len(self.am[-1])-1))
                self.am = self.am[:-1]
                self.mm = self.mm[:-1]
                self.backtrack()
            else:
                self.path.append((len(self.am)-1, len(self.am[-1])-1))
                return
            
        elif self.am[-2][-1] == self.am[-2][-2] == 0:
            if self.am[-1][-1] - self.am[-1][-2] == self.gap_score:
                self.path.append((len(self.am)-1, len(self.am[-1])-1))
                self.am = [[k for k in j[:-1]] for j in self.am]
                self.mm = [[k for k in j[:-1]] for j in self.mm]
                self.backtrack()
            else:
                self.path.append((len(self.am)-1, len(self.am[-1])-1))
                return

        elif self.am[-1][-2] == 0:
            if (self.am[-1][-1] - self.am[-2][-2] == self.match_score and self.mm[-1][-1] == True) or (self.am[-1][-1] - self.am[-2][-2] == self.mismatch_score and self.mm[-1][-1] == False) and (self.am[-2][-2] >= self.am[-2][-1]):
                self.path.append((len(self.am)-1, len(self.am[-1])-1))
                self.am = [[k for k in j[:-1]] for j in self.am[:-1]]
                self.mm = [[k for k in j[:-1]] for j in self.mm[:-1]]
                self.backtrack()
        elif self.am[-2][-1] == 0:
            if (self.am[-1][-1] - self.am[-2][-2] == self.match_score and self.mm[-1][-1] == True) or (self.am[-1][-1] - self.am[-2][-2] == self.mismatch_score and self.mm[-1][-1] == False) and (self.am[-2][-2] >= self.am[-1][-2]):
                self.path.append((len(self.am)-1, len(self.am[-1])-1))
                self.am = [[k for k in j[:-1]] for j in self.am[:-1]]
                self.mm = [[k for k in j[:-1]] for j in self.mm[:-1]]
                self.backtrack()
        elif self.am[-2][-2] == 0:
            if self.am[-1][-1] - self.am[-1][-2] != self.gap_score and self.am[-1][-1] - self.am[-2][-1] != self.gap_score:
                raise(ValueError('BacktrackingError: ambiguous NO-DIAGONAL situation - 1'))
            elif self.am[-1][-1] - self.am[-1][-2] == self.gap_score and self.am[-1][-1] - self.am[-2][-1] != self.gap_score:
                self.path.append((len(self.am)-1, len(self.am[-1])-1))
                self.am = [[k for k in j[:-1]] for j in self.am]
                self.mm = [[k for k in j[:-1]] for j in self.mm]
                self.backtrack()
            elif self.am[-1][-1] - self.am[-1][-2] != self.gap_score and self.am[-1][-1] - self.am[-2][-1] == self.gap_score:
                self.path.append((len(self.am)-1, len(self.am[-1])-1))
                self.am = self.am[:-1]
                self.mm = self.mm[:-1]
                self.backtrack()
            else:
                raise(ValueError('BacktrackingError: ambiguous NO-DIAGONAL situation - 2'))
        else:
            if (self.am[-1][-1] - self.am[-2][-2] == self.match_score and self.mm[-1][-1] == True) or (self.am[-1][-1] - self.am[-2][-2] == self.mismatch_score and self.mm[-1][-1] == False):
                self.path.append((len(self.am)-1, len(self.am[-1])-1))
                self.am = [[k for k in j[:-1]] for j in self.am[:-1]]
                self.mm = [[k for k in j[:-1]] for j in self.mm[:-1]]
                self.backtrack()
            else:
                raise(ValueError('BacktrackingError: attempted diagonal movement'))
    
    def align(self):
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
                    
            
            ls = len(max(self.s1[:i[0][0]], self.s2[:i[0][1]], key=len))
            rs = len(max(self.s1[:i[-1][0]], self.s2[:i[-1][1]], key=len))
            
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
        print('='*50, end='\n\n')
        print('\t'+f'Score: {self.score}', end='\n\n')
        for i in self.alignment:
            print(i, end='\n\n')
        print('='*50, end='\n\n')