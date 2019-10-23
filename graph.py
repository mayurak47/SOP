import time

class AminoAcidNode:
    def __init__(self):
        self.residue = None
        self.secondary = None
        self.kappa = None
        self.alpha = None
        self.phi = None
        self.psi = None
        self.x_ca = None
        self.y_ca = None
        self.z_ca = None
        self.nh_bonded = None
        self.nh_bond_energy = None
        self.o_bonded = None
        self.o_bond_energy = None
        self.predicted_2d = None

class ProteinGraph:
    def __init__(self, filename):
        dsspfile = open(filename, "r")
        lines = dsspfile.readlines()
        self.chains = []
        curr_chain = []
        for i in range(28, len(lines)):
            if lines[i][13]=='!':
                self.chains.append(curr_chain)
                curr_chain = []
            else:
                curr_chain.append(AminoAcidNode())
                curr_chain[-1].residue = lines[i][13]
                if lines[i][16]==' ':
                    curr_chain[-1].secondary = 'C'
                else:
                    curr_chain[-1].secondary = lines[i][16]
                curr_chain[-1].nh_bonded = int(lines[i][41]+lines[i][42]+lines[i][43]+lines[i][44])
                curr_chain[-1].nh_bond_energy = float(lines[i][46]+lines[i][47]+lines[i][48]+lines[i][49])
                curr_chain[-1].o_bonded = int(lines[i][52]+lines[i][53]+lines[i][54]+lines[i][55])
                curr_chain[-1].o_bond_energy = float(lines[i][57]+lines[i][58]+lines[i][59]+lines[i][60])
        self.chains.append(curr_chain)
        self.adjacency_matrices = []
        for chain in self.chains:
            matrix = [[0 for _ in range(len(chain))] for _ in range(len(chain))]
            for i in range(0, len(chain)-1):
                # matrix[i][i+1] = 1
                # matrix[i+1][i] = 1
                if ((i+chain[i].nh_bonded)<len(chain)) and ((i+chain[i].nh_bonded)>=0):
                   matrix[i][i+chain[i].nh_bonded] = 1
                   matrix[i+chain[i].nh_bonded][i] = 1
                if ((i+chain[i].o_bonded)<len(chain)) and ((i+chain[i].o_bonded)>=0):
                   matrix[i][i+chain[i].o_bonded] = 1
                   matrix[i+chain[i].o_bonded][i] = 1
            self.adjacency_matrices.append(matrix)

    def predict2Dhelices(self):
        # for matrix in self.adjacency_matrices:
        for i in range(len(self.chains)):
            len_chain = len(self.chains[i])
            matrix = self.adjacency_matrices[i]
            curr_chain = self.chains[i]
            for j in range(len_chain):
                if j+4<len_chain and matrix[j][j+4]==1:
                    len_run = 0
                    temp_j = j
                    while temp_j+4<len_chain and matrix[temp_j][temp_j+4]==1:
                        len_run += 1
                        temp_j += 1
                    if len_run >= 4:
                        while j+4<len_chain and matrix[j][j+4]==1:
                            matrix[j][j+4] = 1
                            curr_chain[j].predicted_2d = 'H'
                            j += 1
                        curr_chain[j].predicted_2d = 'H'
                        curr_chain[j+1].predicted_2d = 'H'
                        curr_chain[j+2].predicted_2d = 'H'
                        curr_chain[j+3].predicted_2d = 'H'
                        j = j + 4
                elif j+3<len_chain and matrix[j][j+3]==1:
                    len_run = 0
                    temp_j = j
                    while temp_j+3<len_chain and matrix[temp_j][temp_j+3]==1:
                        len_run += 1
                        temp_j += 1
                    if len_run >= 3:
                        while j+3<len_chain and matrix[j][j+3]==1:
                            matrix[j][j+4] = 2
                            curr_chain[j].predicted_2d = 'G'
                            j += 1
                        curr_chain[j].predicted_2d = 'G'
                        curr_chain[j+1].predicted_2d = 'G'
                        curr_chain[j+2].predicted_2d = 'G'
                        j = j + 3
                elif j+5<len_chain and matrix[j][j+5]==1:
                    len_run = 0
                    temp_j = j
                    while temp_j+5<len_chain and matrix[temp_j][temp_j+5]==1:
                        len_run += 1
                        temp_j += 1
                    if len_run >= 5:
                        while j+5<len_chain and matrix[j][j+5]==1:
                            matrix[j][j+5] = 3
                            curr_chain[j].predicted_2d = 'I'
                            j += 1
                        curr_chain[j].predicted_2d = 'I'
                        curr_chain[j+1].predicted_2d = 'I'
                        curr_chain[j+2].predicted_2d = 'I'
                        curr_chain[j+3].predicted_2d = 'I'
                        curr_chain[j+4].predicted_2d = 'I'
                        j = j + 5

    def predict2Dbeta(self):
        for i in range(len(self.chains)):
            len_chain = len(self.chains[i])
            matrix = self.adjacency_matrices[i]
            curr_chain = self.chains[i]
            j = 0
            while j < len_chain:
                if curr_chain[j].predicted_2d==None or curr_chain[j].predicted_2d=='E':
                    if j+5<len_chain:
                        k=j+5
                        while k < len_chain and k < j+7 and (matrix[j][k]==0 or curr_chain[k].predicted_2d!=None):
                            k += 1
                        if not(k<len_chain and k<j+7 and matrix[j][k]==1 and curr_chain[k].predicted_2d==None):
                            j += 1
                            continue
                        print(str(j+1) + " " + str(k+1))
                        temp_j = j
                        temp_k = k
                        while temp_j >= 0 and temp_k < len_chain:
                            if matrix[temp_j][temp_k] == 1:
                                matrix[temp_j][temp_k] = 4
                                curr_chain[temp_j].predicted_2d = 'E'
                                curr_chain[temp_k].predicted_2d = 'E'
                                temp_j -= 1
                                temp_k += 1
                            elif matrix[temp_j][temp_k-1] == 1:
                                matrix[temp_j][temp_k] = 4
                                curr_chain[temp_j].predicted_2d = 'E'
                                curr_chain[temp_k].predicted_2d = 'E'
                                temp_j -= 1
                            elif matrix[temp_j+1][temp_k] == 1:
                                matrix[temp_j][temp_k] = 4
                                curr_chain[temp_j].predicted_2d = 'E'
                                curr_chain[temp_k].predicted_2d = 'E'
                                temp_k += 1
                            elif matrix[temp_j-1][temp_k+1] == 1:
                                matrix[temp_j][temp_k] = 4
                                curr_chain[temp_j].predicted_2d = 'E'
                                curr_chain[temp_k].predicted_2d = 'E'
                                temp_j -= 1
                                temp_k += 1
                            else:
                                j = k
                                break
                    else:
                        break
                else:
                    j+=1






ptn = ProteinGraph("/Users/mayurarvind/ptn.dssp")
# for node in ptn.chains[0]:
#     print(node.residue, end=" ")
# print()
# print()
# for node in ptn.chains[0]:
#     print(node.secondary, end=" ")
#
# print()
# print()

# for i in range(len(ptn.adjacency_matrices[0])):
#     for j in range(len(ptn.adjacency_matrices[0][0])):
#     #     print(ptn.adjacency_matrices[0][i][j], end=" ")
#     # print()
#         if ptn.adjacency_matrices[0][i][j] == 1:
#             print(str(i+1) + " " + str(j+1))
# print(len(ptn.chains[0]))
ptn.predict2Dhelices()
ptn.predict2Dbeta()
f = open("matrix.txt", "w")
for matrix in ptn.adjacency_matrices:
    for row in matrix:
        for col in row:
            f.write(str(col))
        f.write('\n')
    f.write('\n\n')

# f.write(ptn.adjacency_matrices)
f.close()
for i in range(len(ptn.chains[0])):
    print(str(i+1) + " " + str(ptn.chains[0][i].predicted_2d))
# for chain in ptn.chains:
#     for node in chain:
#         print(node.residue, end=" ")
#     print()
