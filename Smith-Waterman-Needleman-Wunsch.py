# Combine Smith Waterman and Needleman-Wunch algorithm as an interactive code using Python
# By Alireza Dantism
# Created on Thursday, June 30, 2022

# Function for reading files and get sequences
def read_sequence_file(filename):
    filereader = open(filename)
    line = filereader.readline()
    line = line.replace("\n", "")
    sequence = ""
    while True:
        line = filereader.readline()
        if not line:
            break
        sequence += line.replace("\n", "")
    return sequence

# 2 dimensional array for calculating Smith–Waterman scores
Matrix       = {}
Matrix[0, 0] = ""
Matrix[0, 1] = "_"

# 2 dimensional array for calculating trace back - each cell data comes from which one (Top or Top Left or Left)
TraceBackMatrix       = {}
TraceBackMatrix[0, 0] = ""
TraceBackMatrix[0, 1] = "_"
TraceBackMatrix[1, 0] = "_"

# at the end of algorithm - to store maximum value and it's coordinate we use a variable as below
FindMaxValueInMainMatrix        = 0
FindMaxValueInMainMatrixAddress = 0

# reading sequence from file and store in variable
sequence1 = read_sequence_file("sequence1.fasta")
sequence2 = read_sequence_file("sequence2.fasta")

print("\n ****** Hi Dear user, welcome to the program we need some score ******  \n")

print("Choose which algorithm do you want?")
ALGORITHM_TYPE   = input("Smith-Waterman or Needleman-Wunsch? | S/N ? ")
MATCH_VALUE      = int(input("Enter a Match score: "))
MISS_MATCH_VALUE = int(input("Enter a Mismatch score: "))
GAP_VALUE        = int(input("Enter a Gap score: "))

# check which sequence has a more amino acid
if len(sequence1) > len(sequence2):
  flag    = len(sequence1)
  Max_STR = len(sequence1)
  Min_STR = len(sequence2)
else:
  flag    = len(sequence2)
  Max_STR = len(sequence2)
  Min_STR = len(sequence1)

# initialize - first row and column of Matrix score and Matrix Coordinate
# in first row and column of value and coordinate matrix we show characters

for x in range(len(sequence1)):
  Matrix[0, x+2] = sequence1[x]
  TraceBackMatrix[0, x+2] = sequence1[x]

Matrix[0, 0] = ""
Matrix[1, 0] = "_"
for z in range(len(sequence2)):
  Matrix[z+2, 0] = sequence2[z]
  TraceBackMatrix[z+2, 0] = sequence2[z]

# initialize - second row and column of Matrix score and Matrix Coordinate
# in second row & column of value and coordinate matrix we show zero value -but for Needleman-Wunsch we should calculate

for x in range(len(sequence1) + 2):
    if ALGORITHM_TYPE == "S":
        # Smith-Waterman
        Matrix[1, x] = 0
        TraceBackMatrix[1, x] = 0
    else:
        # Needleman-Wunsch
        if x == 0 or x == 1:
            Matrix[1, x] = 0
            TraceBackMatrix[1, x] = 0
        else:
            Matrix[1, x] = Matrix[1, x-1] + GAP_VALUE
            TraceBackMatrix[1, x] = TraceBackMatrix[1, x-1] + GAP_VALUE

for x in range(len(sequence2) + 2):
    if ALGORITHM_TYPE == "S":
        # Smith-Waterman
        Matrix[x, 1] = 0
        TraceBackMatrix[x, 1] = 0
    else:
        # Needleman-Wunsch
        if x == 0 or x == 1:
            Matrix[x, 1] = 0
            TraceBackMatrix[x, 1] = 0
        else:
            Matrix[x, 1] = Matrix[x-1, 1] + GAP_VALUE
            TraceBackMatrix[x, 1] = TraceBackMatrix[x-1, 1] + GAP_VALUE


# here Smith-Waterman algorithm start and calculate each cell value
for i in range(len(sequence1)):
    for j in range(len(sequence2)):
        
        # MATCH
        if sequence1[i] == sequence2[j]:
            M = Matrix[j+1, i+1] + MATCH_VALUE
            m = -2000000
        # MIS MATCH
        else:
            m = Matrix[j+1, i+1] + MISS_MATCH_VALUE
            M = -2000000
               
        # Gap Penalty
        gTop  = Matrix[j+2, i+1] + GAP_VALUE
        gLeft = Matrix[j+1, i+2] + GAP_VALUE
        
        # Get max value from scores
        score_list = [M, m, gTop, gLeft]
        max_value = max(score_list)
        max_index = score_list.index(max_value)       
            
        # Check if max value refers to negative value then replace it with zero value
        if ALGORITHM_TYPE == "S":
            if max_value > 0:
                Matrix[j+2, i+2] = max_value
            else:
                Matrix[j+2, i+2] = 0
        else:
            Matrix[j+2, i+2] = max_value
            
        # add index for cell data - here we understand each cell data comes from (Top or Top Left or Left)
        
        # O  -> match
        # OO -> mismatch
        # T  -> gap from top seq
        # L  -> gap from left seq
        
        if max_index == 0:
            TraceBackMatrix[j+2, i+2] = str(j+1) + ',' + str(i+1) + ',' + "O"
        elif max_index == 1:
            TraceBackMatrix[j+2, i+2] = str(j+1) + ',' + str(i+1) + ',' + "OO"
        elif max_index == 2:
            TraceBackMatrix[j+2, i+2] = str(j+2) + ',' + str(i+1) + ',' + "T"
        else:
            TraceBackMatrix[j+2, i+2] = str(j+1) + ',' + str(i+2) + ',' + "L"

           
print("\n")
print("--------------- Here Is Smith-Waterman Matrix ---------------")
print("\n")

# we use temporary variable for check the length of sequences - here we need it for TraceBack Matrix
if len(sequence1) < len(sequence2):
    a = Max_STR
    b = Min_STR
else:
    a = Min_STR
    b = Max_STR

# here we plus 2 because first and second of rows and column is amino acid characters and zero value
for i in range(a+2):
    for j in range(b+2):
        print(Matrix[i, j], end=" \t")
        if i >= 2 and j >=2:
            if FindMaxValueInMainMatrix <= Matrix[i, j]:
                FindMaxValueInMainMatrix        = Matrix[i, j]
                FindMaxValueInMainMatrixAddress = str(i) + ',' + str(j)
    print("\n")
print("\n")

# Let's show backtrack matrix
print("--------------------- Here Is BackTrack Matrix ---------------------")
print("\n")


if len(sequence1) > len(sequence2):
    step_1 = Min_STR
    step_2 = Max_STR
    printAminos = step_2
else: 
    step_1 = Max_STR
    step_2 = Min_STR
    printAminos = step_1


for i in range(step_1+2):
    for j in range(step_2+2):
        print(TraceBackMatrix[i, j], end=" \t")
    print("\n")

print("** Max value in matrix is ", FindMaxValueInMainMatrix, " in coorinates of ", "Matrix[" + FindMaxValueInMainMatrixAddress + "] \n")

# here we show the alignments

firstSeq  = ""
secondSeq = ""

print("** Trace Back: ", Matrix[int(FindMaxValueInMainMatrixAddress.split(',')[0]), int(FindMaxValueInMainMatrixAddress.split(',')[1])],  end = " -> ")

flagR = 0
i     = 1
while i == 1:

    # get index of matrix value in back track matrix
    getPRVAddress = TraceBackMatrix[int(FindMaxValueInMainMatrixAddress.split(',')[0]), int(FindMaxValueInMainMatrixAddress.split(',')[1])]
    
    # get value of cell in matrix value
    # check if we are at the beginning of matrix
    if isinstance(getPRVAddress, int):
        getPRVValue = Matrix[0, 0]
    else:
        getPRVValue = Matrix[int(getPRVAddress.split(',')[0]), int(getPRVAddress.split(',')[1])]
    
    # get aminoacid character
    getPRVIndexTop  = TraceBackMatrix[int(FindMaxValueInMainMatrixAddress.split(',')[0]) - int(FindMaxValueInMainMatrixAddress.split(',')[0]), int(FindMaxValueInMainMatrixAddress.split(',')[1])]
    getPRVIndexLeft = TraceBackMatrix[int(FindMaxValueInMainMatrixAddress.split(',')[0]), int(FindMaxValueInMainMatrixAddress.split(',')[1]) - int(FindMaxValueInMainMatrixAddress.split(',')[1])]

    # each time of this loop we update this address to find previous amino acid
    FindMaxValueInMainMatrixAddress = getPRVAddress

    # O  -> match
    # OO -> mismatch
    # T  -> gap from top seq
    # L  -> gap from left seq

    # first of all check we are at the beginning or not
    if not isinstance(getPRVAddress, int):
        if getPRVAddress.split(',')[2] == "O":
            firstSeq  += getPRVIndexTop
            secondSeq += getPRVIndexLeft
        elif getPRVAddress.split(',')[2] == "OO":
            firstSeq  += getPRVIndexLeft
            secondSeq += getPRVIndexTop
        elif getPRVAddress.split(',')[2] == "L":
            firstSeq  += getPRVIndexLeft
            secondSeq += "-"
        elif getPRVAddress.split(',')[2] == "T":
            firstSeq  += "-"
            secondSeq += getPRVIndexTop
        else:
            firstSeq  += "-"
            secondSeq += "-"
    else:
        if isinstance(getPRVIndexLeft, str):
            firstSeq  += getPRVIndexLeft
            secondSeq += "-"
        else:
            firstSeq  += "-"
            secondSeq += getPRVIndexTop

    if ALGORITHM_TYPE == "N":
        flagR += 1
    
    if ALGORITHM_TYPE == "S":
        if getPRVValue == 0:
            i = 0        
            print(getPRVValue, end="")
        else:
            print(getPRVValue, end=" -> ")
    elif ALGORITHM_TYPE == "N":
        if flagR == printAminos:
            i = 0
            print(getPRVValue, end="")
        else:
            print(getPRVValue, end=" -> ")
        
        
print("\n")
  
# reverse characters and print
print(firstSeq[::-1])
print(secondSeq[::-1])

input()