import numpy as np

def vecNum(vec):
    return int("".join(str(int(n)) for n in vec), base=2)

def printVec(vec, width):
    return "0b{:0{}b},".format(vecNum(vec), width)

def printBinary(mat, width):
    for vec in mat:
        print(printVec(vec, width))

def genToParityCheck(genParity):
    return np.hstack((genParity.transpose(), np.eye(genParity.shape[1])))

# This equals the extended matrix with the rightmost column removed.
parity = np.array([
    [1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0],
    [0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1],
    [1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0],
    [0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0],
    [0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1],
    [1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0],
    [0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0],
    [0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1],
    [1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1],
    [1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1],
    [1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1],
    [1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1],
])

print("core transpose:")
printBinary(parity.transpose(), 12)

parityCheck = genToParityCheck(parity)
print("parity check:")
printBinary(parityCheck, 23)
printBinary(parityCheck.transpose(), 11)

e = np.hstack((np.zeros(11), [1], np.zeros(11)))

syns = []
for i in range(12):
    pat = np.roll(e, -i)
    syns.append((pat @ parityCheck.transpose()) % 2)

print("syndromes:")
for syn in syns:
    print(printVec(syn, width=11))
