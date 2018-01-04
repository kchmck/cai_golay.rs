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

def genToAltParityCheck(genParity):
    return np.hstack((np.eye(genParity.shape[0]), genParity))

# This was pieced together from the P25 and DMR standards, as well as appendix Q of the
# IRIG standard.
parity = np.array([
    [1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1],
    [0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1],
    [1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0],
    [0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0],
    [0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0],
    [1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1],
    [0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1],
    [0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1],
    [1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0],
    [1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1],
    [1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0],
    [1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1],
])
gen = np.hstack((np.eye(12), parity))

print("core transpose:")
printBinary(parity.transpose(), 12)

print("core:")
printBinary(parity, 12)

parityCheck = genToParityCheck(parity)
print("parity check:")
printBinary(parityCheck, 24)

parityCheck = genToAltParityCheck(parity)
print("alt parity check:")
printBinary(parityCheck, 24)

# Verify self-dual property.
for r in range(0, 12):
    for q in range(r, 12):
        dot = (gen[r] @ gen[q]) % 2
        assert dot == 0.0
