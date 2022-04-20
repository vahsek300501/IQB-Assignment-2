import colorama
from termcolor import cprint
colorama.init()
print()
print()


inputSequence = "SGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDTVYCPRHVICTAEDMLNPNYEDLLIRKSNHSFLVQAGNVQLRVIGHSMQNCLLRLKVDTSNPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNHTIKGSFLNGSCGSVGF"
print("Input Sequence->   "+inputSequence)
print()
print()

strideOutput = "TTTT     HHHHHH EEEEEETTEEEEEEEETTEEEEEGGGG  HHHHH   HHHHHHH  GGG EEEETTEEE EEEEEEETTEEEEEE   TTTT        TTTEEEEEEEEETTEEEEEEEEEETTTT B    TTTTTTTEE "
codeOutput = "_HHHHHHHHHHHSSSSSSSSSSSSSSHHHHHHHHHSSSSHHHHHHHHHHHH__HHHHHHHHHHSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS____HSSSSSSSSSSSSSSSSSSSSS__SSSSSSSSSS___HHHHHH_________"

assert len(strideOutput) == len(codeOutput)

countSimilar = 0
foundMismatch = False
mismatchStartInd = -1
mismatchLengths = []

print("Code output->      ", end="")
for i in range(0, len(codeOutput)):
    if strideOutput[i] == codeOutput[i] or (strideOutput[i] == 'E' and codeOutput[i] == 'S'):
        cprint(codeOutput[i], "green", end="")
        countSimilar += 1
        if foundMismatch:
            foundMismatch = False
            mismatchLengths.append((mismatchStartInd, i-1))
    else:
        cprint(codeOutput[i], "red", end="")
        if foundMismatch == False:
            foundMismatch = True
            mismatchStartInd = i

if foundMismatch:
    mismatchLengths.append((mismatchStartInd, len(strideOutput)-1))
    foundMismatch = False

print()
print()

print("Stride output->    ", end="")
for i in range(0, len(codeOutput)):
    if strideOutput[i] == " ":
        cprint(" ", end="")
    elif strideOutput[i] == codeOutput[i] or (strideOutput[i] == 'E' and codeOutput[i] == 'S'):
        cprint(strideOutput[i], "green", end="")
    else:
        cprint(strideOutput[i], "red", end="")

print()
print()
print()

print("% Similarity between stride and code output: " +
      str((countSimilar/len(strideOutput))*100)+str("%"))

print()
print()

print("Dissimilarities indices and protien structure")
for val in mismatchLengths:
  print(str(val)+"   "+str(inputSequence[val[0]:val[1]+1]))
