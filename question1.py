import colorama
from termcolor import cprint
colorama.init()

# Initializing the propensity matrix of alpha helix and beta sheets structure
alphaHelixPropensity = {'E': 1.53, 'A': 1.45, 'L': 1.34, 'H': 1.24, 'M': 1.20, 'Q': 1.17, 'W': 1.14, 'V': 1.14, 'F': 1.12,
                        'K': 1.07, 'I': 1.00, 'D': 0.98, 'T': 0.82, 'S': 0.79, 'R': 0.79, 'C': 0.77, 'N': 0.73, 'Y': 0.61, 'P': 0.59, 'G': 0.53}

betaStrandPropensity = {'M': 1.67, 'V': 1.65, 'I': 1.60, 'C': 1.30, 'Y': 1.29, 'F': 1.28, 'Q': 1.23, 'L': 1.22, 'T': 1.20,
                        'W': 1.19, 'A': 0.97, 'R': 0.90, 'G': 0.81, 'D': 0.80, 'K': 0.74, 'S': 0.72, 'H': 0.71, 'N': 0.65, 'P': 0.62, 'E': 0.26}

# Initializing input sequence
inputSequence = "SGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDTVYCPRHVICTAEDMLNPNYEDLLIRKSNHSFLVQAGNVQLRVIGHSMQNCLLRLKVDTSNPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNHTIKGSFLNGSCGSVGF"


# Function to get the right extension when alpha helix is detected
def getRightExtensionAlpha(endInd):
    global inputSequence, alphaHelixPropensity
    propensityVal = 0
    finalEndInd = endInd
    while(True):
        propensityVal = 0
        if endInd + 1 >= len(inputSequence):
            break
        propensityVal += (alphaHelixPropensity[inputSequence[endInd]] + alphaHelixPropensity[inputSequence[endInd-1]] +
                          alphaHelixPropensity[inputSequence[endInd-2]] + alphaHelixPropensity[inputSequence[endInd+1]])
        if propensityVal < 4:
            break
        finalEndInd = endInd+1
        endInd = endInd+1
    return finalEndInd


# Function to get left extenision when alpha helix is detected
def getLeftExtensionAlpha(startInd):
    global inputSequence, alphaHelixPropensity
    propensityVal = 0
    finalStartInd = startInd
    while(True):
        propensityVal = 0
        if startInd - 1 < 0:
            break
        propensityVal += (alphaHelixPropensity[inputSequence[startInd]] + alphaHelixPropensity[inputSequence[startInd+1]] +
                          alphaHelixPropensity[inputSequence[startInd+2]] + alphaHelixPropensity[inputSequence[startInd-1]])
        if propensityVal < 4:
            break
        finalStartInd = startInd-1
        startInd = startInd-1
    return finalStartInd


# Function to get neucleation sites for alpha helix
# Six consecutive residues are considered and if
# atleast 4 residues have p(H) > 1 it is considered as
# neucleation site
def getAlphaHelixRegions():
    global alphaHelixPropensity, inputSequence
    alphaHelixRanges = []

    for i in range(0, len(inputSequence)):
        totalPropensity = 0
        if i + 5 < len(inputSequence):
            for j in range(0, 6):
                if alphaHelixPropensity[inputSequence[i+j]] >= 1:
                    totalPropensity += 1
            if totalPropensity >= 4:
                alphaHelixRanges.append(
                    (getLeftExtensionAlpha(i), getRightExtensionAlpha(i+5)))
    return alphaHelixRanges


# Function to get the right extension when beta sheets is detected
def getRightExtensionBeta(endInd):
    global inputSequence, betaStrandPropensity
    propensityVal = 0
    finalEndInd = endInd
    while(True):
        propensityVal = 0
        if endInd + 1 >= len(inputSequence):
            break
        propensityVal += (betaStrandPropensity[inputSequence[endInd]] + betaStrandPropensity[inputSequence[endInd-1]] +
                          betaStrandPropensity[inputSequence[endInd-2]] + betaStrandPropensity[inputSequence[endInd+1]])
        if propensityVal < 4:
            break
        finalEndInd = endInd+1
        endInd = endInd+1
    return finalEndInd


# Function to get the left extension when beta sheets is detected
def getLeftExtensionBeta(startInd):
    global inputSequence, betaStrandPropensity
    propensityVal = 0
    finalStartInd = startInd
    while(True):
        propensityVal = 0
        if startInd - 1 < 0:
            break
        propensityVal += (betaStrandPropensity[inputSequence[startInd]] + betaStrandPropensity[inputSequence[startInd+1]] +
                          betaStrandPropensity[inputSequence[startInd+2]] + betaStrandPropensity[inputSequence[startInd-1]])
        if propensityVal < 4:
            break
        finalStartInd = startInd-1
        startInd = startInd-1
    return finalStartInd


# Function to get neucleation sites for beta sheets
# Five consecutive residues are considered and if
# atleast 3 residues have p(H) > 1 it is considered as
# neucleation site
def getBetaHelixRegions():
    global betaStrandPropensity, inputSequence
    betaHelixRanges = []

    for i in range(0, len(inputSequence)):
        totalPropensity = 0
        if i+4 < len(inputSequence):
            for j in range(0, 5):
                if betaStrandPropensity[inputSequence[i+j]] >= 1:
                    totalPropensity += 1
            if totalPropensity >= 3:
                betaHelixRanges.append(
                    (getLeftExtensionBeta(i), getRightExtensionBeta(i+j)))
    return betaHelixRanges


# Function to get propensity score for alpha helix and beta strand for a give range
def getAlphaBetaScoreForRange(startInd, endInd):
    global alphaHelixPropensity, betaStrandPropensity, inputSequence
    alphaScore = 0
    betaScore = 0

    for i in range(startInd, endInd+1):
        alphaScore += alphaHelixPropensity[inputSequence[i]]
        betaScore += betaStrandPropensity[inputSequence[i]]
    return alphaScore, betaScore


# Function to resolve overlaps when both the structures are detected
def resolveOverlaps(overlappingIntervals, finalSeq):
    for interval in overlappingIntervals:
        alphaScore, betaScore = getAlphaBetaScoreForRange(
            interval[0], interval[1])
        finalStruct = None
        if alphaScore > betaScore:
            finalStruct = 'H'
        else:
            finalStruct = 'S'
        for i in range(interval[0], interval[1]+1):
            finalSeq[i] = finalStruct


# Function to return final sequence
def generateFinalSequence():
    global inputSequence
    alphaHelix = getAlphaHelixRegions()
    betaStrand = getBetaHelixRegions()
    tmpAlphaLst = []
    tmpBetaLst = []
    finalSeq = []

    for _ in range(0, len(inputSequence)):
        tmpAlphaLst.append('_')
        tmpBetaLst.append('_')
    for i in range(0, len(alphaHelix)):
        for j in range(alphaHelix[i][0], alphaHelix[i][1]+1):
            tmpAlphaLst[j] = 'H'

    for i in range(0, len(betaStrand)):
        for j in range(betaStrand[i][0], betaStrand[i][1]+1):
            tmpBetaLst[j] = 'S'

    overlappings = []
    overlappingFound = False
    overlapStartInd = -1

    for i in range(0, len(inputSequence)):
        if tmpAlphaLst[i] == 'H' and tmpBetaLst[i] != 'S':
            finalSeq.append('H')
            if overlappingFound:
                overlappings.append((overlapStartInd, i-1))
                overlappingFound = False
                overlapStartInd = -1
        elif tmpAlphaLst[i] != 'H' and tmpBetaLst[i] == 'S':
            finalSeq.append('S')
            if overlappingFound:
                overlappings.append((overlapStartInd, i-1))
                overlappingFound = False
                overlapStartInd = -1

        elif tmpAlphaLst[i] != 'H' and tmpBetaLst[i] != 'S':
            finalSeq.append('_')
            if overlappingFound:
                overlappings.append((overlapStartInd, i-1))
                overlappingFound = False
                overlapStartInd = -1

        else:
            finalSeq.append('H/S')
            if(overlappingFound == True):
                continue
            overlappingFound = True
            overlapStartInd = i

    if(overlappingFound):
        overlappings.append((overlapStartInd, len(inputSequence)-1))
        overlapStartInd = -1
        overlappingFound = False

    resolveOverlaps(overlappings, finalSeq)
    return finalSeq


# Main function used for printing and displaying
def main():
    global inputSequence

    print()
    print()

    print("Input Sequence->     "+str(inputSequence))

    print()
    print()

    print("Predicted Output->   ", end="")
    finalSeq = generateFinalSequence()
    for val in finalSeq:
        if val == 'H':
            cprint(val, 'cyan', end="")
        elif val == 'S':
            cprint(val, 'yellow', end="")
        else:
            cprint(val, 'white', end="")
    print()


main()
