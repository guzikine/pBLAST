import pBlast

kmerLength = 3
thresholdValue = 13
gapScore = -10

if __name__ == "__main__":
    pBlast.runpBlast(kmerLength, thresholdValue, gapScore)