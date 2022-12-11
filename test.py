"""
    This test.py file can be ignored it was only used
    for testing the architecture of this pBLAST algorithm.
"""

#import pBlast
#matrix = mainBlast.getMatrix()

kmerDictionary = {}
kmerDictionary['ISL'] = {'queryIndex': 5, 'databaseIndex': [7], 'initialHSSPValue': 16}
kmerDictionary['MKW'] = {'queryIndex': 0, 'databaseIndex': [1], 'initialHSSPValue': 21}

databaseDictionary = {}
databaseDictionary['Albumin'] = {'start': 0, 'end': 7}
databaseDictionary['Kinase'] = {'start': 7, 'end': 11}

query = 'MKWVTISL'
database = 'RMKWVTFISLLL'

databaseEndIndices = [7, 11]



def scoreBoardInitialization(columns, rows):
    scoreBoard = []
    length = min(columns, rows)

    for i in range(length + 1):
        scoreBoard.append([{'score': 0, 'direction': 'STOP'}])
        for j in range(length):
            scoreBoard[i].append({'score': 0, 'direction': 'STOP'})
    return scoreBoard


def getMatrixFragment(scoreBoard, coordinateArray):
    matrixFragment = [[0, 0], [0, 'x']]
    matrixFragment[0][0] = scoreBoard[coordinateArray[0]-1][coordinateArray[1]-1]
    matrixFragment[0][1] = scoreBoard[coordinateArray[0]-1][coordinateArray[1]]
    matrixFragment[1][0] = scoreBoard[coordinateArray[0]][coordinateArray[1]-1]
    return matrixFragment


def smithWaterman(matrixFragment, scoringMatrix, gapScore, aminoAcidArray):

    vertical = matrixFragment[0][1]['score'] + gapScore
    horizontal = matrixFragment[1][0]['score'] + gapScore
    diagonal = matrixFragment[0][0]['score'] + int(scoringMatrix[aminoAcidArray[0]][aminoAcidArray[1]])

    maxScore = max(vertical, horizontal, diagonal, 0)
    if maxScore == vertical:
        return {'score': vertical, 'direction': 'vertical', 'database': '-', 'query': aminoAcidArray[1]}
    elif maxScore == horizontal:
        return {'score': horizontal, 'direction': 'horizontal', 'database': aminoAcidArray[0], 'query': '-'}
    elif maxScore == diagonal:
        return {'score': diagonal, 'direction': 'diagonal', 'database': aminoAcidArray[0], 'query': aminoAcidArray[1]}
    else:
        return {'score': 0, 'direction': 'STOP'}



def searchAlgorithm(databaseEndIndices, kmerDictionary, databaseDictionary, scoringMatrix, database, query, extensionThreshold = 14, gapScore = -10, kmerLength = 3):

    resultDictionary = {}
    dictionaryIndex = 0

    for kmer in kmerDictionary:
        for dataIndex in kmerDictionary[kmer]['databaseIndex']:
            databaseStartIndex = dataIndex
            databaseEndIndex = 0

            for endIndex in databaseEndIndices:
                if endIndex == databaseStartIndex:
                    continue
                elif endIndex < databaseStartIndex:
                    continue
                else:
                    databaseEndIndex = endIndex
                    break

            if (databaseEndIndex - databaseStartIndex <= kmerLength):
                continue

            currentSequenceHeader = ''
            for header in databaseDictionary:
                if databaseDictionary[header]['end'] == databaseEndIndex:
                    currentSequenceHeader = header
                    break

            trueDatabase = database[dataIndex:databaseEndIndex]
            trueQuery = query[kmerDictionary[kmer]['queryIndex']:]
            queryLength = len(trueQuery)
            databaseLength = len(trueDatabase)
            scoreBoard = scoreBoardInitialization(databaseLength, queryLength)

            indexRange = queryLength + 1 if databaseLength >= queryLength else databaseLength + 1

            count = 0
            maxCoordinates = [0, 0]
            lastMax = 0
            for index in range(1, indexRange):
                maxScore = 0
                aScore, bScore, cScore = {'score': 0}, {'score': 0}, {'score': 0}
                coordinateArray = [i for i in range(1, index + 1)]
                lastElement = coordinateArray.pop()
                for element in coordinateArray:
                    a = getMatrixFragment(scoreBoard, [element, lastElement])
                    b = getMatrixFragment(scoreBoard, [lastElement, element])
                    aScore = smithWaterman(a, scoringMatrix, gapScore, [trueDatabase[element-1], trueQuery[lastElement-1]])
                    bScore = smithWaterman(b, scoringMatrix, gapScore, [trueDatabase[lastElement-1], trueQuery[element-1]])
                    scoreBoard[element][lastElement] = aScore
                    scoreBoard[lastElement][element] = bScore

                    if (index == indexRange-1):
                        if (aScore['score'] > lastMax):
                            lastMax = aScore['score']
                            maxCoordinates = [element, lastElement]
                        if (bScore['score'] > lastMax):
                            lastMax = bScore['score']
                            maxCoordinates = [lastElement, element]

                c = getMatrixFragment(scoreBoard, [lastElement, lastElement])
                cScore = smithWaterman(c, scoringMatrix, gapScore, [trueDatabase[lastElement-1], trueQuery[lastElement-1]])
                scoreBoard[lastElement][lastElement] = cScore

                if (index == indexRange-1):
                    if (cScore['score'] > lastMax):
                        lastMax = cScore['score']
                        maxCoordinates = [lastElement, lastElement]

                maxScore = max(aScore['score'], bScore['score'], cScore['score'])
                count = count + 1

                if (maxScore < extensionThreshold and count >= 3):
                    break



            alignedDatabase = ''
            alignedQuery = ''
            seperator = ''
            globalMaxScore = 0
            currentCoordinates = maxCoordinates

            for reverseIndex in range(0, indexRange-1):
                currentCell = scoreBoard[currentCoordinates[0]][currentCoordinates[1]]
                direction = currentCell['direction']
                globalMaxScore = currentCell['score'] + globalMaxScore

                alignedDatabase = currentCell['database'] + alignedDatabase
                alignedQuery = currentCell['query'] + alignedQuery

                if (currentCell['database'] == currentCell['query']):
                    seperator = '|' + seperator
                else:
                    seperator = ' ' + seperator

                if (direction == 'diagonal'):
                    currentCoordinates = [currentCoordinates[0]-1, currentCoordinates[1]-1]
                elif (direction == 'vertical'):
                    currentCoordinates = [currentCoordinates[0]-1, currentCoordinates[1]]
                elif (direction == 'horizontal'):
                    currentCoordinates = [currentCoordinates[0], currentCoordinates[1]-1]
                else:
                    break

            resultDictionary[dictionaryIndex] = {'Header': currentSequenceHeader, 'Score': globalMaxScore, 'AlignedDatabase': alignedDatabase, 'AlignedQuery': alignedQuery, 'Seperator': seperator}
            dictionaryIndex = dictionaryIndex + 1

    return resultDictionary




#resultDictionary = searchAlgorithm(databaseEndIndices, kmerDictionary, databaseDictionary, matrix, database, query)

#dict = dict(sorted(resultDictionary.items(), key=lambda item: item[1]['Score'], reverse=True))

#for result in dict:
#    object = resultDictionary[result]
#
#    print(object['Header'] + ' with score of:', object['Score'])
#    print('Database hit: ' + object['AlignedDatabase'])
#    print('              ' + object['Seperator'])
#    print('Query hit:    ' + object['AlignedQuery'], '\n')

