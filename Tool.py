import numpy as np
from MSLogging import logGetError

def toolCalCV(list_intensity: np.ndarray):

    list_intensity = list_intensity[list_intensity != 0]
    if len(list_intensity) < 2:
        return -1

    std = np.std(list_intensity)
    mean = np.mean(list_intensity)
    cv = std / mean

    return cv


def toolGetNameFromPath(path):

    lenStr = len(path)
    iStart = 0
    iEnd = -1

    for i in range(lenStr):

        j = lenStr - 1 - i

        if path[j] == '.':
            iEnd = j
            break

    for i in range(lenStr):

        j = lenStr - 1 - i

        if path[j] == '\\' or path[j] == '/':
            iStart = j + 1
            break

    return path[iStart:iEnd]


def toolGetWord(inputString, index, d):

    if inputString[0] != d:
        inputString = d + inputString

    if inputString[-1] != d:
        inputString = inputString + d

    p_d = []

    i = 0
    for c in inputString:

        if c == d:

            p_d.append(i)

        i = i+ 1

    result = inputString[p_d[index] + 1:p_d[index + 1]]

    return result


def toolGetWord1(inputString, d1, d2):

    start = 0
    end = len(inputString)

    for i in range(len(inputString)):

        if inputString[i] == d1:
            start = i + 1

        if inputString[i] == d2:
            end = i

    return inputString[start:end]



def toolCountCharInString(inputStr, inputChar):

    result = 0

    for c in inputStr:
        if c == inputChar:
            result = result + 1

    return result


def toolStr2List(inputStr, inputSeparator):

    outputList = []

    word = ''

    if inputStr[-1] != inputSeparator:
        inputStr = inputStr + inputSeparator

    for c in inputStr:

        if c == inputSeparator:

            number = float(word)
            outputList.append(number)
            word = ''

        else:

            word = word + c

    return outputList


def toolLinearRegression(X: np.ndarray, Y: np.ndarray):

    Array_len = X.shape[0]
    X_mean = X.sum() / Array_len
    Y_mean = Y.sum() / Array_len
    slope = np.dot(X, Y) / np.dot(X, X)
    pearson = np.dot(X - X_mean, Y - Y_mean) / np.sqrt(np.dot(X - X_mean, X - X_mean) * np.dot(Y - Y_mean, Y - Y_mean))

    return slope, pearson


def soldierBinarySearch(inputList, start, end, number):

    if end >= len(inputList):
        logGetError("MSTool,The length of inputList is " + str(len(inputList)) + ", but you want to get " + str(end) + "?")


    if number == inputList[end]:

        return end

    while start < end:
        mid = (start + end) // 2
        if number < inputList[mid]:
            end = mid
        elif number > inputList[mid]:
            start = mid + 1
        else:
            return mid

    start = start - 1

    return start


def toolFindNeighborFromSortedList0(inputList, start, end, number,):

    if number <= inputList[start]:
        return start

    if number >= inputList[end]:
        return end

    neighbor = soldierBinarySearch(inputList, start, end, number)

    disLeft = abs(inputList[neighbor] - number)
    disRight = abs(inputList[neighbor+1] - number)

    if disLeft < disRight:
        return neighbor
    else:
        return neighbor + 1


def toolFindNeighborFromSortedList1(inputList, number):

    return toolFindNeighborFromSortedList0(inputList, 0, len(inputList)-1, number)


def toolFindNeighborFromDisorderdList(inputList, number):

    minGap = abs(inputList[0] - number)
    minIndex = 0

    for i in range(1, len(inputList)):
        tmpGap = abs(inputList[i] - number)
        if tmpGap <= minGap:
            minGap = tmpGap
            minIndex = i

    return minIndex