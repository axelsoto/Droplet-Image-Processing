#-------------------------------------------------------------------------------
# Name:    dropImageProcessing. Version 1.0
# Purpose: Estimate changes in the droplet volume over time using different image frames. This is a proof of concept version. More advanced techniques are under development.

# Citation: Alejandro Martin Cohen, Axel Juan Soto, and James P Fawcett. "Determination of Flow Rates in Capillary Liquid Chromatography Coupled to a Nano Electrospray Source using Droplet Image Analysis Software" Analytical Chemistry
#
# Author:      Axel J. Soto
#
# Created:     20/05/2014
# Licence:     BSD
#-------------------------------------------------------------------------------

from PIL import Image, ImageFilter
import numpy
import time
import pdb
import sys

#These roughly indicate the boundary where the droplet should be. If set to the image size, it will take longer to compute
RELEVANTUPPERLIMIT = 68
RELEVANTLOWERLIMIT = 260
RELEVANTLEFTLIMIT = 0
RELEVANTRIGHTLIMIT = 475

#tolerance to identify difference in pixels intensity.
#toleranceContour = 230#230 for darker lighting (eg ./flow_400_c)
toleranceContour = 130#130 for regular lighting

#Constants from the image
pixelCapillary = 79
pixelDiameter = 360
umPerPixel = 4.55


def main(directory= './Sample videos/flow_300_c/',filePrefix = 'sample-',nFiles=47):
    start = time.clock()

    volumeArray,volumeDiffArray, volumeDiffArrayNeg, volumeArrayDropOnly,volumeDiffArrayDropOnly,volumeDiffArrayNegDropOnly, volumeArrayDropOnlyAsEllipsoid, volumeDiffArrayDropOnlyAsEllipsoid ,volumeDiffArrayNegDropOnlyAsEllipsoid = processFiles(directory, filePrefix, nFiles)
    end = time.clock()
    print "Time: " + str(end-start)
    saveVolumeDiffArrayOnFile(volumeDiffArray,volumeDiffArrayNeg,volumeArray,"dropAndCapillar")
    saveVolumeDiffArrayOnFile(volumeDiffArrayDropOnly,volumeDiffArrayNegDropOnly,volumeArrayDropOnly,"dropOnly")
    saveVolumeDiffArrayOnFile(volumeDiffArrayDropOnlyAsEllipsoid,volumeDiffArrayNegDropOnlyAsEllipsoid,volumeArrayDropOnlyAsEllipsoid,"dropOnlyAsEllipsoid")

    for interval in [2,4,8]:
        saveIntervalDifferences(volumeArray,interval,"volumeAtInterval"+str(interval)+".csv")
        saveIntervalDifferences(volumeArrayDropOnly,interval,"volumeAtInterval_dropOnly"+str(interval)+".csv")
        saveIntervalDifferences(volumeArrayDropOnlyAsEllipsoid,interval,"volumeAtInterval_dropOnlyAsEllipsoid"+str(interval)+".csv")

    #return array[0]


def processFiles(directory,filePrefix,nFiles):
    volumeArray = []
    volumeDiffArray = []
    volumeArrayDropOnly = []
    volumeDiffArrayDropOnly = []
    volumeArrayDropOnlyAsEllipsoid = []
    volumeDiffArrayDropOnlyAsEllipsoid = []
    for numb in range(1,nFiles+1):
        if numb !=1:
            #Skip first one because it is usually noisy

            #Calculate using drop and capillar
            vol, volAsEllipsoid = processBMP(directory + filePrefix + str(numb).zfill(3) + '.bmp',directory + 'out'+str(numb).zfill(3)+'.bmp', False)
            volumeArray.append(vol)
            if numb>2:
                currentFlow = volumeArray[numb-2]-volumeArray[numb-3]
                print "Current Flow: " + str(currentFlow)
                volumeDiffArray.append(currentFlow)

            #Calculate using drop only
            vol, volAsEllipsoid = processBMP(directory + filePrefix + str(numb).zfill(3) + '.bmp',directory + 'out'+str(numb).zfill(3)+'.bmp', True)
            volumeArrayDropOnly.append(vol)
            volumeArrayDropOnlyAsEllipsoid.append(volAsEllipsoid)
            if numb>2:
                currentFlow = volumeArrayDropOnly[numb-2]-volumeArrayDropOnly[numb-3]
                print "Current Flow (drop only): " + str(currentFlow)
                volumeDiffArrayDropOnly.append(currentFlow)

                currentFlow = volumeArrayDropOnlyAsEllipsoid[numb-2]-volumeArrayDropOnlyAsEllipsoid[numb-3]
                print "Current Flow (drop only as ellipsoid): " + str(currentFlow)
                volumeDiffArrayDropOnlyAsEllipsoid.append(currentFlow)


    volumeDiffArrayNeg = volumeDiffArray[:]
    volumeDiffArrayNegDropOnly = volumeDiffArrayDropOnly[:]
    volumeDiffArrayNegDropOnlyAsEllipsoid = volumeDiffArrayDropOnlyAsEllipsoid[:]

    volumeDiffArray = removeNegatives(volumeDiffArray)
    volumeDiffArrayDropOnly = removeNegatives(volumeDiffArrayDropOnly)
    volumeDiffArrayDropOnlyAsEllipsoid = removeNegatives(volumeDiffArrayDropOnlyAsEllipsoid)


    averageFlow = numpy.sum(numpy.array(volumeDiffArray))/len(volumeDiffArray)
    averageFlowDropOnly = numpy.sum(numpy.array(volumeDiffArrayDropOnly))/len(volumeDiffArrayDropOnly)
    averageFlowDropOnlyAsEllipsoid = numpy.sum(numpy.array(volumeDiffArrayDropOnlyAsEllipsoid))/len(volumeDiffArrayDropOnlyAsEllipsoid)

    print "Average Flow: " + str(averageFlow)
    print "Average Flow (drop only): " + str(averageFlowDropOnly)
    print "Average Flow (drop only as ellipsoid): " + str(averageFlowDropOnlyAsEllipsoid)
    return volumeArray,volumeDiffArray, volumeDiffArrayNeg, volumeArrayDropOnly,volumeDiffArrayDropOnly, volumeDiffArrayNegDropOnly, volumeArrayDropOnlyAsEllipsoid,volumeDiffArrayDropOnlyAsEllipsoid, volumeDiffArrayNegDropOnlyAsEllipsoid


def processBMP(filePath,filePathOut,onlyDrop):
    #onlyDrop True identifies the area of the drop only
    #onlyDrop false identifies drop and capillar

    original = Image.open(filePath)
    #rotated_img = original.rotate(5)
    converted_img = original.convert('L')
    ndarr = numpy.array(converted_img)
    ndarr3 = numpy.array(original)

    relevantZone = ndarr[RELEVANTUPPERLIMIT:RELEVANTLOWERLIMIT,RELEVANTLEFTLIMIT:RELEVANTRIGHTLIMIT]


    rightPoint = detectRightMostPoint(relevantZone,toleranceContour)

    if (onlyDrop):
        arrayOfTop, arrayOfBottom,beginOfDrop = detectTopAndBottomContourSmallDrop(relevantZone,rightPoint[1],rightPoint[0])
    else:
        arrayOfTop, arrayOfBottom,beginOfDrop = detectTopAndBottomContour(relevantZone,rightPoint[1],rightPoint[0])

    volumeAsEllipsoid = 0
    if (onlyDrop):
        volumeAsEllipsoid =  calculateVolumeAsEllipsoid(arrayOfTop,arrayOfBottom,umPerPixel)
    volume =  calculateVolume(arrayOfTop,arrayOfBottom,umPerPixel)

    ndarr3 = colorPoints(ndarr3,rightPoint[1] - beginOfDrop,arrayOfTop,arrayOfBottom,[237,24,72])

    converted_img = Image.fromarray(ndarr3)
    converted_img.save(filePathOut,original.format)

    return volume, volumeAsEllipsoid

def calculateVolume(arrayTop,arrayBottom,umPerPixel):
    differenceVector = numpy.diff([numpy.array(arrayBottom),numpy.array(arrayTop)],axis=0)
    differenceVectorRadious = differenceVector*0.5*umPerPixel

    volume = numpy.sum(differenceVectorRadious**2 * numpy.pi * umPerPixel)
    #change from um/s to fl/m
    volume = (volume/1000000) *60

    return volume

def calculateVolumeAsEllipsoid(arrayTop,arrayBottom,umPerPixel):
    differenceVector = numpy.diff([numpy.array(arrayBottom),numpy.array(arrayTop)],axis=0)

    radious1 = len(arrayTop) * 0.5 * umPerPixel
    radious2 = (max(arrayBottom) - min(arrayTop)) * 0.5 *umPerPixel
    radious3 = radious2

    volume = (4/3) * numpy.pi * radious1 * radious2 * radious3
    #change from um/s to fl/m
    volume = (volume/1000000) *60
    return volume

def colorPoints(matrix,xr,arrTop,arrBot,color):
    cont = 0
    for value in arrTop:
        matrix[RELEVANTUPPERLIMIT + value:RELEVANTUPPERLIMIT + arrBot[cont], RELEVANTLEFTLIMIT + xr]=[237,24,72]
        cont = cont + 1
        xr = xr - 1

    return matrix

def detectTopAndBottomContour(matrix,xr,yr):
    epsilon = 25
    arrayOutTop = []
    arrayOutBottom = []
    finishCondition = False

    yrTop = yr
    yrBottom = yr
    while not(finishCondition):
        yrTop = lowestWhite(matrix,xr,toleranceContour,max([yrTop-epsilon,0]))
        arrayOutTop.append(yrTop)

        yrBottom = highestWhite(matrix,xr,toleranceContour,min([yrBottom+epsilon,matrix.shape[0]-1]))
        arrayOutBottom.append(yrBottom)
        if xr == 0:
            finishCondition = True
        else:
            xr = xr - 1

    return arrayOutTop,arrayOutBottom, 0

def detectTopAndBottomContourSmallDrop(matrix,xr,yr):
    epsilon = 25
    context = 5 #context is the window size that is considered when evaluating different phases of the drop contour
    arrayOutTop = []
    arrayOutBottom = []
    finishCondition = False

    yrTop = yr
    yrBottom = yr

    topContourPhase = 0 #phase 1 goes up, phase 2 goes down, phase 3 goes up again. Phase 0 is the capillar before the drop is formed
    currentBPTop = 0
    bottomContourPhase = 0 #phase 1 goes down, phase 2 goes up, phase 3 goes down again
    currentBPBottom = 0
    topContourChangeIndex = []
    bottomContourChangeIndex = []
    while not(finishCondition):
        yrTop = lowestWhite(matrix,xr,toleranceContour,max([yrTop-epsilon,0]))
        arrayOutTop.append(yrTop)

        yrBottom = highestWhite(matrix,xr,toleranceContour,min([yrBottom+epsilon,matrix.shape[0]-1]))
        arrayOutBottom.append(yrBottom)

        if (currentBPTop==0):
            currentBPTop = yrTop
            currentBPBottom = yrBottom

        topContourPhase,currentBPTop = nextPhaseTop(arrayOutTop,topContourPhase,currentBPTop,bottomContourPhase,topContourChangeIndex)
        bottomContourPhase,currentBPBottom = nextPhaseBottom(arrayOutBottom,bottomContourPhase,currentBPBottom,topContourPhase,bottomContourChangeIndex)

        #if didn't find a proper drop
        #print(topContourPhase)
        if (xr<=1) and (topContourPhase==2):
            #print(arrayOutTop)
            topCopy = arrayOutTop[:]
            topCopy.reverse()
            topContourChangeIndex.append(len(topCopy)-topCopy.index(currentBPTop))
            bottomContourChangeIndex.append(len(topCopy)-topCopy.index(currentBPTop))
            topContourPhase = 3
            bottomContourPhase = 3

        #Assuming that we don't know what part is more notorious (would be nice to detect that when changing from phase 0 to 1)
        #if (topContourPhase==3) and (bottomContourPhase==3):
        #Assuming that the top part is always more notorious
        if (topContourPhase==3):
            if (bottomContourChangeIndex[-1]=='foo'):
                #Most likely bottomContour does not have a clear contour and hence didn't detect the change
                bottomCopy = arrayOutBottom[:]
                bottomCopy.reverse()
                bottomContourChangeIndex[-1] = (len(bottomCopy) - bottomCopy.index(currentBPBottom))

            finishCondition = True
            init = min(topContourChangeIndex[0], bottomContourChangeIndex[0])
            fin = max(topContourChangeIndex[-1], bottomContourChangeIndex[-1])
            arrayOutTop =  arrayOutTop[init:fin]
            arrayOutBottom = arrayOutBottom[init:fin]
        else:
            xr = xr - 1

    arrayOutTop,arrayOutBottom = modifyLeftEndingOfDrop(arrayOutTop, arrayOutBottom,topContourChangeIndex,bottomContourChangeIndex)
    return arrayOutTop, arrayOutBottom, init

def modifyLeftEndingOfDrop(arrayOutTop, arrayOutBottom,topContourChangeIndex,bottomContourChangeIndex):
    """ It tries to calculate how the ending is when this is distorted by the capillar"""
    diameter = max(arrayOutBottom) - min(arrayOutTop)
    vertCenter = min(arrayOutTop) + diameter/2
    topIsLater = (topContourChangeIndex[-1]>=bottomContourChangeIndex[-1])>0
    if (topIsLater):
        #The top part completes the drop after the bottom part
        #Left part
        for q in range((topContourChangeIndex[-1] - bottomContourChangeIndex[-1])):
            index = q + 1
            #This would mimic the top contour
            arrayOutBottom[-index] = max(vertCenter - arrayOutTop[-index],0) + vertCenter

        #Right part
        #Try to identify where the drop differentiates from the capillar
        beginningOfCircleFound = False
        approxValue = arrayOutBottom[0] + 2
        while not(beginningOfCircleFound):
            beginningOfCircle = find_element_in_list(approxValue, arrayOutBottom)
            if (beginningOfCircle>=0):
                beginningOfCircleFound = True
            approxValue = approxValue + 1
        #From the reconstructed part try to make the beginning similar to the top part but at the end, it should connect with the bottom part
        for q in range(beginningOfCircle):
            if (beginningOfCircle==1):
                arrayOutBottom[q] = ((max(vertCenter - arrayOutTop[q],0) + vertCenter) + arrayOutBottom[q])/2
            else:
                #arrayOutBottom[q] = ((max(vertCenter - arrayOutTop[q],0) + vertCenter) + arrayOutBottom[q])/2
                arrayOutBottom[q] = ((max(vertCenter - arrayOutTop[q],0) + vertCenter) * ((beginningOfCircle-1)-q)/(beginningOfCircle-1) + arrayOutBottom[q]*(q)/(beginningOfCircle-1))
    else:
        #The bottom part completes the drop after the top part
        #Left part
        for q in range((bottomContourChangeIndex[-1] - topContourChangeIndex[-1])):
            index = q + 1
            #This would mimic the top contour
            arrayOutTop[-index] = vertCenter - max(arrayOutBottom[-index] - vertCenter,0)

        #Right part
        #Try to identify where the drop differentiates from the capillar

        beginningOfCircleFound = False
        approxValue = arrayOutTop[0] - 2
        while not(beginningOfCircleFound):
            beginningOfCircle = find_element_in_list(approxValue, arrayOutTop)
            if (beginningOfCircle>=0):
                beginningOfCircleFound = True
            approxValue = approxValue - 1

        #From the reconstructed part try to make the beginning similar to the bottom part but at the end, it should connect with the top part

        for q in range(beginningOfCircle):
            if (beginningOfCircle==1):
                arrayOutTop[q] = (vertCenter - max(arrayOutBottom[q] - vertCenter,0) + arrayOutTop[q]) / 2
            else:
                #arrayOutTop[q] = (vertCenter - max(arrayOutBottom[q] - vertCenter,0) + arrayOutTop[q]) / 2
                arrayOutTop[q] = (vertCenter - max(arrayOutBottom[q] - vertCenter,0)) * ((beginningOfCircle-1)-q)/(beginningOfCircle-1) + arrayOutTop[q]*(q)/(beginningOfCircle-1)


    return arrayOutTop, arrayOutBottom


def find_element_in_list(element,list_element):
    try:
        index_element=list_element.index(element)
        return index_element
    except ValueError:
        return -1

def nextPhaseTop(arrayOutTop,topContourPhase,currentBreakPoint,bottomContourPhase,changeIndex):
    context = min(14,len(arrayOutTop))
    if (topContourPhase==2):
        if (currentBreakPoint < arrayOutTop[-1]):
            currentBreakPoint = arrayOutTop[-1]
            changeIndex[-1] = len(arrayOutTop)-1
        else:
            finishPhase = True
            cantSup = 0
            for q in range(1,context):
                if arrayOutTop[-q-1] < arrayOutTop[-q]:
                    finishPhase = False
                elif arrayOutTop[-q-1] > arrayOutTop[-q]:
                    cantSup = cantSup + 1
            if (finishPhase) and (cantSup>=2):
                topContourPhase = topContourPhase + 1

            #if (arrayOutTop[-2] >= arrayOutTop[-1]) and (arrayOutTop[-3] >= arrayOutTop[-2]) and (arrayOutTop[-4] >= arrayOutTop[-3]) and ((arrayOutTop[-4] > arrayOutTop[-1])) and (((len(arrayOutTop)-1)-changeIndex[-1])>=3):
            #    topContourPhase = topContourPhase + 1
    if (topContourPhase==1):
        if (currentBreakPoint > arrayOutTop[-1]):
            currentBreakPoint = arrayOutTop[-1]
            changeIndex[-1] = len(arrayOutTop)-1
        elif len(arrayOutTop)>8:
            finishPhase = True
            cantSup = 0
            for q in range(1,context):
                if arrayOutTop[-q-1] > arrayOutTop[-q]:
                    finishPhase = False
                elif arrayOutTop[-q-1] < arrayOutTop[-q]:
                    cantSup = cantSup + 1
            if (finishPhase) and (cantSup>=2):
                topContourPhase = topContourPhase + 1
                changeIndex.append('foo')

            #if (arrayOutTop[-2] <= arrayOutTop[-1]) and (arrayOutTop[-3] <= arrayOutTop[-2]) and (arrayOutTop[-4] <= arrayOutTop[-3]) and ((arrayOutTop[-4] < arrayOutTop[-1])) and (((len(arrayOutTop)-1)-changeIndex[-1])>=3):
            #    topContourPhase = topContourPhase + 1
            #    changeIndex.append('foo')
    if (topContourPhase==0):
        if (bottomContourPhase>0) or ((len(arrayOutTop)>=2) and ((arrayOutTop[-2] - arrayOutTop[-1])>3)):
            topContourPhase = topContourPhase + 1
            changeIndex.append(len(arrayOutTop)-1)
            changeIndex.append('foo')

    return topContourPhase,currentBreakPoint

def nextPhaseBottom(arrayOutBottom,bottomContourPhase,currentBreakPoint,topContourPhase,changeIndex):
    context = min(14,len(arrayOutBottom))
    if (bottomContourPhase==2):
        if (currentBreakPoint > arrayOutBottom[-1]):
            currentBreakPoint = arrayOutBottom[-1]
            changeIndex[-1] = len(arrayOutBottom)-1
        else:
            finishPhase = True
            cantSup = 0
            for q in range(1,context):
                if arrayOutBottom[-q-1] > arrayOutBottom[-q]:
                    finishPhase = False
                elif arrayOutBottom[-q-1] < arrayOutBottom[-q]:
                    cantSup = cantSup + 1
            if (finishPhase) and (cantSup>=1):
                bottomContourPhase = bottomContourPhase + 1
            #if (arrayOutBottom[-2] <= arrayOutBottom[-1]) and (arrayOutBottom[-3] <= arrayOutBottom[-2]) and (arrayOutBottom[-4] <= arrayOutBottom[-3]) and ((arrayOutBottom[-4] < arrayOutBottom[-1])) and (((len(arrayOutBottom)-1)-changeIndex[-1])>=3):
            #    bottomContourPhase = bottomContourPhase + 1
    if (bottomContourPhase==1):
        if (currentBreakPoint < arrayOutBottom[-1]):
            currentBreakPoint = arrayOutBottom[-1]
            changeIndex[-1] = len(arrayOutBottom)-1
        elif len(arrayOutBottom)>8:
            finishPhase = True
            cantSup = 0
            for q in range(1,context):
                if arrayOutBottom[-q-1] < arrayOutBottom[-q]:
                    finishPhase = False
                elif arrayOutBottom[-q-1] > arrayOutBottom[-q]:
                    cantSup = cantSup + 1
            if (finishPhase) and (cantSup>=1):
                bottomContourPhase = bottomContourPhase + 1
                changeIndex.append('foo')
            #if (arrayOutBottom[-2] >= arrayOutBottom[-1]) and (arrayOutBottom[-3] >= arrayOutBottom[-2]) and (arrayOutBottom[-4] >= arrayOutBottom[-3]) and ((arrayOutBottom[-4] > arrayOutBottom[-1])) and (((len(arrayOutBottom)-1)-changeIndex[-1])>=3):
            #    bottomContourPhase = bottomContourPhase + 1
            #    changeIndex.append('foo')
    if (bottomContourPhase==0):
        if (topContourPhase>0) or ((len(arrayOutBottom)>=2) and ((arrayOutBottom[-2] - arrayOutBottom[-1])>3)) :
            bottomContourPhase = bottomContourPhase + 1
            changeIndex.append(len(arrayOutBottom)-1)
            changeIndex.append('foo')

    return bottomContourPhase, currentBreakPoint


def highestWhite(matrix,x,tol,startY):
    for rowIndex in range(startY,-1,-1):

        if (matrix[rowIndex,x] < (255-tol)):
            axisRow = rowIndex
            break

    return axisRow

def lowestWhite(matrix,x,tol,startY):
    for rowIndex in range(startY,matrix.shape[0]):
        if (matrix[rowIndex,x] < (255-tol)):
            axisRow = rowIndex
            break

    return axisRow

def detectRightMostPoint(npArray,tol):
    #find the upper border of the upper part

    for colIndex in range(npArray.shape[1]-1,-1,-1):
        if (npArray[:,colIndex].min() < (250-tol)):
            axisRow = npArray[:,colIndex].argmin(axis=0)
            axisCol = colIndex

            break

    return axisRow, axisCol

def removeNegatives(values):
    newArray=[]
    for v in values:
        if v>0:
            newArray.append(v)
    return newArray

def saveIntervalDifferences(volume,interval,fileName):
    volumeDiff = []
    f = open(fileName,'w')
    for q in range(0, len(volume)-interval):
        volumeDiff.append((volume[q+interval] - volume[q])/interval)

    for elem in volumeDiff:
        f.write(str(elem) + '\n')
    f.close()

def saveVolumeDiffArrayOnFile(arr,arr2,arr3,prefix):
    f = open(prefix + '_DiffVolume.csv','w')
    for elem in arr:
        f.write(str(elem) + '\n')
    f.close()

    f = open(prefix + '_DiffNegVolume.csv','w')
    for elem in arr2:
        f.write(str(elem) + '\n')
    f.close()

    f = open(prefix + '_Volume.csv','w')
    for elem in arr3:
        f.write(str(elem) + '\n')
    f.close()

if __name__ == '__main__':
    if len(sys.argv)==1:
        main()
    elif len(sys.argv)==5:
        toleranceContour = int(sys.argv[4])
        main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
    else:
        print "Provide as arguments [input directory], [file prefix], [# files] and [tolerance contour]. No arguments for default behaviour"
