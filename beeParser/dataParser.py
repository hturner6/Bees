import pandas as pd

# Positional offset used to determine the up and downstream feature positions
# Is also used when determining the intergeneic regions
updownstreamOffset = 5000

def GetUpstreamAndDownStreamRegions(dataframe):
    # In format (chromeID, (starts, ends))
    geneData = GetGenesStartAndEnd(dataframe)
    updownstreamRegions = CreateUpdownstreamRegions(geneData)
    return updownstreamRegions

def GetIntergeneicRegions(dataframe):
    # In format (chromeID, (starts, ends))
    geneData = GetGenesStartAndEnd(dataframe)
    intergeneticRegions = CreateIntergeneicRegions(geneData)
    return intergeneticRegions

def GetIntronDataFromExternalFile(filePath):
    return pd.read_csv(filePath, sep='\t', header=None, names=['seqname', 'start', 'end', 'sign', 'gene_id1', 'gene_id2', 'length'])

def CreateUpdownstreamRegions(geneData):
    list = []
    for chromeIndex in range(0, len(geneData)):
        chromeID = geneData[chromeIndex][0]
        starts = geneData[chromeIndex][1]
        ends = geneData[chromeIndex][2]
        # Do for entire chrome
        for index in range(0, len(starts)):
            currentGeneStart = starts[index]
            currentGeneEnd = ends[index]

            (newDownstreamRegionStart, newDownstreamRegionEnd) = (currentGeneStart - updownstreamOffset, currentGeneStart - 1)
            list.append((chromeID, newDownstreamRegionStart, newDownstreamRegionEnd))

            (newUpstreamRegionStart, newUpstreamRegionEnd) = (currentGeneEnd + 1, currentGeneEnd + updownstreamOffset)
            list.append((chromeID, newUpstreamRegionStart, newUpstreamRegionEnd))
    return list

def CreateIntergeneicRegions(geneData):
    list = []
    for chromeIndex in range(0, len(geneData)):
        chromeID = geneData[chromeIndex][0]
        starts = geneData[chromeIndex][1]
        ends = geneData[chromeIndex][2]

        list.append((chromeID, 1, starts[0] - 1))
        for index in range(1, len(starts)):
            currentGeneEnd = ends[index - 1]
            nextGeneStart = starts[index]

            noIntergeneticRegion = GeneIsOverlapping(currentGeneEnd, nextGeneStart)
            if (noIntergeneticRegion):
                continue

            # start is current gene end + upstream value + 1
            (newIntergeneticRegionStart, newIntergeneticRegionEnd) = (currentGeneEnd + updownstreamOffset + 1, nextGeneStart - updownstreamOffset - 1)
            list.append((chromeID, newIntergeneticRegionStart, newIntergeneticRegionEnd))
    
    return list

# Returns a list in the format of [(chromeID, (starts, ends)), ...]
# Where the postions arrays are sorted by start positions
# (sorting is required by the other functions)
def GetGenesStartAndEnd(dataframe):
    geneData = dataframe[dataframe['feature'] == 'gene']

    chromes = geneData['seqname'].tolist()
    starts = geneData['start'].tolist()
    ends = geneData['end'].tolist()

    resultDict = AggregatePositionsByChromeId(chromes, starts, ends)
    resultData = []
    for chromeId, chromePositions in resultDict.items():
        sortedPositions = sorted(chromePositions, key=lambda x: x[0], reverse=False)
        unzipped = [[i for i, j in sortedPositions],
            [j for i, j in sortedPositions]]
        resultData.append((chromeId, unzipped[0], unzipped[1]))

    return resultData

def AggregatePositionsByChromeId(chromes, starts, ends):
    resultDict = {}
    for index in range(0, len(chromes)):
        chromeID = chromes[index]
        # If first time chrome is touched, create a new list
        if chromeID not in resultDict:
            resultDict[chromeID] = []
        
        resultDict[chromeID].append((starts[index], ends[index]))
    return resultDict

def GeneIsOverlapping(currentGeneEnd, nextGeneStart):
    return nextGeneStart <= currentGeneEnd