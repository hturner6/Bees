def GetIntronRegions(dataframe):
    features = GetGeneExonStartAndEnd(dataframe)
    intronRegions = CreateIntronRegions(features)
    return intronRegions

#TODO: Implement
def CreateIntronRegions(features):
    introns = []
    currentGene = next((x for x in features if x.type == "gene"), None)
    previousBoundary = currentGene.start
    # How do we enforce that the 0th element is always a gene (it should be)
    for index in range(1, len(features)):
        currentFeature = features[index]
        if currentFeature.type == 'gene':
            # Get final intron from final exon end to gene end before starting on a new gene
            intronSize = previousBoundary - currentGene.end - 1
            intronExists = intronSize >= 1
            if (intronExists):
                newIntronRegion = (previousBoundary + 1, currentGene.end - 1, intronSize)
                introns.append(newIntronRegion)

            currentGene = currentFeature
            # Update boundary if current gene is further right (it should always be)
            # if (currentGene.end >= previousBoundary):
            previousBoundary = currentGene.start
            continue

        exonIsEntirelyWithinGene = currentFeature.start >= currentGene.start and currentFeature.end <= currentGene.end
        if (not exonIsEntirelyWithinGene):
            # Update boundary if current exon is further right
            if (currentFeature.end >= previousBoundary):
                previousBoundary = currentFeature.end
            continue

        # Make sure exons have gap
        # exonIsOverlapping = previousBoundary

        # Make sure there is a gap
        # We -1 as consecutive positions still do not have a gap
        intronSize = currentFeature.start - previousBoundary - 1
        intronExists = intronSize >= 1
        if (not intronExists):
            # Update boundary if current exon is further right
            if (currentFeature.end >= previousBoundary):
                previousBoundary = currentFeature.end
            continue

        newIntronRegion = (previousBoundary + 1, currentFeature.start - 1, intronSize)
        introns.append(newIntronRegion)
        
        # Update boundary if current exon is further right
        if (currentFeature.end >= previousBoundary):
            previousBoundary = currentFeature.end

    return introns
        
def GetGeneExonStartAndEnd(dataframe):
    featureNames = ['gene', 'exon']
    geneData = dataframe[dataframe['feature'].isin(featureNames)]

    # Order genes by start to ensure future elements positions always after previous
    geneData.sort_values('start')
    starts = geneData['start'].tolist()
    ends = geneData['end'].tolist()
    features = geneData['feature'].tolist()
    featureList = []
    for index in range(0, len(starts)):
       newFeature = Feature(features[index], starts[index], ends[index])
       featureList.append(newFeature)
    # return sorted(featureList, key=lambda x: x.start, reverse=False)
    return featureList


class Feature:
    def __init__(self, type, start, end):
        self.type = type
        self.start = start
        self.end = end