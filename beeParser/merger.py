from gtfParser import gtf2df
from dataParser import GetIntergeneicRegions, GetUpstreamAndDownStreamRegions, GetIntronDataFromExternalFile
import pandas as pd

# Path to the existing .gtf file
gtfInputFilePath = ""
# Path to the output file (.csv). File must already be created
outputFilePath = ""
# To create the intron.bed file install gtftools and run:
# gtftools -i introns.bed [inputFileName].gtf
# Add the path to the created "introns.bed" file here
intronDataPath = ""

gtfDataframe = gtf2df(gtfInputFilePath)
updownstreamRegions = GetUpstreamAndDownStreamRegions(gtfDataframe)
intergeneicRegions = GetIntergeneicRegions(gtfDataframe)

updownstreamRegionsDataframe = pd.DataFrame(updownstreamRegions, columns=['seqname', 'start', 'end'])
updownstreamRegionsDataframe['feature'] = 'stream'
intergeneicRegionsDataframe = pd.DataFrame(intergeneicRegions, columns=['seqname', 'start', 'end'])
intergeneicRegionsDataframe['feature'] = 'intergenic'
intronDataFrame = GetIntronDataFromExternalFile(intronDataPath)
intronDataFrame['feature'] = 'intron'

wantedColumnsGTF = ['seqname', 'start', 'end', 'feature']

# Combines and adds the data to a new csv
# Assumes there is no overlapping data
allData = pd.concat([
    gtfDataframe.filter(wantedColumnsGTF),
    intergeneicRegionsDataframe.filter(wantedColumnsGTF),
    updownstreamRegionsDataframe.filter(wantedColumnsGTF),
    intronDataFrame.filter(wantedColumnsGTF)
], axis=0)

allData.sort_values(by=['seqname', 'start'], ascending=[True, True], inplace=True)
allData.to_csv(outputFilePath, sep='\t', index=False)