#!/usr/bin/env python

################################################################################
## This script takes outputs from NUMT_detection.sh to generates NUMT clusters
################################################################################

import fileinput
import sys, os
import pandas as pd
import glob
import scipy.stats as stats
import numpy as np

#################### Extract cluster from mtDNA discordant sam files  #############################

#clustering function
def cluster(data, maxgap):
	data.sort() #sort by name
	groups = [[data[0]]] #groups is nested list of first column of dataframe
	for x in data[1:]: #for each value 
		if abs(x - groups[-1][-1]) <= maxgap:
            		groups[-1].append(x)
		else:
			groups.append([x])
	return groups  #returns nested list with positions that are all within at least 500 bp of each other

#pass arguments from command line
sampleID,wgsBAM,input1 = sys.argv[1:]

#rename sample ID to remove last part
sampleID = sampleID.replace(".mt.disc","")

#sets discordant read file to dataframe
df0 = pd.read_csv(input1, names =['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL','SM','RG','NM','BC','OC','ZX','ZY','SA','NA1','NA2'], on_bad_lines='skip', sep="\t", engine='python', comment="@")

#filtering rows that contain RNAME "Un_|random|\."
df = df0[~df0.RNAME.str.contains("Un_|random|\.")]
#make new dataframe with MAPQ (mapping quality) > 0, RNAME = "MT" (mitochondrial) OR RNAME = chrM (mitochondrial), or the mate pair/ next read in pair is mitochondrial
df1 = df[(df["MAPQ"].astype(int) > 0) & ((df['RNAME'] == "MT") | (df['RNAME'] == "chrM") | (df['RNEXT'] == "MT") | (df['RNEXT'] == "chrM"))]

##### remove reads map to mtDNA #####
df3 = df1[~df1.RNAME.isin(['chrM','MT'])]

##### order by chromosome and pos #####
df3.sort_values(['RNAME','POS'], inplace=True)

##### look for the clutser by mapgap on pos #####
##### extract the clusters with number of elements no less than 5 #####

#et empty dataframe
output1 = pd.DataFrame([])

#create grouped dataframe by name of chromosome
df_chr = df3.groupby(['RNAME'])

for clusterID, myclusters in df_chr:
    clusterID = clusterID[0]
    myclusters['POS'] = myclusters['POS'].astype(int)
    sub_cluster = cluster(myclusters['POS'].tolist(), maxgap=500)
    output = pd.DataFrame([])
    for x in sub_cluster: #goes through and adds clusters that are >5 read pairs to df
        if len(x)>=5:   #make sure list is list of lists is >5
            mycluster = x  #name cluster
            df_cluster = df3[df3.POS.isin(mycluster)] #make df that has not only positions but also all other info
            df_cluster_pairMT = df[df.QNAME.isin(df_cluster['QNAME'])]  #go back to original sorted dataframe and get pairs that match QNAME
            mt_cluster = cluster(df_cluster['PNEXT'].tolist(), maxgap=500)  #make nested list of clusters of next read in mate pair
            for y in mt_cluster:  
                df_cluster_pairMT_out = df_cluster_pairMT[df_cluster_pairMT.PNEXT.isin(y)]
                df_cluster_pairMT_out['subCluster_No'] = len(y)
                df_cluster_pairMT_out['Cluster_No'] = len(x)			
                df_cluster_pairMT_out['Cluster_ID'] = clusterID + "_" + str(min(df_cluster['POS'])) + "_" + \
                                                      str(max(df_cluster['POS'])) + "_MTboth_"  + \
                                                      str(min(df_cluster_pairMT_out['PNEXT'])) + "_" + \
                                                      str(max(df_cluster_pairMT_out['PNEXT']))
                #output = output.append(df_cluster_pairMT_out, ignore_index=True)
                output = pd.concat([output, df_cluster_pairMT_out], ignore_index=True)
            #output1 = output1.append(output)
            output1 = pd.concat([output1, output], ignore_index=True)
	        #output2 = df[df.QNAME.isin(output1['QNAME'])]
            output1["IndividualID"] = sampleID
            if len(output1) != 0:
                output1['clusterID'] = output1['RNAME'].astype(str) + "_" + output1['Cluster_No'].astype(str)
                output2 = output1.groupby(['IndividualID','Cluster_ID','Cluster_No','subCluster_No']).size().to_frame('size').reset_index()
            
            else: 
                next
#output2 = pd.DataFrame(output2)
#output2.columns=['IndividualID','Cluster_ID','Cluster_No','subCluster_No','size']
output1.to_csv(input1 + '.cluster.tsv', sep = '\t',header=True, index=False)
output2.to_csv(input1 + '.cluster.summary.tsv', sep = '\t',header=False, index=True)

##### generate cluster range for defining the breakpoints ######

#output2['Cluster_ID'].replace("_MT","",inplace=True)
#output2 = output2.drop_duplicates(["Cluster_ID"])
output2['disFile'] = input1
output2['splitFile'] = input1.replace('disc','split')
output2['wgsBAM'] = wgsBAM
df_pos = pd.DataFrame(output2['Cluster_ID'].str.split('_').tolist(),columns = ['chr','start','end','chrM','mtstart','mtend'])
del output2['Cluster_ID']
del output2['subCluster_No']
del output2['size']
output3 = pd.concat([output2, df_pos[['chr','start','end']]], axis=1)
output3 = output3.drop_duplicates(['chr','start','end'])
output3['start'] = output3['start'].astype(int) - 500
output3['end'] = output3['end'].astype(int) + 500 + 150
output3['Cluster_No'] = output3['Cluster_No'].astype(int)

##### output to file #####
output3.to_csv(input1 + '.breakpointINPUT.tsv', sep = '\t',header=False, index=False)
