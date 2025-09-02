#!/usr/bin/env python
import numpy as np 
import pandas as pd 
import sys, os
# read in breakpoints data which is output of NUMT pipeleine
Input_Breakpoints_File=sys.argv[1]

# make sure file has content
if (os.stat(Input_Breakpoints_File).st_size == 0):
    print('File is empty!')
    print('Exiting program...')
    sys.exit(0)
NUMT_output=pd.read_csv(Input_Breakpoints_File, sep = '\t', header= None)
# get only Mito breakpoints
filtered = NUMT_output[NUMT_output[0].isin(['mt_Tstart','mt_Tend'])]
# reformat 
reformatted_columns = ['Type','Chromosome','Pos','Strand','Sample']
rf = pd.DataFrame(columns = reformatted_columns)
for row in range(len(filtered)):
    if filtered.iloc[row][0] == 'mt_Tstart':
        vals = filtered.iloc[row]
        Type = vals[0]
        Chromosome = vals[2]
        Pos = vals[6]
        Strand = vals[4]
        Sample = vals[7]
        rf.loc[row] = [Type, Chromosome, Pos, Strand, Sample]
    elif filtered.iloc[row][0] == 'mt_Tend':
        vals = filtered.iloc[row]
        Type = vals[0]
        Chromosome = vals[2]
        Pos = vals[3]
        Strand = vals[4]
        Sample = vals[7]
        rf.loc[row] = [Type, Chromosome, Pos, Strand, Sample]
    else:
        print("Missing Start or End Value")

# group values by start or end, position, and strand
grouped = filtered.groupby([0,6,4]).first()
grouped.reset_index(inplace=True)

# first sort by strand then type[start|end]in case there are multiuple NUMTs in different regions that are not overlapping
grouped = grouped.sort_values(by=[4,0], ascending=[True,False])

# group by strand and type to make sure there are not multiple start and end sites for each position
name_counts=grouped[0].value_counts()

#check if keys exist
key1='mt_Tstart'
key2='mt_Tend'
if key1 in name_counts and key2 in name_counts:
    pass
else:
    print("Error: No NUMT Startpoints identified")
    print("Exiting program...")
    sys.exit(0)  


# when both keys exist make sure there are not multiple
if name_counts['mt_Tstart'] > 1:
    print("WARNING: MULTIPLE NUMTS Startpoints IDENTIFIED, please curate manually")
elif name_counts['mt_Tend'] > 1:
    print("WARNING: MULTIPLE NUMTS Endpoints IDENTIFIED, please curate manually")
elif name_counts['mt_Tstart'] < 1:
    print("Error: No NUMT Startpoints identified")
    print("Exiting program...")
    sys.exit(0)
elif name_counts['mt_Tend'] < 1:
    print("Error: No NUMT Endpoints identified")
    print("Exiting program...")
    sys.exit(0)  
# else:
#     print("Error: No NUMT Endpoints identified")
#     print("Exiting program...")
#     sys.exit(0)  

# write to bed file
bed_input = ['Chromosome','Start','End','Name','Score','Strand']
bed = pd.DataFrame(columns = bed_input)
Chromosome = 'chrM'
Start = grouped.iloc[0,1]
End = grouped.iloc[1,5]
Name = grouped.iloc[0,7]
Score = 0
Strand = grouped.iloc[0,2]
bed.loc[1] = [Chromosome, Start, End, Name, Score, Strand]

nm = str.split(os.path.basename(Input_Breakpoints_File),'.')[0]

#insert check intofilename
for val in name_counts.values:
    if val > 1:
        ext = '_NUMT_regions_ERROR_raised_compared_to_input_file.bed'
    else:
        ext = '_NUMT_regions.bed'
filename = f"{nm}_{Chromosome}_{Start}_{End}{ext}"
bed.to_csv(filename, sep='\t', index=False)
print("SUCCESS!!!")
