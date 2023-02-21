######## Convert interacting beads in txt file to csv. ########
####### ------>"May the Force serve u well..." <------#########
###############################################################

############# One above all #############
##-------------------------------------##
import pandas as pd
import sys
import os

file = sys.argv[1]

with open( file, "r" ) as f:
	txt = f.readlines()

prot1, prot2 = [], []
for i in range( len( txt )-1 ):
	if "\n" == txt[i]:
		continue
	elif txt[i+1].startswith( " [" ):
		
		p1 = txt[i].strip().split( ":" )
		line = txt[i+1].strip().split( "(" )
		for j in range( 1, len( line ) ):
			p2 = line[j].replace( "'", "").replace( ")", "").replace( "(", "").replace( "]", "").replace( "[", "").replace( " ", "").split( "," )
			prot1.append( p1[2] + "-" + p1[3] )
			prot2.append( p2[0] + "-" + p2[1] )

df = pd.DataFrame()
df["Protein1"] = prot1
df["Protein2"] = prot2

file_name = file.split( ".txt" )[0]
df.to_csv( f"{file_name}.csv", index = False )
