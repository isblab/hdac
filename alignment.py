import numpy as np
import pandas as pd

import Bio
from Bio import SeqIO
from Bio.Emboss.Applications import NeedleCommandline
from Bio import AlignIO



class Align_Paralogs():
	def __init__(self, hdac, exclude, outfile):
		self.hdac = hdac
		self.exclude = exclude    # If True, excludes the other paralog.
		self.remove = "HDAC1" if self.hdac == "HDAC2" else "HDAC2"
		self.data = pd.read_csv("../Production_runs/NuRD Complex Interactions.xlsx - NuRD Member XL Sites.csv")
		self.fasta_path = "../Uniprot/"
		self.output_path = "../EMBOSS/Needle/"
		self.outputfile_name = f"Mapped_XLs_{outfile}"   # Change the outputfile name accordingly.
		self.included_proteins = ["CHD", "GATA", "HDAC", "MBD", "MTA", "RBBP"]
		self.mapping = {i:[] for i in ["Protein1", "Residue1", "Protein2", "Residue2"]}


	def forward( self ):
		print("No. of XLs for each protein = ", self.data["Protein 1"].value_counts() )

		alignments_dict = self.map_aligned_pos()
		
		# exit()
		self.mapped_XLs( alignments_dict )
		df = pd.DataFrame()
		for key in self.mapping.keys():
			df[key] = self.mapping[key]

		gap_align = df[( df["Residue1"] == 0 ) | ( df["Residue2"] == 0 )].index
		df.drop( gap_align, inplace=True )
		print( "Dropped XLs containing residues aligned to a gap." )
		df.to_csv( f"{self.outputfile_name}.csv", index = False )
		print( "No. of XLs = ", len( df["Protein1"] ))
		# with open( f"{self.outputfile_name}.dat", "w" ) as w:
		# 	[w.writelines( "DSSO,{},{},{},{}\n".format(df.iloc[i,0], str(df.iloc[i,1]), df.iloc[i,2], str(df.iloc[i,3]) ) )  for i in range( len( df["Protein1"] ) )]


	def get_aligned_pos(self, target, query, position ):
		# Return the residue position from the aligned seq for protein2.
		# target --> sequence for which the aligned position is required (dominant paralog).
		# query --> sequence to which target is mapped (non-dominant paralog).
		# position --> residue position of query seq for which the aligned target residue position is required.

		iter1 = 1
		for i, pos in enumerate( query ):
			if pos != "-":
				if iter1 == position:
					idx = i
					break
				iter1 += 1
		
		iter2 = 1
		for j in range( len( target[:idx] ) ):
			# For residues mapped to a gap, 0 will be returned.
			if target[idx] == "-":
				return 0
				break
			elif target[j] != "-":
				iter2 += 1
			else:
				continue
		# print(idx, "  ", iter1, "  ", iter2, "  ", position)
		return iter2

	def Needle_alignment(self, seqa_path, seqb_path, outputfile ):
		# Sequence alignment using EMBOSS Needle.
		# seqa_path --> PATH for seqa fatsa file.
		# seqb_path --> PATH for seqb fatsa file.
		# outputfile --> Output file PATH.

		# Read the seq from fasra file.
		for record in SeqIO.parse( seqa_path, "fasta"):
			seqa = record.seq
		for record in SeqIO.parse( seqb_path, "fasta"):
			seqb = record.seq

		# Define parameters for the alignment.
		needle_cline = NeedleCommandline()
		needle_cline.asequence = "asis:" + seqa
		needle_cline.bsequence = "asis:" + seqb
		needle_cline.gapopen = 10
		needle_cline.gapextend = 0.5
		needle_cline.endopen = 10
		needle_cline.endextend = 0.5
		needle_cline.outfile = outputfile

		stdout, stderr = needle_cline()
		alignments = AlignIO.read( outputfile, "emboss" )

		# if "HDAC" in seqa_path:
		# 	print("".join([a for a in alignments[0]] )[89])
		# 	print("".join([a for a in alignments[1]] )[89])

		return "".join( [a for a in alignments[0]] ), "".join([b for b in alignments[1]] )


	def get_paralog(self, protein ):
		# Returns all the paralogs of a protein.
		if protein == "CHD":
			return ["CHD4", "CHD3"]
		if protein == "GATA":
			return ["GATAD2B", "GATAD2A"]
		if protein == "HDAC":
			return ["HDAC1", "HDAC2"]
		if protein == "MBD":
			return ["MBD3", "MBD2"]
		if protein == "MTA":
			return ["MTA1", "MTA2", "MTA3"]
		if protein == "RBBP":
			return ["RBBP4", "RBBP7"]
		else:
			print("Please enter the correct protein name...")
			exit()
		

	def return_dominant_paralog(self, prt ):
		# Returns the dominant paralog for a protein.
		if prt in ["CHD", "CHD3", "CHD4"]:
			return "CHD4"
		elif prt in ["GATA", "GATAD2A", "GATAD2B"]:
			return "GATAD2B"
		elif prt in ["HDAC1", "HDAC2"]:
			return "HDAC1"
		elif prt in ["MBD", "MBD2", "MBD3"]:
			return "MBD3"
		elif prt in ["MTA", "MTA1", "MTA2", "MTA3"]:
			return "MTA1"
		elif prt in ["RBBP", "RBBP4", "RBBP7"]:
			return "RBBP4"
		else:
			print(" --> ", prt)
			print("Please enter the correct protein name...")
			exit()


	def map_aligned_pos( self ):
		# Align the protein sequences and map the aligned positions to dominant isoform.
		alignments_dict = {}

		for protein in self.included_proteins:
			# Get all the paralogs for the protein.
			all_proteins = self.get_paralog( protein )
			
			# The first paralog has the most XLs (hardcoded as of now).
			p1 = all_proteins[0]

			alignments_dict[p1] = {}
			
			# Path for the protein fasta files.
			seqa_path = self.fasta_path + f"{protein}/{p1}_human.fasta"
			# For all the paralogs of a protein.
			for p2 in all_proteins[1:]:
				seqb_path = self.fasta_path + f"{protein}/{p2}_human.fasta"
				outputfile = self.output_path + f"{p1}-{p2}.txt"
				alignments_dict[p1][p2] = []
				seqa, seqb = self.Needle_alignment(seqa_path, seqb_path, outputfile)
				
				prot = []
				# Obtain all positions to be mapped for non-dominant paralog.
				[prot.append( pos ) for res in ["1", "2"] for i, pos in enumerate( self.data[f"Residue {res}"] ) if p2 in self.data[f"Protein {res}"][i]]
				# [original pos, aligned pos]
				[alignments_dict[p1][p2].append( [position, self.get_aligned_pos( seqa, seqb, position )] ) for position in prot]
				# if "HDAC" in protein:
				# 	print(alignments_dict["HDAC1"])

		return alignments_dict


	def mapped_XLs( self, alignments_dict ):
		for i in range( len( self.data["Protein 1"] ) ):
			if self.exclude and ( self.remove in self.data["Protein 1"][i] or self.remove in self.data["Protein 2"][i] ):
				continue
			else:
				prot1 = self.data["Protein 1"][i].replace( " ", "" )
				dom_prot1 = self.return_dominant_paralog( prot1 )
				
				if prot1 in dom_prot1:
					self.mapping["Protein1"].append( prot1 )
					self.mapping["Residue1"].append( self.data["Residue 1"][i] )
				# Consider either of HDAC1/HDAC2 XL and exclude the other.
				# If exclude = False, consider both.
				elif self.hdac in prot1 and self.exclude:
					self.mapping["Protein1"].append( dom_prot1 )
					self.mapping["Residue1"].append( self.data["Residue 1"][i] )			
				else:
					self.mapping["Protein1"].append( dom_prot1 )
					self.mapping["Residue1"].append( [pos[1] for pos in alignments_dict[dom_prot1][prot1] if pos[0] == self.data["Residue 1"][i]][0] )

				prot2 = self.data["Protein 2"][i].replace(" ", "")
				dom_prot2 = self.return_dominant_paralog( prot2 )

				if prot2 in dom_prot2:
					self.mapping["Protein2"].append( dom_prot2 )
					self.mapping["Residue2"].append( self.data["Residue 2"][i] )
				# Consider either of HDAC1/HDAC2 XL and exclude the other.
				# If exclude = False, consider both.
				elif self.hdac in prot2 and self.exclude:
					self.mapping["Protein2"].append( prot2 )
					self.mapping["Residue2"].append( self.data["Residue 2"][i] )
				else:
					self.mapping["Protein2"].append( dom_prot2 )
					self.mapping["Residue2"].append( [pos[1] for pos in alignments_dict[dom_prot2][prot2] if pos[0] == self.data["Residue 2"][i]][0] )


if __name__ == "__main__":
	print( "XL for HDAC1 only..." )
	Align_Paralogs("HDAC1", True, "HDAC1").forward()
	print( "\nXL for HDAC2 only..." )
	Align_Paralogs("HDAC2", True, "HDAC2").forward()
	print( "\nXL for HDAC1 and HDAC2 (with HDAC2 mapped to HDAC1)..." )
	Align_Paralogs("HDAC", False, "HDAC1-2").forward()


