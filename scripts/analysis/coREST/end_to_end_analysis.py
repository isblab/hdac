################ Perform end-to-end analysis for IMP models ################
############## ------>"May the Force serve u well..." <------###############
############################################################################

############# One above all #############
##-------------------------------------##
import numpy as np
import sys
import os
import subprocess

############ This is the way ############
##-------------------------------------##
# Check all the path variables before proceeding.
# Following scripts must be in the same directory --
	# run_analysis_trajectories.py
	# run_extract_models.py
	# variable_filter_v1.py
	# get_xl_viol_validation_set.py
	# contact_maps_all_pairs_surface.py
# Define the protein pairs and protein lengths for contact map script.
class IMPAnalysis():
	def __init__( self ):
		self.imp_path = "/home/kartik/imp-clean/build/setup_environment.sh"
		self.modeling_dir_path = "../modeling/"
		self.output_dir = "Output_"
		self.density_file = "./density_coREST_hdac1.txt"
		self.var_filter_path = "./Variable_filter"
		self.sampcon_path = "/home/kartik/imp-clean/imp/modules/sampcon/pyext/src/exhaust.py"
		self.sampcon_cluster = 0
		self.xl_viol_path = "./XL_Violations"
		self.XLs_path = "../../Data/XL/coREST_Mapped_XLs_HDAC1.csv"
		self.XL_cutoff = 35
		self.cmap_path = "./Contact_Maps"
		self.prism_path = "./prism"
		self.prism_script = "/home/kartik/Documents/PrISM/prism/src/main.py"
		self.prism_colour_script = "/home/kartik/Documents/PrISM/prism/src/color_precision.py"


	def forward( self ):
		self.cluster = self.return_major_cluster()
		self.pmi_analysis()
		self.model_extraction()
		self.sampcon()
		self.fit_to_inputXLs()
		self.contact_maps()
		self.prism()

		print( "Alas, the journey comes to an end...\n" )
		print( "May the Force be with you..." )


	def return_major_cluster( self, cluster_num = None ):
		# cluster_num --> if None returns models in the major cluster 
		# 					else returns models in the speciffied cluster.
		file = "./model_analysis/summary_hdbscan_clustering.dat"

		with open (file, "r") as f:
			cluster_summary = f.readlines()

		models_count = {}
		#Get the index for the N_models header
		index = cluster_summary[0].split(",").index("N_models") 

		# Get cluster no. and model counts for all clusters.
		for i in range(1,len(cluster_summary)):
			line = cluster_summary[i].split(",")
			if line[0] != "-1":
				models_count[int(line[0])] = int(line[index])
			else:
				continue                              # "Jaisa chal raha h, chalne do!!!"

		if cluster_num == None:
			clust = [(i,models_count[i]) for i in range(len(models_count)) if models_count[i] == max( models_count.values() )][0]
			print("Max models in cluster: ", clust[0], " no. of models = ", clust[1])
		else:
			clust = ( cluster_num, models_count[cluster_num] )
			print("Cluster: ", clust[0], " no. of models = ", clust[1])
		return clust


	def pmi_analysis( self ):
		############# PMI Analysis ##############
		##-------------------------------------##
		print("\n<-----------PMI Analysis----------->")
		subprocess.call( [f"{self.imp_path}", 
							"python", 
							"run_analysis_trajectories.py", 
							f"{self.modeling_dir_path}", 
							f"{self.output_dir}"] )
		self.cluster = return_major_cluster( cluster_num = None )



	def model_extraction( self ):
		############ Model Extraction ###########
		##-------------------------------------##
		print("\n<-----------Model Extraction----------->")
		if int( cluster[1] ) >= 30000:	
			print(">30000 models....\n using Variable filter \n")
			if os.path.exists( self.var_filter_path ):
				os.rename( f"{self.var_filter_path}/variable_filter_v1.py", 
								"./variable_filter_v1.py" )
			else:
				os.mkdir( self.var_filter_path )

			subprocess.call( ["python", "variable_filter_v1.py", "-c", f"{self.cluster[0]}"] )
			subprocess.call( [f"{self.imp_path}", 
								"python", 
								"run_extract_models.py", 
								f"{self.modeling_dir_path}", 
								f"{self.output_dir}", 
								f"{self.cluster[0]}", "True"] )

			os.rename( "ignore.KS_Test.txt", f"{self.var_filter_path}/ignore.KS_Test.txt" )
			os.rename( "ignore.Score_Hist_A.txt", f"{self.var_filter_path}/ignore.Score_Hist_A.txt" )
			os.rename( "ignore.Score_Hist_B.txt", f"{self.var_filter_path}/ignore.Score_Hist_B.txt" )
			os.rename( "var_filt_out.log", f"{self.var_filter_path}/var_filt_out.log" )
			os.rename( "var_filt_out.png", f"{self.var_filter_path}/var_filt_out.png" )
			os.rename( "variable_filter_v1.py", f"{self.var_filter_path}/variable_filter_v1.py" )
		else:
			subprocess.call( [f"{self.imp_path}", 
								"python", 
								"run_extract_models.py", 
								f"{self.modeling_dir_path}", 
								f"{self.output_dir}", 
								f"{self.cluster[0]}", "False"] )



	def sampcon( self ):
		########### Sampcon analysis ############
		##-------------------------------------##
		print("\n<-----------Sampcon----------->")
		subprocess.call([
			f"{self.imp_path}",
			"python",
			f"{self.sampcon_path}",
			"-n", "coREST",
			"-m", "cpu_omp",
			"-c", "4",
			"-d", f"{self.density_file}",
			"-gp",
			"-g", "5",
			"-sa", f"model_analysis/A_models_clust{self.cluster[0]}.txt",
			"-sb", f"model_analysis/B_models_clust{self.cluster[0]}.txt",
			"-ra", f"model_analysis/A_models_clust{self.cluster[0]}.rmf3",
			"-rb", f"model_analysis/B_models_clust{self.cluster[0]}.rmf3",
			"--align",
			"--prism",
			"--skip",
			"--cluster_threshold", "15"
			])


	def fit_to_inputXLs( self ):
		########### Fit to Input XLs ############
		##-------------------------------------##
		print("\n<-----------Fit to input XLs----------->")
		# Create XL_violations directory if not already existing.
		if not os.path.exists( self.xl_viol_path ):
			os.mkdir( self.xl_viol_path )
		# Move the xl_violation script inside the XL_violations dir if not already there.
		if not os.path.exists( f"{self.xl_viol_path}/get_xl_viol_validation_set.py" ):
			os.rename( "get_xl_viol_validation_set.py", f"{self.xl_viol_path}/get_xl_viol_validation_set.py" )
		os.chdir( self.xl_viol_path )

		# os.system(
		# 	f"{imp_path} python get_xl_viol_validation_set.py \
		# 	-ia ../cluster.0.sample_A.txt \
		# 	-ib ../cluster.0.sample_B.txt \
		# 	-ra ../model_analysis/A_models_clust{cluster[0]}.rmf3 \
		# 	-rb ../model_analysis/B_models_clust{cluster[0]}.rmf3 \
		# 	-c ../cluster.0/cluster_center_model.rmf3 \
		# 	-ta ../model_analysis/A_models_clust{cluster[0]}.txt \
		# 	-x {XLs_path} \
		# 	-t {XL_cutoff}"
		# 	)

		subprocess.call([
			f"{self.imp_path}",
			"python",
			"get_xl_viol_validation_set.py",
			"-ia", "../cluster.0.sample_A.txt",
			"-ib", "../cluster.0.sample_B.txt",
			"-ra", f"../model_analysis/A_models_clust{self.cluster[0]}.rmf3",
			"-rb", f"../model_analysis/B_models_clust{self.cluster[0]}.rmf3",
			"-c", f"../cluster.{self.sampcon_cluster}/cluster_center_model.rmf3",
			"-ta", f"../model_analysis/A_models_clust{self.cluster[0]}.txt",
			"-x", f"{self.XLs_path}",
			"-t", f"{self.XL_cutoff}"
			])
		os.chdir( "../" )

		# exit()


	def contact_maps( self ):
		############# Contact Maps ##############
		##-------------------------------------##
		print("\n<-----------Creating Contact maps----------->")
		# Create the contact map directory if not already existing.
		if not os.path.exists( self.cmap_path ):
			os.mkdir( self.cmap_path )
		# Move the contact map script inside the contact map dir if not already there.
		if not os.path.exists( f"{self.cmap_path}/contact_maps_all_pairs_surface.py" ):
			os.rename( "./contact_maps_all_pairs_surface.py", f"{self.cmap_path}/contact_maps_all_pairs_surface.py" )
		os.chdir( self.cmap_path )

		protein1 = ["HDAC1.0", "RCOR1.0"] 
		protein2 = ["RCOR1.0", "KDM1A.0"] 

		prots = []

		for index1 in range(len(protein1)):
			for index2 in range(len(protein2)):
				prot1 = protein1[index1]
				prot2 = protein2[index2]
				if prot1 == prot2:
					continue
				elif f"{prot2},{prot1}" in prots:
					continue
				else:
					prots.append(f"{prot1},{prot2}")

		prots.append( "RCOR1.0,RCOR1.0" )

		print( f"Contact maps will be created for the following protein pairs -- \n{prots}" )
		print( f"No. of Contact maps created = {len( prots )}" )
		for prot in prots:
			os.system(
				f"~/imp-clean/build/setup_environment.sh \
				python contact_maps_all_pairs_surface.py \
				-ia ../cluster.{self.sampcon_cluster}.sample_A.txt \
				-ib ../cluster.{self.sampcon_cluster}.sample_B.txt \
				-ra ../model_analysis/A_models_clust{self.cluster[0]}.rmf3 \
				-rb ../model_analysis/B_models_clust{self.cluster[0]}.rmf3 \
				-c ../cluster.{self.sampcon_cluster}/cluster_center_model.rmf3 \
				-ta ../model_analysis/A_models_clust{self.cluster[0]}.txt \
				-p {prot} &" 
				)
		os.chdir( "../" )


	def prism( self ):
		################ PrISM ##################
		##-------------------------------------##
		print("\n<-----------Precision Analysis: PrISM----------->")
		# Create teh Prism directory if not already existing.
		if not os.path.exists( self.prism_path ):
			os.mkdir( self.prism_path )
		os.chdir( self.prism_path )

		subprocess.call( [
			f"{self.imp_path}",
			"python", f"{self.prism_script}",
			"--input", f"../cluster.{self.sampcon_cluster}.prism.npz", 
			"--input_type", "npz",
			"--output", f"output{self.sampcon_cluster}",
			"--voxel_size", "2",
			"--return_spread",
			"--classes", "2",
			"--cores", "50",
			"--models", "1.0",
			"--n_breaks", "50",
			"--resolution", "1"
			] )

		print( "Colouring the PrISM output..." )
		subprocess.call( 
			[f"{self.imp_path}",
			"python", f"{self.prism_colour_script}",
			"--resolution", "1",
			"--annotations_file", f"./output{self.sampcon_cluster}/annotations_cl2.txt",
			"--input", f"../cluster.{self.sampcon_cluster}/cluster_center_model.rmf3",
			"--output", "./patch_colored_cluster_center_model.rmf3"]
			)
		os.rename( "./patch_colored_cluster_center_model.rmf3", f"./output{self.sampcon_cluster}/patch_colored_cluster_center_model.rmf3" )


if __name__ == "__main__":
	IMPAnalysis().forward()


