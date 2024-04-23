#!/usr/bin/env python
import unittest
import os
import shutil
import sys
import subprocess

TOPDIR = os.path.abspath( os.path.join(os.path.dirname( sys.argv[0]), "../" ) )


class Tests( unittest.TestCase ):
    def test_sims( self ):
        """Test simulations for HDAC corepressor complexes (NuRD, coREST, SIn3A)"""
        for cplex in ["NuRD", "coREST", "Sin3A"]:
            # os.chdir( f"{TOPDIR}/scripts/modeling/{cplex}" )
            os.chdir( f"../scripts/modeling/{cplex}/" )

            input_dir = f"../../../input/{cplex}/"
            if cplex == "NuRD":
                script = "hdac1_modeling.py"
                # Need to unzip the compressed mrc file.
                if not os.path.exists( f"{input_dir}/GMMs/EMD_22895.mrc.gmm.15.mrc" ):
                    subprocess.call( ["unzip", f"{input_dir}/GMMs/EMD_22895.mrc.gmm.15.zip",
                                        "-d", f"{input_dir}/GMMs/"] )
            
            elif cplex == "coREST":
                script = "coREST_hdac1_modeling.py"
                # Need to unzip the compressed mrc file.
                if not os.path.exists( f"{input_dir}/GMMs/emd_10627.mrc.gmm.15.mrc" ):
                    subprocess.call( ["unzip", f"{input_dir}/GMMs/emd_10627.mrc.gmm.15.zip",
                                    "-d", f"{input_dir}/GMMs/"] )
            
            elif cplex == "Sin3A":
                script = "sin3a_hdac1_modeling.py"

            else:
                raise Exception( "No such complex modeled..." )

            p = subprocess.check_call(
                ["python", f"{script}", "test", "test"]
            )

            # Require that the number of frames is present.
            total_num_lines_stat_files = 0
            for i in range( 1 ):
                with open( "Output_test/stat." + str( i ) + ".out", "r" ) as statf:
                    total_num_lines_stat_files += len( statf.readlines() )
            self.assertEqual( total_num_lines_stat_files, 11 )

            # Require that output files were produced
            for i in range(1):
                os.unlink( "Output_test/rmfs/" + str( i ) + ".rmf3" )
                os.unlink( "Output_test/stat." + str( i ) + ".out" )
                os.unlink( "Output_test/stat_replica." + str( i ) + ".out" )
            
            # Clean up.
            if cplex == "NuRD":
                print( "Removing NuRD mrc file...\n" )
                subprocess.call( ["rm", "-r",
                                    "Output_test",
                                    f"{input_dir}/GMMs/EMD_22895.mrc.gmm.15.mrc",
                                    "excluded.None.xl.db", 
                                    "included.None.xl.db", 
                                    "missing.None.xl.db"] )

            elif cplex == "coREST":
                subprocess.call( ["rm", "-r",
                                    "Output_test",
                                    f"{input_dir}/GMMs/emd_10627.mrc.gmm.15.mrc",
                                    "excluded.None.xl.db", 
                                    "included.None.xl.db", 
                                    "missing.None.xl.db"] )
            elif cplex == "Sin3A":
                subprocess.call( ["rm", "-r",
                                    "Output_test",
                                    "excluded.None.xl.db", 
                                    "included.None.xl.db", 
                                    "missing.None.xl.db"] )
            # back to base.
            os.chdir( f"../../../test" )
            print( "\n------------------------------------------------\n" )


if __name__ == "__main__":
    unittest.main()
