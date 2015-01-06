#! /usr/bin/env python
""" 
    Routines to performe the analysis of virtual screening
    These routines were developed by:
    Rodrigo Antonio Faccioli - rodrigo.faccioli@usp.br / rodrigo.faccioli@gmail.com  
    Leandro Oliveira Bortot  - leandro.bortot@usp.br / leandro.obt@gmail.com 
"""

import analysis 
import analysisio as ana_io
import xvg_histogram_energy_values as xvghist

def call_vs_analysis(path_analysis, path_log):
	"""
	Executes the analysis of Virtual Screening. 
	The analysis consists of: 
	1) Creating a txt file that contains receptor, ligand and mode sorted by energies.
	2) Creating a xvg file that is a histogram of energies	
    Example:
        >>> call_vs_analysis(path_analysis, path_log)
    @param path_analysis: place where analysis files will be saved
    @type path_analysis: string        
    @param path_log: place where log files were saved
    @type path_log: string            
	"""		
	log_dict = analysis.log_files_by_energy( path_log )

	log_sorted_dict = ana_io.create_file_by_sorted_energy(path_analysis, log_dict) 

	xvghist.create_xvg_histogram_energy_values(path_analysis, log_sorted_dict)	
	