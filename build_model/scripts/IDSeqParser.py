import pandas as pd
import numpy as np


# matrix transformations

def log_transform_rpm(matrix):
	""" log transform rpm values """
    
	matrix['log(NR_rpm)'] = [np.log10(i + 1) for i in matrix['NR_rpm']]
	matrix['log(NT_rpm)'] = [np.log10(i + 1) for i in matrix['NT_rpm']]
	matrix['log(genus_NR_rpm)'] = [np.log10(i + 1) for i in matrix['genus_NR_rpm']]
	matrix['log(genus_NT_rpm)'] = [np.log10(i + 1) for i in matrix['genus_NT_rpm']]
    
def set_percent(matrix):
	""" get the percentage of non-host reads represented by this pathogen via non-host read value in metadata """

	nonhost_rpm = matrix['nonhost_reads_percent'] * 10000
	matrix['perc_NRr'] = matrix['NR_rpm'] / nonhost_rpm
	matrix['perc_NTr'] = matrix['NT_rpm'] / nonhost_rpm
	matrix['perc_genus_NRr'] = matrix['genus_NR_rpm'] / nonhost_rpm
	matrix['perc_genus_NTr'] = matrix['genus_NT_rpm'] / nonhost_rpm

	#matrix['perc_NRr'] = matrix['NR_r'] / matrix['nonhost_reads']
	#matrix['perc_NTr'] = matrix['NT_r'] / matrix['nonhost_reads']
	#matrix['perc_genus_NRr'] = matrix['genus_NR_r'] / matrix['nonhost_reads']
	#matrix['perc_genus_NTr'] = matrix['genus_NT_r'] / matrix['nonhost_reads']
    
def set_nucl_type(matrix):
	""" indicate whether this is RNA- or DNA-seq (by parsing sample ID """
	matrix['RNA'] = ['RNA' in i for i in matrix['sampleID']]
    
    
def set_percent_local(matrix):
	""" set the percentage of non-host reads represented by this pathogen via sum of all pathogens """
	
	matrix['perc_NR2'] = None
	matrix['perc_NT2'] = None

	for i in list(set(matrix['sampleID'])):
		total_NR = matrix[matrix['sampleID']==i]['NR_r'].sum()
		total_NT = matrix[matrix['sampleID']==i]['NT_r'].sum()
        
		for j in matrix[matrix['sampleID']==i].index:
			matrix.loc[j,'perc_NR2'] = matrix.loc[j,'NR_r']/total_NR
			matrix.loc[j,'perc_NT2'] = matrix.loc[j,'NT_r']/total_NT


def create_idseq_project(project_directory, metadata_file, output_directory, reference_pathogen_file):
	""" Check that the project directory contains required data, parse into microbe form """

	merged_metadata = parse_metadata(metadata_file)

	list_of_microbe_matrices = []
	for i in merged_metadata.index:

		parsed_data_matrix = parse_data(i, project_directory)

		try:
			full_matrix = pd.merge(parsed_data_matrix, merged_metadata, how = 'left', on = 'sampleID')
			set_positivity(full_matrix)
			set_pathogenicity(full_matrix, reference_pathogen_file)
			list_of_microbe_matrices.append(full_matrix)
		except:
			continue


	all_microbes = pd.concat(list_of_microbe_matrices, ignore_index = True)


	# Full training dataset exists, do some post-modification of values
	all_microbes.fillna(0,inplace=True)
	log_transform_rpm(all_microbes)
	set_percent(all_microbes)
	set_nucl_type(all_microbes)
	set_percent_local(all_microbes)   #### NOTE!!! HAVE NOT TESTED THAT THIS IS WORKING TOTALLY CORRECTLY
	all_microbes['is_virus'] = all_microbes['category_name'] == "Viruses"

	# write final output to file
	all_microbes.to_csv(output_directory  + "/microbe_data.csv")
	return(all_microbes)


def parse_notes(note):
	""" Read in the "Notes" column from individual samples in IDSeq metadata 
	    and parse into separate fields """

	_split_note_ = note.splitlines()#('\r')
	_parsed_note_ = {}

	extra_info = []

	for i in _split_note_:
		try:
			_parsed_note_[i.split(':')[0].replace("-","").strip()] = i.split(':')[1].strip()
		except:
			extra_info.append(i.strip())

	_parsed_note_['extra_info'] = '-'.join(extra_info)

	return(_parsed_note_)

def parse_metadata(metadata_file, relevant_columns = ['total_reads', 'nonhost_reads', 'nonhost_reads_percent', 'compression_ratio', 'tissue_type', 'nucleotide_type', 'known_organisms','effective_group']):
	""" Loop through the input metadata file to extract relevant info """

	metadata = pd.read_csv(metadata_file, sep = '\t', index_col = 0, header = 0)

	if 'notes' in metadata.columns:
		notes_parsed = {}
		for sample in metadata.index:

			try:
				#if True:
				notes_parsed[sample] = parse_notes(metadata.loc[sample]['notes'])
			except Exception as e:
				print("WARNING : " + str(e) +" : " + sample + 
					" - failed to parse metadata, check that this sample contains data in the notes section")
				notes_parsed[sample] = metadata.loc[sample]['notes']
				
		notes_parsed_df = pd.DataFrame(notes_parsed).transpose()
		merged_metadata = metadata[list(set(relevant_columns).intersection(set(metadata.columns)))].join(notes_parsed_df, how = 'outer')
	else:
		merged_metadata = metadata[list(set(relevant_columns).intersection(set(metadata.columns)))]
	merged_metadata['sampleID'] = merged_metadata.index
	return(merged_metadata)


def parse_data(samplename, project_directory):

	""" Read IDSeq .csv report file and parse the relevant microbes out, appending sample info to each microbe """



	# input variables
	filename = project_directory + "/" + samplename.lower() + ".csv"
	keep_kingdoms = ['Bacteria','Fungi','Viruses','Escherichia coli','Eukaryota']

	# filter data based on thresholds / qualifications pre-established
	try:
		data = pd.read_csv(filename, sep = ",")
	except pd.errors.EmptyDataError as E:
		print("ERROR   : " +str(E) + " : " + samplename)
		return None
	#data = data[data.tax_id  > 0]  
	data = data[data.category_name.isin(keep_kingdoms)]  #remove kingdoms not in the list of kingdoms
	data = data[data.family_taxid != -300]
	data = data[data.is_phage != 1]

	# split into genus and species rows
	genus_df = data[data.tax_level == 2]
	species_df = data[data.tax_level < 2]

	species_df = species_df[species_df.species_taxid > 0]  # added 7/30 to remove "non-species specific..." as the most common species ID

	# rename the genera and merge lines with species data
	genus_df = genus_df[['name', 'genus_taxid', 'NT_rpm', 'NT_zscore', 'NR_rpm', 'NR_zscore', 'NT_r', 'NR_r']]
	genus_df.columns = ['genus_name', 'genus_taxid', 'genus_NT_rpm', 'genus_NT_zscore', 'genus_NR_rpm', 'genus_NR_zscore', 'genus_NT_r', 'genus_NR_r']
	merged = pd.merge(species_df, genus_df, how='left', on = 'genus_taxid')

	# filters prior to selecting most abundant species 
	merged = merged[merged.NT_rpm > 0]
	merged = merged[merged.NR_rpm > 0]

	# added 7/31 to keep low-level viruses but not dilute bacterial signal
	#viruses = merged[merged.category_name == 'Viruses']
	#bacteria = merged[merged.category_name != 'Viruses']
	#bacteria.sort_values(by='NT_rpm', inplace = True, ascending = False)
	#merged = pd.concat([viruses, bacteria.head(n=15)])

	# OVERRIDE cases where the genus doesn't match (ie genus = -200 for S. aureus in sample ID 212)
	idx_fix = merged[pd.isnull(merged['genus_name'])].index
	for j in idx_fix:
		merged.loc[j, 'genus_name'] = str(merged.loc[j, 'name']).split()[0]
		merged.loc[j, 'genus_NT_rpm'] = merged.loc[j, 'NT_rpm'] # override genus rpm (NaN) with species rpm


	# select the most abundant species per genus
	merged.sort_values(by='NT_rpm', inplace=True, ascending = False)
	idx = merged.groupby(['genus_name'])['NT_rpm'].transform(max) == merged['NT_rpm']  # group by Genus, collapsing to take the species with the greatest species-level rM
	merged = merged.loc[idx]   # keep the species with max rM for each genus
	
	# filters after selecting most abundant species
	#merged = merged[merged.NT_zscore > .5]
	#merged = merged[merged.NR_zscore > .5]

	merged['sampleID'] = samplename

	merged.sort_values(by='NT_rpm', inplace=True, ascending = False)
	merged['rank'] = [i for i in range(len(merged.index))]
	return(merged.head(n=15))


def set_positivity(matrix):
	""" Take in matrix of pathogens and append column for "positive" indicating if this microbe was identified by standard clinical microbiology """
	matrix['positive'] = False  # default
	for i in matrix.index:
		
		try:
			known_organisms = matrix.loc[i, 'known_organisms'].split('-')
			positive = (matrix.loc[i,'genus_name'] in known_organisms) or (matrix.loc[i,'name'] in known_organisms)
			matrix.loc[i, 'positive'] = positive
		except AttributeError:  # no known organism for this patient (known_organisms == nan)
			continue


def set_pathogenicity(matrix, filename):
	""" Take in matrix of pathogens and append column for "reference_pathogen" indicating if this microbe was in the list of reference pathogens """
	reference_pathogens = open(filename,'r').readlines()
	ref_pathogens = [l.strip() for l in reference_pathogens]
	for i in matrix.index:
		pathogen = (matrix.loc[i,'genus_name'] in ref_pathogens) or (matrix.loc[i,'name'] in 
			ref_pathogens) or (matrix.loc[i,'name'].split(' ')[0] in 
			ref_pathogens) or (' '.join(matrix.loc[i,'name'].split(' ')[0:2]) in ref_pathogens)
		matrix.loc[i, 'reference_pathogen'] = pathogen


## Run the script ## 

metadata_filename = '/Users/kkalantar/Downloads/project-rapid_response_007_sample-table_BN.tsv'
project_directory_name = '/Users/kkalantar/Desktop/pathoModel/data/070518/rapid-response-007_reports'
output_directory_name = '/Users/kkalantar/Desktop/pathoModel/build_model/output/080218_BN_2'
reference_pathogen_filename = '/Users/kkalantar/Desktop/pathoModel/reference/pathogens_bangladesh_official.txt'

'''
metadata_filename = '/Users/kkalantar/Desktop/pathoModel/build_model/benchmark_data/project-mbal_study_sample-table_073018.csv'
project_directory_name = '/Users/kkalantar/Desktop/pathoModel/build_model/benchmark_data/mbal_study_reports_073018'
output_directory_name = '/Users/kkalantar/Desktop/pathoModel/build_model/output/073018'
reference_pathogen_filename = '/Users/kkalantar/Desktop/pathoModel/build_model/reference/known_respiratory_pathogens.txt'
'''


create_idseq_project(project_directory_name, metadata_filename, output_directory_name, reference_pathogen_filename)


