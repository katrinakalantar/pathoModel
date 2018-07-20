import pandas as pd
import sys, os.path
import argparse
import pandas as pd
import sys
import sklearn
from sklearn import metrics
import numpy as np
import pickle

parser = argparse.ArgumentParser()
parser.add_argument("filename", help="name input IDSeq table")
parser.add_argument("pathogens", help="name of reference pathogen list")
parser.add_argument("model", help="logistic regression model file")
args = parser.parse_args()


scripts_dir = './scripts/'
sys.path.append(scripts_dir)


def convert_input(df):
    #print("Input counts:" + str(df.shape))
    df = df[df.tax_id > 0]

    # added 7/5 to remove phage
    df = df[df.family_taxid != -300]
    df = df[df.is_phage != 1]

    genus_df = df[df.tax_level == 2]
    species_df = df[df.tax_level < 2]

    genus_df = genus_df[['name', 'genus_taxid', 'NT_rpm', 'NT_zscore']]
    genus_df.columns = ['genus_name', 'genus_taxid', 'genus_NT_rpm', 'genus_NT_zscore']
    
    merged = pd.merge(species_df, genus_df, how='left', on = 'genus_taxid')
    return(merged)


def calculate_scores(df_rna, df_dna):

	keep_kingdoms = ['Bacteria','Fungi','Viruses','Escherichia coli']

	### prepare the RNA dataframe
	df_rna = df_rna[df_rna['category_name'].isin(keep_kingdoms)]  #remove kingdoms not in the list of kingdoms
	
	# RNA must have non-zero expression for all microbes, including DNA viruses;
	# require concordance on NR and NT at genus level
	df_rna = df_rna[df_rna['NR_rpm'] > 0]
	df_rna = df_rna[df_rna['NT_rpm'] > 0]
	
	df_rna = df_rna[['category_name','genus_taxid','NT_rpm','NT_zscore','genus_NT_rpm','genus_NT_zscore','name','NR_rpm']]   #'NT Species rM', subset the original dataframe to contain only the relevant columns
	idx = df_rna.groupby(['genus_taxid'])['NT_rpm'].transform(max) == df_rna['NT_rpm']				  # group by Genus, collapsing to take the species with the greatest species-level rM
	df_rna = df_rna[['category_name','genus_taxid','NT_rpm','NT_zscore','genus_NT_rpm','genus_NT_zscore','name','NR_rpm']][idx]			   # keep the species with max rM for each genus
	df_rna = df_rna[df_rna['NT_zscore'] > 1.0]  # ADDED 7/05 to filter out SPECIES z < 1 
	df_rna.columns = ['category_name','genus_taxid','NT_rpm RNA','NT_zscore RNA','genus_NT_rpm RNA','genus_NT_zscore RNA','name','NR_rpm RNA']		# re-name the columns
	
	if df_dna != None:
		
		### prepare the DNA dataframe
		df_dna = df_dna[df_dna['category_name'].isin(keep_kingdoms)]

		df_rna = df_dna[df_dna['NR_rpm'] > 0]
		df_rna = df_dna[df_dna['NT_rpm'] > 0]
		
		df_dna = df_dna[['category_name','genus_taxid','NT_rpm','NT_zscore','genus_NT_rpm','genus_NT_zscore','name', 'NR_rpm']]   #'NT Species rM', subset the original dataframe to contain only the relevant columns
		idx = df_dna.groupby(['genus_taxid'])['NT_rpm'].transform(max) == df_dna['NT_rpm']				  # group by Genus, collapsing to take the species with the greatest species-level rM
		df_dna = df_dna[['category_name','genus_taxid','NT_rpm','NT_zscore','genus_NT_rpm','genus_NT_zscore','name','NR_rpm']][idx]			   # keep the species with max rM for each genus
		df_dna = df_dna[df_dna['NT_zscore'] > 1.0]  # ADDED 7/05 to filter out SPECIES z < 1 
		df_dna.columns = ['category_name','genus_taxid','NT_rpm DNA','NT_zscore DNA','genus_NT_rpm DNA','genus_NT_zscore DNA','name','NR_rpm DNA']		# re-name the columns

		if(type(subtract_dna) == type(pd.Series())):
			#subtract_dna['genus_taxid'] = subtract_dna.index
			subtract_dna = pd.Series.to_frame(subtract_dna)
			subtract_dna['genus_taxid'] = list(subtract_dna.index)
			df_dna_temp = df_dna.merge(pd.DataFrame(subtract_dna), how = 'left', on='genus_taxid')
			df_dna_temp['subtracted_values'] = df_dna_temp['genus_NT_rpm DNA'] - df_dna_temp[0]
			df_dna_temp['genus_NT_rpm DNA'] = df_dna_temp['subtracted_values']
			df_dna_temp.drop([0, 'subtracted_values'], axis=1, inplace=True)
			df_dna = df_dna_temp
			df_dna = df_dna[df_dna['genus_NT_rpm DNA'] > 0]

		# write to file for manual verification
		df_dna.sort_values(by='NT_rpm DNA', inplace=True, ascending=False)   # sort by Genus rM

		# merge DNA and RNA into one matrix with DNA rM and RNA rM values
		master_df = pd.merge(df_dna, df_rna, how='outer', on=['category_name','genus_taxid'])   # merge on Genus ID
		master_df.fillna(value=0, inplace=True)	# convert NA values -> 0

		master_df.drop_duplicates(subset=['genus_taxid'],inplace=True)

		# Filter for presence in both RNA and DNA-seq
		no_vir = master_df[master_df['category_name']!='Viruses']		   # get the subset of microbes that are not viruses (bacterial/fungal rows)
		no_vir_rna = no_vir[no_vir['NT_rpm RNA'] > 0]			 # for bacterial/fungal rows, remove microbes with NT_rpm RNA = 0
		no_vir_dna = no_vir_rna[no_vir_rna['NT_rpm DNA'] > 0]	 # for bacterial/fungal rows, remove microbes with NT_rpm DNA = 0
		no_vir_dna = no_vir_dna[no_vir_dna['NR_rpm RNA'] > 0]			 # added 2/13 - for bacterial/fungal rows, remove microbes with NT_rpm RNA = 0
		no_vir_dna = no_vir_dna[no_vir_dna['NR_rpm DNA'] > 0]	 # added 2/13 - for bacterial/fungal rows, remove microbes with NT_rpm DNA = 0
		v = master_df[master_df['category_name'] == 'Viruses']			  # get the subset of microbes that are viruses (no restrictions on RNA and DNA presence)

		dna_viruses_list = ['Mastadenovirus ( 10509 )','Iridovirus ( 10487 )','Roseolovirus ( 40272 )','Betapapillomavirus ( 333922 )','Molluscipoxvirus ( 10278 )','Gammapapillomavirus ( 325455 )',
							'Alphatorquevirus ( 687331 )','Cytomegalovirus ( 10358 )','Simplexvirus ( 10294 )','Lymphocryptovirus ( 10375 )','Polyomavirus ( 10624 )','Varicellovirus ( 10319 )']
		dna_virs = [i for i in range(len(v.index)) if v.iloc[i]['genus_taxid'] in dna_viruses_list ]
		not_dna_virs = [i for i in range(len(v.index)) if v.iloc[i]['genus_taxid'] not in dna_viruses_list ]
		v_dna_vir = v.iloc[dna_virs]
		v_rna_vir = v.iloc[not_dna_virs]
		v_dna_vir = v_dna_vir[v_dna_vir['NT_rpm DNA'] > 0]
		v_dna_vir = v_dna_vir[v_dna_vir['NT_rpm RNA'] > 0]
		v_dna_vir = v_dna_vir[v_dna_vir['NR_rpm DNA'] > 0]   # added 2/13
		v_dna_vir = v_dna_vir[v_dna_vir['NR_rpm RNA'] > 0]   # added 2/13
		
		#temp = pd.concat([no_vir_dna,v], axis=0)					   # concatonate filtered bacterial/fungal with viral reads
		temp = pd.concat([no_vir_dna,v_dna_vir, v_rna_vir])
		master_df = temp

		# create separate output DataFrames for each bacterial, viral, and combined
		bacterial_df = master_df[master_df['category_name'] != 'Viruses']
		viral_df = master_df[master_df['category_name'] == 'Viruses']
		combined = master_df

		combined.sort_values(by='NT_rpm RNA', inplace=True, ascending=False)  # sort combined data by NT rM RNA value
		

		assigned_species_names = []
		for i in combined.index:
			curr = combined.loc[i]['name_x'] # default is DNA
			if(combined.loc[i]['genus_NT_rpm DNA'] < combined.loc[i]['genus_NT_rpm RNA']): #if DNA < RNA
				curr = combined.loc[i]['name_y']
			assigned_species_names.append(curr)
		combined['Species_Assignment'] = assigned_species_names

	else: #df_dna is none...need to work with just RNA

		master_df = df_rna
		master_df.fillna(value=0, inplace=True)	# convert NA values -> 0
		master_df.drop_duplicates(subset=['genus_taxid'],inplace=True)

		# create separate output DataFrames for each bacterial, viral, and combined
		bacterial_df = master_df[master_df['category_name'] != 'Viruses']
		viral_df = master_df[master_df['category_name'] == 'Viruses']
		combined = master_df
		combined.columns = ['category_name','genus_taxid','NT_rpm RNA','NT_zscore RNA','genus_NT_rpm RNA','genus_NT_zscore RNA','name_y','NR_rpm RNA']		# re-name the columns
	

		combined.sort_values(by='NT_rpm RNA', inplace=True, ascending=False)  # sort combined data by NT rM RNA value

		assigned_species_names = []
		for i in combined.index:
			curr = combined.loc[i]['name_y']
			assigned_species_names.append(curr)
		combined['Species_Assignment'] = assigned_species_names
		#combined[['category_name','genus_taxid','Species_Assignment','NT_rpm RNA','NR_rpm RNA','NT_zscore RNA']].head(n=20).to_csv(output_file, sep='\t', mode='a',index=False)

	return [bacterial_df, viral_df, combined]




v = open(args.pathogens, 'r').readlines()
list_of_viruses_in_dictionary = [i.strip() for i in v]

# Curated list of top respiratory pathogens
resp_mic = open(args.pathogens, 'r').readlines()
full_respiratory_microbes_list = [i.strip() for i in resp_mic]


df_rna = pd.read_csv(args.filename, error_bad_lines=False)
df_rna = convert_input(df_rna)


scores = calculate_scores(df_rna, None)


bac = scores[2]
if(bac.shape[0] > 0):
	bac = bac.reset_index()

	# sort the operating matrix by column RPM
	bac.sort_values(by='genus_NT_rpm RNA', inplace=True, ascending=False)   
	bac = bac.reset_index() #reindex([i for i in range(len(bac.index))])


	bac = bac.head(n=15) # only consider the top 15 most abundant  
	bac.fillna(0)

	grew_in_culture = []   # fill color
	is_pathogen = []	   # edge color
	is_virus = []		  # shape (square = virus, circle = bacteria/fungi)

	bac = bac.dropna(axis=0, how='all')

	rank_number_identified = 0


	RNA_value = []
	DNA_value = []
	pathogenic = []
	virus = []
	microbe_id = []
	microbe_id_genus = []
	ranks = []

	for n in bac.index:
		rank_number_identified += 1   # increment with each microbe

		curr_genus = None
		species_colID = None

		if 'name_x' in bac.columns:				
			# implementing a heuristic for assigning most likely 
			# species to genus-level rpM
			curr_genus = bac.loc[n]['name_x']	# default is DNA species...
			species_colID = 'name_x'
						
			# ...but if RNA rM is greater, then set curr_genus to the RNA species
			if bac.loc[n]['genus_NT_rpm DNA'] < bac.loc[n]['genus_NT_rpm RNA']:   
				curr_genus = bac.loc[n]['name_y']
				species_colID = 'name_y'
		else:
			curr_genus = bac.loc[n]['name_y']
			species_colID = 'name_y'

		#
		# set the facecolor of points based on classificaiton of microbe. Overall... 
		# blue = non-pathogenic bacteria/fungi
		# green = non-pathogenic virus
		# red = pathogenic (bacteria, virus, or fungi)
		#
					
		set_color = "darkblue"	 # default blue
		if (curr_genus in list_of_viruses_in_dictionary) or ('irus' in str(curr_genus)):  #if genus is any virus, green
			set_color = "green"
			is_virus.append("green")   
		else:
			is_virus.append('none')
			# set red if pathogen, regardless of whether it was previously 
			# set to green for "virus"; this means that viruses on list of pathogens
			# will appear red, not green.

		try:
			if (curr_genus in full_respiratory_microbes_list) or (curr_genus.split(' ')[0] in full_respiratory_microbes_list):  
				set_color = "red"
		except:
			if (curr_genus in full_respiratory_microbes_list):
				set_color = "red"
					
		# to the vector of pathogenicity, append the color; red == pathogen
		is_pathogen.append(set_color)  # "red") ADDED 7/6 to make all microbes "pathogens" 
			
		g = bac.loc[n]['genus_taxid']   # genus-level microbe ID



	correct_species_names = []
	for n in bac.index:
		species_colID = None
		if 'name_x' in bac.columns:				
			# implementing a heuristic for assigning most likely 
			# species to genus-level rpM
			curr_genus = bac.loc[n]['name_x']	# default is DNA species...
			species_colID = 'name_x'
						
			# ...but if RNA rM is greater, then set curr_genus to the RNA species
			if bac.loc[n]['genus_NT_rpm DNA'] < bac.loc[n]['genus_NT_rpm RNA']:   
				curr_genus = bac.loc[n]['name_y']
				species_colID = 'name_y'
		else:
			curr_genus = bac.loc[n]['name_y']
			species_colID = 'name_y'

		correct_species_names.append(bac.loc[n][species_colID])


	#for sensitivity analysis with real model 
	RNA_value = RNA_value + list(np.log10((bac['genus_NT_rpm RNA']+1)))  
	#DNA_value = DNA_value + list(np.log10((bac['genus_NT_rpm DNA']+1)))  
	pathogenic = pathogenic + is_pathogen
	virus = virus + is_virus
	microbe_id = microbe_id  + correct_species_names
	microbe_id_genus = microbe_id_genus + [bac['genus_taxid'][i] for i in range(len(bac['genus_taxid']))]
	ranks = ranks + [j for j in range(len(bac.index))]
					

# this is the matrix that we will be using for predictions
sensitivity_DF = pd.DataFrame.from_dict({'RNAvalue':RNA_value,
										 'pathogenic':pathogenic,
										 'ranks':ranks,
										 'virus':virus,
										 'microbe':microbe_id,
										 'microbe_genus':microbe_id_genus}, orient='columns')


X = sensitivity_DF[['RNAvalue','ranks','microbe','microbe_genus']]
X[['pathogenic_red']] = sensitivity_DF[['pathogenic']] == 'red'

lr = pickle.load( open( args.model, "rb" ) )
training_variables = ['RNAvalue','pathogenic_red','ranks'] #'DNAvalue','pathogenic_green',
PROBABILITY_THRESHOLD = .2

predicted = lr.predict_proba(X[training_variables])
X['score'] = predicted[:,1]
print(X)


