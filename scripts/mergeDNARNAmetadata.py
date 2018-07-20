import pandas as pd
import sys
import math
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-RNA", help="name of RNA metadata file")
parser.add_argument("-DNA", help="name of DNA metadata file")
parser.add_argument("-o", help="name of output metadata file")
parser.add_argument("-parseRNA", help="parse the RNA metadata notes", action = "store_true")
parser.add_argument("-parseDNA", help="parse the DNA metadata notes", action = "store_true")
args = parser.parse_args()


def parse_notes(notes):
	total_notes = 0
	text_df = {}
	for n in notes:

		note_split = n.split("\n")
		total_notes += 1 
		random_text = []
		class_text = {}
		for subnote in note_split:
			s = subnote.split(':')
			try:
				class_text[s[0].replace("-","").strip()] = s[1]
			except:
				random_text.append(subnote)
		class_text['random_text'] = '; '.join(random_text)
		text_df[total_notes] = class_text
	df = pd.DataFrame(text_df).transpose()
	return df
	#print(df)



RNAmetadata = pd.read_csv(args.RNA)

if(args.parseRNA):
	notes = parse_notes(list(RNAmetadata['notes']))
	temp = pd.concat([RNAmetadata.reset_index(drop=True), notes], axis=1)
	RNAmetadata = temp

new_colnames = [i+'_RNA' for i in RNAmetadata.columns]
RNAmetadata.columns = new_colnames
RNA_sampleName = RNAmetadata['sample_name_RNA']
#print(list(RNA_sampleName))
#print([ for i in list(RNA_sampleName)])
RNA_sampleID = ['_'.join([i.split('_')[0],i.split('_')[2]]) if (type(i) != type(1.0)) else '' for i in list(RNA_sampleName)]

RNAmetadata['ID'] = pd.Series(RNA_sampleID, index = RNAmetadata.index)
#
if(args.parseRNA):
	RNAmetadata.drop('notes_RNA', inplace=True, axis=1)

print(RNAmetadata.shape)

#print(RNAmetadata.head())
#print(RNAmetadata.shape)

#parse_notes(list(RNAmetadata['notes_RNA']))


DNAmetadata = pd.read_csv(args.DNA)
DNA_sampleName = DNAmetadata['sample_name']

if(args.parseDNA):
	notes = parse_notes(list(DNAmetadata['notes']))
	temp = pd.concat([DNAmetadata.reset_index(drop=True), notes], axis=1)
	DNAmetadata = temp

new_colnames = [i+'_DNA' for i in DNAmetadata.columns]
DNAmetadata.columns = new_colnames
DNA_sampleID = ['_'.join([i.split('_')[0], i.split('_')[2]]) for i in list(DNA_sampleName)]
DNAmetadata['ID'] = pd.Series(DNA_sampleID, index = DNAmetadata.index)

if(args.parseDNA):
	DNAmetadata.drop('notes_RNA', inplace=True, axis=1)


print(DNAmetadata.shape)
#that 'output.csv' is a SYLK file, but cannot load it. Either the file has 
#merged_metadata = RNAmetadata.join(DNAmetadata.set_index('ID'), on = 'ID')
merged_metadata = pd.merge(RNAmetadata, DNAmetadata, how='outer', on = 'ID')
#merged_metadata = RNAmetadata.append(DNAmetadata, ignore_index=True)
merged_metadata.index = merged_metadata['ID']
#print(merged_metadata.head())
merged_metadata.to_csv(args.o)


