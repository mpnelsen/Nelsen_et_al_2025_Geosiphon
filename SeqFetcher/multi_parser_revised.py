#! /usr/bin/env python
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os, sys
loci=["nuSSU"]
#pathtofolders="..."
pathtofolders=".../accs"
locus_genenames={}
locus_genenames["nuSSU"]=["18S ribosomal RNA","small subunit ribosomal RNA","18S small subunit ribosomal RNA", "small subunit 18S ribosomal RNA", "sequence contains 18S rRNA gene", "18S rRNA", "nuclear encoded 18S ribosomal RNA, small subunit", "nuclear encoded small subunit ribosomal RNA", "small subunit (18S) ribosomal RNA", "18S ribosomal RNA", "16S small subunit ribosomal RNA", "18S small ribosomal RNA subunit", "nuclear small subunit ribosomal RNA", "contains 18S ribosomal RNA, internal transcribed spacer 1, 5.8S ribosomal RNA, internal transcribed spacer 2, and 28S ribosomal RNA", "16S ribosomal RNA", "SSU ribosomal RNA","contains 18S ribosomal RNA, internal transcribed spacer 1, 5.8S ribosomal RNA, internal transcribed spacer 2, and 26S ribosomal RNA", "contains 18S ribosomal RNA, internal transcribed spacer 1, 5.8S ribosomal RNA, internal transcribed spacer 2, and 28S ribosomal RNA", "small subunit ribosomal RNA","18S SSU ribosomal RNA", "18S small subunit ribosomal RNA", "ribosomal RNA small subunit","contains 18S ribosomal RNA, internal transcribed spacer 1, 5.8S ribosomal RNA, internal transcribed spacer 2, and 28S ribosomal RNA","small ribosomal RNA subunit RNA","SSU ribosomal RNA","18S rRNA","nuclear 18S ribosomal RNA","18S ribosomal RNA gene","put. 18S ribosomal RNA","18S ribosomal RNA, small subunit","ribosomal RNA; small subunit ribosomal RNA","genomic DNA"]

def parse(file):
	for record in SeqIO.parse(open(file,"r"),"genbank"):
		accession=record.id
		species=record.annotations['organism']
		species2=species.replace(" ","_")
		species3=species2.replace(".","")
		if locus=="nuSSU":
			print locus
			for i, feature in enumerate (record.features):
				if locus=="nuSSU":
				#print('{0}_genenames'.format(locus))
				#genenames=["18S ribosomal RNA","small subunit ribosomal RNA","18S small subunit ribosomal RNA", "small subunit 18S ribosomal RNA", "sequence contains 18S rRNA gene", "18S rRNA", "nuclear encoded 18S ribosomal RNA, small subunit", "nuclear encoded small subunit ribosomal RNA", "small subunit (18S) ribosomal RNA", "18S ribosomal RNA", "16S small subunit ribosomal RNA", "18S small ribosomal RNA subunit", "nuclear small subunit ribosomal RNA"]
					if any(feature.type=='rRNA' for i, feature in enumerate (record.features)):
						if feature.type=='rRNA':
							if 'product' in feature.qualifiers:
								geneofinterest=feature.qualifiers['product'][0]
								seq = feature.extract(record.seq)
								if geneofinterest in locus_genenames[locus]:
									print '>%s %s\n%s' % (species3, geneofinterest, seq)
									parsedseqfile.write('>{0}\n{1}\n'.format(accession,seq))
									acclistfile.write('%s\t%s\n' % (species3,accession))
									break
							else:
								if 'note' in feature.qualifiers:
									geneofinterest=feature.qualifiers['note'][0]
									seq = feature.extract(record.seq)
									if geneofinterest in locus_genenames[locus]:
										print '>%s %s\n%s' % (species3, geneofinterest, seq)
										parsedseqfile.write('>{0}\n{1}\n'.format(accession,seq))
										acclistfile.write('%s\t%s\n' % (species3, accession))
										break		
					else:
						if any(feature.type=='misc_feature' for i, feature in enumerate (record.features)):
							if feature.type=='misc_feature':
								if 'note' in feature.qualifiers:
									geneofinterest=feature.qualifiers['note'][0]
									seq = feature.extract(record.seq)
									if geneofinterest in locus_genenames[locus]:
										print '>%s %s\n%s' % (species3, geneofinterest, seq)
										parsedseqfile.write('>{0}\n{1}\n'.format(accession,seq))
										acclistfile.write('%s\t%s\n' % (species3, accession))
										break		
						else:
							if any(feature.type=='misc_RNA' for i, feature in enumerate (record.features)):
								if feature.type=='misc_RNA':
									if 'product' in feature.qualifiers:
										geneofinterest=feature.qualifiers['product'][0]
										seq = feature.extract(record.seq)
										if geneofinterest in locus_genenames[locus]:
											print '>%s %s\n%s' % (species3, geneofinterest, seq)
											parsedseqfile.write('>{0}\n{1}\n'.format(accession,seq))
											acclistfile.write('%s\t%s\n' % (species3, accession))
											break		
									else:
										if 'note' in feature.qualifiers:
											geneofinterest=feature.qualifiers['note'][0]
											seq = feature.extract(record.seq)
											if geneofinterest in locus_genenames[locus]:
												print '>%s %s\n%s' % (species3, geneofinterest, seq)
												parsedseqfile.write('>{0}\n{1}\n'.format(accession,seq))
												acclistfile.write('%s\t%s\n' % (species3,accession))
												break
							else:
								if any(feature.type=='precursor_RNA' for i, feature in enumerate (record.features)):
									if feature.type=='precursor_RNA':
										if 'product' in feature.qualifiers:
											geneofinterest=feature.qualifiers['product'][0]
											seq = feature.extract(record.seq)
											if geneofinterest in locus_genenames[locus]:
												print '>%s %s\n%s' % (species3, geneofinterest, seq)
												parsedseqfile.write('>{0}\n{1}\n'.format(accession,seq))
												acclistfile.write('%s\t%s\n' % (species3,accession))
												break
										else:
											if 'note' in feature.qualifiers:
												geneofinterest=feature.qualifiers['note'][0]
												seq = feature.extract(record.seq)
												if geneofinterest in locus_genenames[locus]:
													print '>%s %s\n%s' % (species3, geneofinterest, seq)
													parsedseqfile.write('>{0}\n{1}\n'.format(accession,seq))
													acclistfile.write('%s\t%s\n' % (species3, accession))
													break		
								else:
									if any(feature.type=='source' for i, feature in enumerate (record.features)):					
										if feature.type=='source':
											if 'mol_type' in feature.qualifiers:
												geneofinterest=feature.qualifiers['mol_type'][0]
												seq = feature.extract(record.seq)
												if geneofinterest in locus_genenames[locus]:
													print '>%s %s\n%s' % (species3, geneofinterest, seq)
													parsedseqfile.write('>{0}\n{1}\n'.format(accession,seq))
													acclistfile.write('%s\t%s\n' % (species3, accession))
													break		
											else:
												break		


for locus in loci:
	if locus=="nuSSU":
		parsedseqfile = open("{0}/{1}/{1}_Parsed.fasta".format(pathtofolders,locus),"wa")
		acclistfile = open("{0}/{1}/{1}_accession_list.txt".format(pathtofolders,locus),"wa")
		parse("{0}/{1}/{1}.gb".format(pathtofolders,locus))
		parsedseqfile.close()
		acclistfile.close();
		