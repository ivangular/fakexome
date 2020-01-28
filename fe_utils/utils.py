#! /usr/bin/python3
#
# This source code is part of fakexome, a simulated exome generator.
#
# Fakexome is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Fakexome is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see<http://www.gnu.org/licenses/>.
#
# Contact: ivana.mihalek@gmail.com, ivana.mihalek@uniri.hr

from fe_utils.CrossMap import *
from fe_utils.mysql import *

#########################################
def translate_positions(positions, chromosome, from_assembly, to_assembly, rootname=None):

	if type(positions)==set: positions = list(positions)
	int_positions = [int(p) for p in positions]
	if from_assembly == to_assembly:
		return dict(zip(int_positions, int_positions))
	# GRCh37 and hg19 only differ for MT
	if (from_assembly.lower() in ['grch37', 'hg19']) and (to_assembly.lower() in ['grch37', 'hg19']) and (chromosome != "MT"):
		return dict(zip(int_positions, int_positions))
	if "grch" in from_assembly.lower():  from_assembly = from_assembly.lower().replace("grch", "GRCh")
	if "grch" in to_assembly.lower():  to_assembly = to_assembly.lower().replace("grch", "GRCh")

	if not rootname: rootname="{}.{}.{}.{}".format(chromosome, from_assembly, to_assembly, str(os.getpid()))
	# otherwise we'll need  tools to translate
	chain_file ="/storage/databases/liftover/{}To{}.over.chain".format(from_assembly, to_assembly.capitalize())
	if not os.path.exists(chain_file):
		print(chain_file, "not found")
		exit()

	chrom = chromosome if 'chr' in chromosome else 'chr'+chromosome
	outfile = "%s.tsv"%rootname
	outf = open (outfile,"w")
	for p in int_positions:
		outf.write("\t".join([chrom, str(p), str(p), str(p)]) + "\n")
	outf.close()

	# this is CrossMap now
	outfile_translated  =  "%s.translated.tsv"%rootname
	(map_tree, target_chrom_sizes, source_chrom_sizes) = read_chain_file(chain_file, print_table=False)
	crossmap_bed_file(map_tree, outfile, outfile_translated)

	#read  regions back in - note that some positions might end upr untranslatable (the unmap file below)
	#others might be on entirely different chromosome - deal with that some other time
	wrongchrom = open(rootname+'.wrong_chrom', "w")
	translation = {}
	with open(outfile_translated,"r") as inf:
		for line in inf:
			if len(line.replace(" ",""))==0: continue
			f = line.rstrip().split("\t")
			if f[0]!=chrom:
				wrongchrom.write(line)
			else:
				translation[int(f[3])] = int(f[2])
	wrongchrom.close()
	# remove aux files
	os.remove(outfile)
	os.remove(outfile_translated)
	if os.stat(outfile_translated+".unmap").st_size != 0:
		print("Warning: {}.unmap in {}".format( outfile_translated, os.getcwd()))
		print("        is not empty - check this file for untranslated positions.")
	else:
		os.remove(outfile_translated+".unmap")
	if os.stat(wrongchrom.name).st_size != 0:
		print("Warning: {} in {}".format( wrongchrom.name, os.getcwd()))
		print("        is not empty - check this file for [psitions translating to antoher chromosome.")
	else:
		os.remove(wrongchrom.name)

	return translation


def read_omim(omim_path):
	if not os.path.exists(omim_path):
		print(omim_path, "not found")
		return {}
	hgnc2omim = {}
	with open(omim_path, "r") as inf:
		for line in inf:
			fields = line.strip().split()
			if len(fields)<2: continue
			[omim, hgnc] = fields[:2]
			hgnc2omim[hgnc] = omim
	return hgnc2omim

###############################
'''
mysql> select count(hgnc_approved) from uniprot_hgnc;
+----------------------+
| count(hgnc_approved) |
+----------------------+
|                20368 |
+----------------------+
mysql> select hgnc_approved, count(uniprot_id) as ct from uniprot_hgnc group by hgnc_approved having ct>1
....
47 rows in set (0.04 sec)
mysql> select count(hgnc_approved) from uniprot_hgnc where uniprot_id is null;
+----------------------+
| count(hgnc_approved) |
+----------------------+
|                    0 |
+----------------------+
1 row in set (0.00 sec)

'''
def hgnc2uniprot(cursor, approved):
	# identifier_maps.uniprot_hgnc has the current uniprot as the key
	# (there are many old uniprot id's that map to the new one in the use around)
	qry  = "select uniprot_id from identifier_maps.uniprot_hgnc "
	qry += "where hgnc_approved='%s'" % approved
	ret = error_intolerant_search(cursor, qry)
	# we are ignoring the fact there might be 2 uniprots
	# the estimate (by couting the intries in the table) is it happens in 0.2% of cases (see above)
	return ret[0][0] if ret else "x"


###############################
'''
388 kegg_id s do not have uniprot id
mysql> select uniprot,  count(kegg_id) as ct from kegg_human group by uniprot having ct>1;
+---------+-----+
| uniprot | ct  |
+---------+-----+
| NULL    | 388 |
+---------+-----+
1 row in set (0.03 sec)
no kegg ids have two uniprots
mysql> select kegg_id, count(uniprot) as ct from kegg_human group by  kegg_id having ct>1;
Empty set (0.01 sec)

'''
def uniprot2kegg(cursor, uniprot):
	qry  = "select kegg_id, kegg_pathways from identifier_maps.kegg_human "
	qry += "where uniprot='%s'" % uniprot
	ret = error_intolerant_search(cursor, qry)
	# we are counting on the fact there might be 2 uniprots (see above)
	return ret[0] if ret else ["x","x"]


###############################
def get_kegg_pthwy_name(cursor, kegg_pathway_id):
	# kegg pathway_id uses padding, as in 00030,  which i chose not to strip:
	qry = "select name from identifier_maps.kegg_pathway_name "
	qry += "where kegg_pathway_id='%s'" % kegg_pathway_id
	ret = error_intolerant_search(cursor, qry)
	return ret[0][0] if ret else "x"


###############################
# flags
HOMOZYGOTE   =  1
COMMON       =  2
EXONIC       =  4
SILENT       =  8
DE_NOVO      = 16
PARENT_HOMOZYGOTE = 32

###############################
def get_variants_per_chrom(fpath):
	vars_per_chrom = {}
	with open(fpath, "r") as inf:
		for line in inf:
			if line[0] == "#": continue
			fields = line.strip().split()
			if len(fields) <3: continue
			chrom = fields[0]
			pos = int(fields[1])
			if not chrom in vars_per_chrom: vars_per_chrom[chrom] = {}
			vars_per_chrom[chrom][pos] = fields[2:]
	return vars_per_chrom


def variants_output(fnm, hdr, vars_per_chrom):
	genes_affected = set([])
	outf = open(fnm, "w")
	outf.write(hdr)
	# output the variants sorted
	for chrom, vars in vars_per_chrom.items():
		for pos in sorted(vars.keys()):
			outf.write("\t".join([str(x) for x in [chrom, pos] + vars[pos]]) + "\n")
			if ":" in vars[pos][2]:
				gene_symbol = vars[pos][2].split(":")[0]
				genes_affected.add(gene_symbol)
	outf.close()
	return genes_affected


def output_affected_genes(cursor, fnm, genes_affected, hgnc2omim):
	# TODO - hgnc2omim to indentifier_maps mysql
	outf = open(fnm,"w")
	for approved in genes_affected:
		if approved=="anon": continue
		uniprot = hgnc2uniprot(cursor, approved)
		[kegg_geneid, kegg_pathways] = uniprot2kegg(cursor, uniprot)
		if not kegg_pathways: kegg_pathways="x"
		omim = hgnc2omim.get(approved, "x")
		outf.write("\t".join([approved, uniprot, omim, str(kegg_geneid), kegg_pathways]) + "\n")
		# the "x" pathwys might be disease pathways which I chose not to use
		# if kegg_pathways and kegg_pathways != "x":
		# 	kegg_pathway_ids.update(set(kegg_pathways.split(";")))
	outf.close()

