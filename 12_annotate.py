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

import os

from fe_utils.mysql import *
from fe_utils.ensembl import *

# there could be multiple genes separated by ';'
# esi - one of exonic, splice, intronic, upstream, downstream
# e can be followed by another: and aaFromPOSaaTo
# chrom pos  ref alt1|alt2  mom1|mom2  pop1|pop2  gene:[esiud] freq

# TODO: reporting pseudogenes? (possible misalignment of fragments)


stable2approved = {}


####################################
def get_seq_region_ids(cursor):
	# coord_system with id 4 is chromosome for GRCh38
	qry = "select name, seq_region_id from seq_region where coord_system_id=4 "
	qry += "and length(name)<3"
	seq_region_id = dict(hard_landing_search(cursor,qry))
	return seq_region_id


######
def find_approved(cursor, stable_id):
	ret = error_intolerant_search(cursor, "select approved_symbol from icgc.hgnc where  ensembl_gene_id='%s'"%stable_id)
	stable2approved[stable_id] = (ret[0][0] if ret else "anon")


######
def annotate(cursor, seq_region_id, line):
	if line[0]=='#': return ""
	[chrom,  variant_pos,  gt,  gt_mom,  gt_dad] = line.strip().split("\t")
	variant_pos = int(variant_pos)
	# print(chrom,  pos,  gt,  gt_mom,  gt_dad)
	qry  = "select gene_id, stable_id, seq_region_start, seq_region_end from gene where seq_region_id=%d "%seq_region_id[chrom]
	qry += "and seq_region_start-100<=%d and %d<=seq_region_end+100" % (variant_pos, variant_pos)
	ret  = error_intolerant_search(cursor, qry)
	if not ret:
		print("no return for", qry)
		return ""

	print("============     %d" % variant_pos)
	gene_annotation_fields = []
	for [gene_id, stable_id, seq_region_start, seq_region_end] in ret:
		if stable_id not in stable2approved: # see if we have it stored by any chance
			find_approved(cursor, stable_id)
		print(chrom,  variant_pos, gene_id, stable_id, stable2approved[stable_id])
		# find all canonical exons in this gene
		qry = "select start_in_gene, end_in_gene from gene2exon  "
		qry += "where gene_id=%d and is_canonical=1" % gene_id
		ret2 = error_intolerant_search(cursor, qry)
		if  not ret2:
			print("no exons found")
			continue
		variant_pos_in_gene = variant_pos - seq_region_start
		exon_intervals =  sorted(ret2, key=lambda x: x[0])


		placed = False
		location_code = "?"
		if variant_pos_in_gene < exon_intervals[0][0]:
			# print("   %d  <---- upstream" % variant_pos_in_gene)
			placed = True
			location_code = "u"
		exon_id = 0
		while not placed and exon_id<len(exon_intervals):
			[start, end] = exon_intervals[exon_id]
			#print(start, end)
			if not placed and variant_pos_in_gene < start:
				#print("   %d  <---- intronic" % variant_pos_in_gene)
				location_code = "i"
				placed = True
			elif not placed and variant_pos_in_gene <= end:
				#print("   %d  <---- exonic" % variant_pos_in_gene)
				location_code = "e"
				placed = True
			exon_id += 1
		if not placed and variant_pos_in_gene < exon_intervals[-1][1]:
			#print("   %d  <---- downstream" % variant_pos_in_gene)
			location_code = "d"
		gene_annotation_fields.append("%s:%s"%(stable2approved[stable_id], location_code))
	print(";".join(gene_annotation_fields))


############################
def main():
	fnm = "test1.hg38.vcf"

	expected_assembly = "hg38"
	if not expected_assembly in fnm:
		print("The expected assembly %s does not appear in the name '%s'." % (expected_assembly, fnm))
		print("Note that, as of this writing, it is the version used by Ensembl.")
		exit()
	vcfdir = "/storage/sequencing/openhumans/fakexomes/progeny"
	fpath = "{}/{}".format(vcfdir, fnm)
	mysql_conf_file = "/home/ivana/.tcga_conf"
	for dep in [vcfdir, fpath, mysql_conf_file]:
		if not os.path.exists(dep):
			print(dep, "not found")
			exit()

	db = connect_to_mysql(mysql_conf_file)
	cursor = db.cursor()
	[all_species, ensembl_db_name] = get_species(cursor)
	switch_to_db(cursor, ensembl_db_name['homo_sapiens'])
	# coord_system with id 4 is chromosome for GRCh38
	seq_region_id = get_seq_region_ids(cursor)

	inf = open(fpath,"r")
	for line in inf:
		line_annotated = annotate(cursor, seq_region_id, line)
		# print(line_annotated)
	inf.close()

	cursor.close()
	db.close()
	return


#########################################
if __name__ == '__main__':
	main()
