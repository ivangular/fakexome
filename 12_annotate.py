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


####################################
def get_seq_region_ids(cursor):
	# coord_system with id 4 is chromosome for GRCh38
	qry = "select name, seq_region_id from seq_region where coord_system_id=4 "
	qry += "and length(name)<3"
	seq_region_id = dict(hard_landing_search(cursor,qry))
	return seq_region_id


stable2approved = {}


def find_approved(cursor, stable_id):
	ret = error_intolerant_search(cursor, "select approved_symbol from icgc.hgnc where  ensembl_gene_id='%s'"%stable_id)
	stable2approved[stable_id] = (ret[0][0] if ret else "anon")


def annotate(cursor, seq_region_id, line):
	if line[0]=='#': return ""
	[chrom,  pos,  gt,  gt_mom,  gt_dad] = line.strip().split("\t")
	pos = int(pos)
	#rint(chrom,  pos,  gt,  gt_mom,  gt_dad)
	qry  = "select gene_id, stable_id from gene where seq_region_id=%d "%seq_region_id[chrom]
	qry += "and seq_region_start-100<=%d and %d<=seq_region_end+100" % (pos, pos)
	ret  = error_intolerant_search(cursor, qry)
	if not ret:
		print("no return for", qry)
		return ""
	for [gene_id, stable_id] in ret:
		if not stable_id in stable2approved:
			find_approved(cursor, stable_id)
		print(chrom,  pos, gene_id, stable_id, stable2approved[stable_id])
	#xit()

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
		#print(line_annotated)
	inf.close()

	cursor.close()
	db.close()
	return


#########################################
if __name__ == '__main__':
	main()
