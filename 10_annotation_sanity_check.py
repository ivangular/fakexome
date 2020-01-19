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

from random import random

from fe_utils.mysql import *
from fe_utils.ensembl import *
stable2approved = {}


####################################
def get_seq_region_ids(cursor):
	# coord_system with id 4 is chromosome for GRCh38
	qry  = "select name, seq_region_id from seq_region where coord_system_id=4 "
	qry += "and length(name)<3"
	seq_region_id = dict(hard_landing_search(cursor,qry))
	return seq_region_id


def find_approved(cursor, stable_id):
	ret = error_intolerant_search(cursor, "select approved_symbol from icgc.hgnc where  ensembl_gene_id='%s'"%stable_id)
	stable2approved[stable_id] = (ret[0][0] if ret else "anon")


def annotation_check(cursor, seq_region_id, line):
	if line[:3] != 'chr': return ""
	[chrom,  start, end, annot] = line.strip().split()
	chr_short = chrom.replace("chr","")
	diagnostics = "##########################\n"
	expected_symbol = annot.split(",")[0].split("|")[1] if "," in annot else annot
	diagnostics += expected_symbol + "\n"
	expected = {"anon", expected_symbol} # set literal
	discrepancy = False
	for pos in [start,end]:
		qry  = "select gene_id, stable_id from gene where seq_region_id=%d "%seq_region_id[chr_short]
		qry += "and seq_region_start-100<=%s and %s<=seq_region_end+100" % (pos, pos)
		ret  = error_intolerant_search(cursor, qry)
		if not ret:
			diagnostics += "no return for %s\n" % qry
			discrepancy = True
		else:
			approved = {"-"}
			for [gene_id, stable_id] in ret:
				if not stable_id in stable2approved:
					find_approved(cursor, stable_id)
				diagnostics += " %s %s  %s  %s   %s \n" % (chrom,  pos, str(gene_id), stable_id, stable2approved[stable_id])
				approved.add(stable2approved[stable_id])
			if not expected & approved:
				discrepancy = True
	if discrepancy:
		print(diagnostics)
		return False

	return True

############################
def main():
	# this file is already annotated
	fnm = "agilent_target_regions.S07084713.hg38.bed"

	expected_assembly = "hg38"
	if not expected_assembly in fnm:
		print("The expected assembly %s does not appear in the name '%s'." % (expected_assembly, fnm))
		print("Note that, as of this writing, it is the version used by Ensembl.")
		exit()
	vcfdir = "/storage/sequencing/openhumans/"
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
	ok = 0
	fail = 0
	for line in inf:
		#if random() < 0.01:
		success = annotation_check(cursor, seq_region_id, line)
		if success:
			ok += 1
		else:
			fail += 1

	print(ok, fail)

	inf.close()
	cursor.close()
	db.close()
	return


#########################################
if __name__ == '__main__':
	main()
