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



from fe_utils.annotation import *


############################
def outer_loop(base_name):

	fnm = "{}.hg38.vcf".format(base_name)
	print("annotating", fnm)
	time0 = time()

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


	inf = open(fpath,"r")
	outf = open("{}.annotated.vcf".format(base_name),"w")
	line_ct = 0
	for line in inf:
		if line[0]=='#': continue
		annotation_fields = annotate(cursor,  line)
		line_annotated = add_flags(annotation_fields)
		#print(line_annotated, "       {0:b}".format(int(line_annotated.split("\t")[-1])))
		outf.write(line_annotated + "\n")
		line_ct += 1
	inf.close()
	outf.close()


	cursor.close()
	db.close()

	print(fnm, "done, %2.1f min" % ((time()-time0)/60) )
	
	return

############################
def main():
	for i in range(15):
		outer_loop("child%d"%i)
	return

#########################################
if __name__ == '__main__':
	main()
