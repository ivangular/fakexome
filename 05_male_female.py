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

# divide males and females so we don't have to do it all over again each time

import os, subprocess


############################
def main():
	vcfdir  = "/storage/sequencing/openhumans/fakexomes"
	for dep in [vcfdir]:
		if not os.path.exists(dep):
			print(dep, "not found")
			exit()
	maledir   = "{}/males".format(vcfdir)
	femaledir = "{}/females".format(vcfdir)
	for d in [maledir, femaledir]:
		if not os.path.exists(d): os.mkdir(d)

	for fnm in os.listdir(vcfdir):
		if fnm[-4:]!='.vcf': continue
		cmd = "grep chrY %s/%s | wc -l" % (vcfdir, fnm)
		retval = subprocess.run([cmd], shell=True, stdout=subprocess.PIPE)
		y_count = int(retval.stdout.decode('utf-8').strip())
		print(fnm, y_count)
		if y_count==0:
			os.rename("{}/{}".format(vcfdir, fnm), "{}/{}".format(femaledir, fnm))
		else:
			os.rename("{}/{}".format(vcfdir, fnm), "{}/{}".format(maledir, fnm))

	return


#########################################
if __name__ == '__main__':
	main()
