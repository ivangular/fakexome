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
#
import subprocess, os, sys
import urllib.request, json


############################
def parseOpenHumansInfo(data, downloaddir):

	# note that many files are not really available,
	# even though they are listed as so
	# in such case the server will return

	for k in ['count', 'next', 'previous']:
		print(data[k])

	count = 0
	for k in data['results']:
		if not 'metadata' in k or not 'tags' in k['metadata']: continue
		# gvcf is compressed vcf
		if not ("vcf" in k['metadata']['tags'] or 'gvcf' in k['metadata']['tags']): continue
		if 'user_notes' in k['metadata']['tags']: continue
		filename = k['basename']
		if 'imputed' in filename: continue
		if not (filename[-3:]=='vcf' or filename[-6:]=='vcf.gz' or filename[-7:]=='vcf.bz2'): continue
		# files with the names PGP seem the most promising
		# looks like it it the Hrvarde Pesonal Genome Project https://pgp.med.harvard.edu/
		if not filename[:3]=="PGP": continue
		count += 1
		print(count)
		print(filename)
		print(k['metadata'])
		url = k['download_url']
		downloadfile = "{}/{}".format(downloaddir, filename)
		zipped = False
		if downloadfile[-4:]==".bz2":
			unzipped = downloadfile[0:-4]
			zipped = "bzip"
		elif downloadfile[-3:]==".gz":
			unzipped = downloadfile[0:-3]
			zipped = "gzip"
		else:
			unzipped = downloadfile

		if os.path.exists(unzipped) and os.path.getsize(unzipped)>0:
			print(unzipped, "found\n")
			continue

		if not os.path.exists(downloadfile) or os.path.getsize(downloadfile)==0:
			cmd = "curl -L '{}' -o {} ".format(url, downloadfile)
			subprocess.call(cmd, shell=True)
		else:
			print("found", downloadfile)

		if zipped:
			print("unzipping")
			if zipped == "bzip":
				subprocess.call("bzip2 -d %s"%downloadfile, shell=True)
				# if unzipped not created, remove - it is probably one of the unavailable files
			elif zipped == "gzip":
				subprocess.call("gunzip %s"%downloadfile, shell=True)
			if not os.path.exists(unzipped):
				print("unzip failed, removing", downloadfile)
				os.remove(downloadfile)

		print()

############################
def main():

	# check out what's available on openhumans
	# currently there ar 822 entries of type 11 (human genome data)
	# checkout count firld in the downloaded json file
	# if you as for too many files the server will throttle your request
	# if needed write some more intelligent download file than this
	url = "https://www.openhumans.org/api/public/datatype/11/datafiles/?format=json&limit=900&offset=0"
	jsonfile = "openhumans.json"
	downloaddir = "/storage/sequencing/openhumans/orig"

	# this way it takes a bit too long to do repatedly each time
	# with urllib.request.urlopen(url) as dwld:
	# 	data = json.loads(dwld.read().decode())
	# 	print(data)
	if not os.path.exists(jsonfile) or os.path.getsize(jsonfile)==0:
		cmd = "curl -L '{}' -o {} ".format(url, jsonfile)
		subprocess.call(cmd, shell=True)
	with open(jsonfile,"r") as inf:
		data = json.loads(inf.read())
		parseOpenHumansInfo(data, downloaddir)
	return


#########################################
if __name__ == '__main__':
	main()
