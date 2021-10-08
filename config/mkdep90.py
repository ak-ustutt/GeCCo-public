#!/usr/bin/python3
#
# $Id: mkdep90.py 1.1 03/10/31 15:14:43-00:00 jonas@ $
#

import sys, string, re

use=re.compile('\s*use\s+(\w+)[^!,]*', re.I)
mod=re.compile('\s*module(?!\s+procedure)\s+(\w+)[^!,]*', re.I)
include=re.compile('\s*#?include\s+["<\']?([\w\._]+)[">\']?[^!,]*', re.I)

class depfile:
	def __init__(self, name):
		self.name=name
		self.oname=""
		self.uses={}
		self.includes={}
		ri=name.rindex('.')
		self.oname=name[:ri]+'.o'
	
	def adduses(self, u):
		self.uses[u]=''

	def addincludes(self, u):
		self.includes[u]=''

def main():
	allmods={}
	dfiles=[]
	for ff in sys.argv[1:]:
		fd=open(ff, 'r')
		buf=fd.readlines()
		fd.close()
		dfiles.append(depfile(ff))

		for ln in buf:
			m=mod.match(ln)
			if m is not None:
				allmods[m.group(1).lower()]=ff
				continue

			m=use.match(ln)
			if m is not None:
				dfiles[-1].adduses(m.group(1).lower())
				continue

			m=include.match(ln)
			if m is not None:
				dfiles[-1].addincludes(m.group(1).lower())
				continue

	for df in dfiles:
		deps=[]
		for dd in list(df.uses.keys()):
			try:
				if (allmods[dd] != df.name):
					ri=allmods[dd].rindex('.')
					omod='$(ARCH)/'+allmods[dd][:ri]+'.o'
					deps.append(omod)
			except:
				print('Missing dependency for', dd, 'in',\
				df.name, file=sys.stderr)

		for dd in list(df.includes.keys()):
				deps.append(dd)

		if deps:
			dstr='$(ARCH)/'+df.oname+': '
			for i in deps:
				dstr=dstr+i+' '
			print(dstr)


if __name__ == '__main__':
	main()
