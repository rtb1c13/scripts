#!/usr/bin/env python

# Script using Biopython objects to import and adapt PDB files
from Bio.PDB import *
import numpy as np
import requests
import argparse
parser=PDBParser(QUIET=True)

#Protocol
# 1) Read in arguments for downloading PDB files
# 2) Read in table of PDB entries and CDR ranges
# 3) Delete atoms in CDR ranges
# 4) Delete Hetatms
# 5) Print out CA atoms only
def main():
  ftp=False
  options=argparse.ArgumentParser(description='Processes and edits PDB files ready for overlay and PCA.')
  options.add_argument('-f','--ftp',action='store_true',help="Uses ftp downloads of structures from RCSB PDB")
  args=options.parse_args()
  ftp=args.ftp
  table=rangeimport("../Abysis/FR_CDR_table.txt")
  heavyobs=[]
  for line in table:
    pdbid=line[0]
    if line[0].startswith('#'):
       break
    structure1=getpdb(ftp,pdbid)
    structure2=abnum(pdbid)
    cdrlist = []
    cdrlist = line[2].split('-')
    cdrlist.extend(line[4].split('-'))
    cdrlist.extend(line[6].split('-'))
    cdrlist=np.array(map(int, cdrlist))
    if structure1[0].has_id('H'):
      heavy1=structure1[0]['H']
    elif structure1[0].has_id('B'):
      heavy1=structure1[0]['B']
    else:
      print "PDB %s has no chain H or chain B. Skipping editing of heavy chain." % pdbid
      break
    try: heavy2=structure2[0]['H']
    except KeyError: print "Cannot find chain H in Kabat numbered structure of %s.\n This PDB should be edited separately!" % pdbid
# Delete Hetatm records
    for j in [ i for i in heavy1 if i.id[0] is not ' ' ]:
      heavy1.__delitem__(j.get_id())
    heavy1=renum_heavy(heavy1,heavy2)
    heavy1=delheavy(cdrlist,heavy1)
    heavyobs.append([pdbid,heavy1])
# Find consensus start point of heavys
  
  heavyobs=chopstart(heavyobs)
  for i in heavyobs:
    io=PDBIO()
    io.set_structure(i[1])
    io.save("Created/%s_H_edited.pdb" % i[0], CASelect())
# Delete residues, need to work out how to check for chain H or B etc.    
    

def getpdb(ftp,pdbid):
  if ftp is True:
    pdbl=PDBList()
    structure1=parser.get_structure(pdbid,pdbl.retrieve_pdb_file(pdbid,pdir="Downloaded"))
    return structure1
  else:
    pdbfile = pdbid + ".pdb"
    structure1=parser.get_structure(pdbid,pdbfile)
    return structure1

def abnum(pdbid):
  oldfn = "pdb" + pdbid.lower() + ".ent"
  newfn = "Created/" + pdbid + "_V_Kabat.pdb"
  files={'pdb': open('./Downloaded/%s' % oldfn)}
  data={'scheme': '-k'}
  url='http://www.bioinf.org.uk/cgi-bin/abnum/abnumpdb.pl'
  r = requests.post(url,files=files,data=data)
  f = open(newfn, 'w')
  f.write(r.text)
  f.close()
  structure2=parser.get_structure(pdbid + "2",newfn)
  return structure2

def renum_heavy(heavy1,heavy2):
# Function to renumber chain using kabat numbered chain
# Must first renumber original chain to avoid duplicate resids
# Then insert Kabat numbered chain
# Then renumber residues sequentially following Kabat numbering

# Step 1, renumber heavy1 chain. Must renumber child_list and child_dict
  count = 5001
#  print len(heavy1.child_dict)
  for residue in heavy1:
    origid = residue.id
#    print residue.id
    newid = (' ',count,' ')
    residue.id = newid
    heavy1.child_dict[newid] = heavy1.child_dict.pop(origid)
#    print heavy1.child_dict[newid]
    count += 1
# Step 2, detach old numbers, insert Kabat numbers
  maxindex=len(heavy2)
  count = 0
  resids=[]
#  print len(heavy1.child_dict)
  for resid in heavy2.get_list():
    resids.append(resid.id)
  while count < len(resids):
    heavy1.detach_child(heavy1.get_list()[0].id)
    count += 1
  index = 0
  for residue in heavy2:
    heavy1.insert(index,residue)
    index += 1
# Step 3, renumber res following Kabat numbers
  count = int(heavy1.get_list()[index-1].id[1]) + 1
  while index < len(heavy1):
    origid = heavy1.get_list()[index].id
    newid = (' ',count,' ')
    heavy1.get_list()[index].id = newid
    heavy1.child_dict[newid] = heavy1.child_dict.pop(origid)
    index += 1
    count += 1 
  return heavy1

def rangeimport(filename):
#Import table of CDR regions to delete
  table=np.genfromtxt(filename, skip_header=9, dtype=None)
  print "Now importing CDR ranges from file %s \n" % filename
  print "(The first entry looks like this):"
  print str(table[0])
  return table

def delheavy(cdrlist,heavy1):
#Delete CDR regions in heavy chain
  for delres in range(cdrlist[0],cdrlist[1]+1,1):
    if heavy1.has_id(delres):
      for j in [ i for i in heavy1 if i.id[1] == delres ]:
        heavy1.__delitem__(j.get_id())
  for delres in range(cdrlist[2],cdrlist[3]+1,1):
    if heavy1.has_id(delres):
      for j in [ i for i in heavy1 if i.id[1] == delres ]:
        heavy1.__delitem__(j.get_id())
  for delres in range(cdrlist[4],cdrlist[5]+1,1):
    if heavy1.has_id(delres):
      for j in [ i for i in heavy1 if i.id[1] == delres ]:
        heavy1.__delitem__(j.get_id())
#if heavy1.has_id(82):
#  for j in [ i for i in heavy1 if i.id[1] == 82 ]:
#    heavy1.__delitem__(j.get_id())
# Delete Hetatm records
#  for j in [ i for i in heavy1 if i.id[0] is not ' ' ]:
#    heavy1.__delitem__(j.get_id())
  return heavy1

def chopstart(heavyobs):
#Finds the chain with the highest number start point, chops all others to this length
  resid_list=[]
  min_n=int(1)
  for obj in heavyobs:
    resids=[]
    for resid in obj[1].get_list():
      resids.append(resid.id[1])
    resid_list.append(resids)
  for j in resid_list:
    if np.min(j) > min_n:
      min_n = np.min(j)
  print "Now chopping up to residue %i " % min_n
  for obj in heavyobs:
    for delres in range(1,min_n,1):
      if obj[1].has_id(delres):
        for j in [ i for i in obj[1] if i.id[1] == delres ]:
          obj[1].__delitem__(j.get_id())
# Now finds lengths of chains, chops to shortest one
  minlen = int(100000)
  maxres = int(1)
  for j in resid_list:
    if np.max(j) > maxres:
      maxres = np.max(j)
  for obj in heavyobs:
    if minlen > len(obj[1]):
      minlen = len(obj[1])
#      minres = obj[1][minlen].id[1]
  print "Now chopping all chains to be %i long " % minlen
  for obj in heavyobs:
    for index in range(len(obj[1]),minlen-1,-1):
      try: 
#        print "Now deleting %s residue %s" % (obj[0], obj[1].get_list()[index].id[1])
        obj[1].__delitem__(obj[1].get_list()[index].id[1])
      except IndexError:
#        print "Not editing ",obj[0]
        continue
  
  return heavyobs
  


# Only print out CA atoms. Need to define a Select class for this
class CASelect(Select):
  def accept_atom(self, atom):
    if atom.get_name()=='CA':
      if not atom.is_disordered() or atom.get_altloc()=='A':
        return 1
      else:
        return 0
    else:
      return 0

# Main runs here
main()
