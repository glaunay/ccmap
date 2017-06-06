import sys, os
sys.path.append("/Users/guillaumelaunay/work/DVL/pyproteinsExt/src")
sys.path.append("/Users/guillaumelaunay/work/DVL/pyproteins/src")
import pyproteinsExt.structure.coordinates as PDB
import pyproteinsExt.structure.operations as PDBop

parser = PDB.Parser()
pdbObj = parser.load(file="/Users/guillaumelaunay/work/mesh/pdbFiles/2FJU.pdb")
pdbObjRec = pdbObj.chain("A")
pdbObjLig = pdbObj.chain("X")

import ccmap
atomDatumRec = pdbObjRec.atomDictorize
atomDatumLig = pdbObjLig.atomDictorize
jsonResults = ccmap.duals([(atomDatumRec, atomDatumLig)], 4.5)


print jsonResults
