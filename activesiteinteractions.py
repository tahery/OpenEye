#!/usr/bin/env python
# (C) 2017 OpenEye Scientific Software Inc. All rights reserved.
#
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of OpenEye products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms. OpenEye claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable OpenEye offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall OpenEye be
# liable for any damages or liability in connection with the Sample Code
# or its use.

#############################################################################
# Derived script of 'printinteractions.py' of openeye
#############################################################################
import sys
from openeye import oechem
import argparse
import csv
import json

def main(argv=[__name__]):

    itf = oechem.OEInterface()
    oechem.OEConfigure(itf, InterfaceData)
    oechem.OEConfigureSplitMolComplexOptions(itf, oechem.OESplitMolComplexSetup_LigName)

    #if not oechem.OEParseCommandLine(itf, argv):
    #    return 0
    
    #iname = itf.GetString("-complex")

    argParser = argparse.ArgumentParser()
    argParser.add_argument("--complex", help="PDB Complex", dest="iname")
    argParser.add_argument("--hbond", help="Maximum Distance of Hydrogen Bond (default=3.2A)", type=float, default=3.2)
    argParser.add_argument("--hbondni", help="Maximum Distance of Non-Ideal Hydrogen Bond (default=3.7A)", type=float, default=3.8)
    argParser.add_argument("--hbondca", help="Maximum Distance of Charge-Aided Hydrogen Bond (default=3.5A)", type=float, default=3.5)
    argParser.add_argument("--halogenbond", help="Maximum Distance of Halogen Bond (default=3.2A)", type=float, default=3.2)
    argParser.add_argument("--pistack", help="Maximum Distance of Pi-Stacking (default=5.0A)", type=float, default=5.0)
    argParser.add_argument("--tstack", help="Maximum Distance of T-Stacking (default=5.35A)", type=float, default=5.35)
    argParser.add_argument("--saltbridge", help="Maximum Distance of Salt Bridge (default=5.0A)", type=float, default=5.0)
    argParser.add_argument("--cationpi", help="Maximum Distance of Cation-Pi (default=5.5A)", type=float, default=5.5)
    argParser.add_argument("--contact", help="Maximum Distance of Contacts (default=1.2A)", type=float, default=1.2)
    argParser.add_argument("--clashcontact", help="Minimum Distance of Contacts (default=0.8A)", type=float, default=0.8)
    argParser.add_argument("--subsites", help="1 (True) or 0 (False) (default=1)", type=int, choices=[0,1], default=1)
    args = argParser.parse_args()


    iname = args.iname
    ifs = oechem.oemolistream()
    if not ifs.open(iname):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % iname)

    complexmol = oechem.OEGraphMol()
    if not oechem.OEReadMolecule(ifs, complexmol):
        oechem.OEThrow.Fatal("Unable to read molecule from %s" % iname)

    if not oechem.OEHasResidues(complexmol):
        oechem.OEPerceiveResidues(complexmol, oechem.OEPreserveResInfo_All)

    # Separate ligand and protein

    sopts = oechem.OESplitMolComplexOptions()
    oechem.OESetupSplitMolComplexOptions(sopts, itf)

    ligand = oechem.OEGraphMol()
    protein = oechem.OEGraphMol()
    water = oechem.OEGraphMol()
    other = oechem.OEGraphMol()

    pfilter = sopts.GetProteinFilter()
    wfilter = sopts.GetWaterFilter()
    sopts.SetProteinFilter(oechem.OEOrRoleSet(pfilter, wfilter))
    sopts.SetWaterFilter(
          oechem.OEMolComplexFilterFactory(oechem.OEMolComplexFilterCategory_Nothing))

    oechem.OESplitMolComplex(ligand, protein, water, other, complexmol, sopts)

    if ligand.NumAtoms() == 0:
        oechem.OEThrow.Fatal("Cannot separate complex!")

    # Perceive interactions

    asite = perceive_interaction_hints_user_def_params(protein,ligand,args)
   
    list_type_contacts= ["hbond", "halogen","stacking","sbridge","cation-pi","clashcontact","contact"]
    get_protein_interactions(asite, list_type_contacts,args.subsites) 

def perceive_interaction_hints_user_def_params(protein, ligand, args):
    """
    :type protein: oechem.OEMol
    :type ligand: oechem.OEMol
    :rtype: oechem.OEInteractionHintContainer
    """
    asite = oechem.OEInteractionHintContainer()
    asite.AddMolecule(protein, oechem.OEProteinInteractionHintComponent())
    asite.AddMolecule(ligand, oechem.OELigandInteractionHintComponent())
    if not oechem.OEIsValidActiveSite(asite):
        oechem.OEThrow.Fatal("Cannot initialize active site!")

    opts = oechem.OEPerceiveInteractionOptions()
    opts.SetMaxHBondDistance(args.hbond) # Default : 3.2A
    opts.SetMaxNonIdealHBondDistance(args.hbondni) # Default : 3.7A
    opts.SetMaxChargeAidedHBondDistance(args.hbondca) # Default : 3.5A
    opts.SetMaxContactFraction(args.contact) # Default : 1.2A
    opts.SetMaxHalogenBondDistance(args.halogenbond) # Default : 3.2A   
    opts.SetMaxPiStackDistance(args.pistack) # Default : 5.0
    opts.SetMaxTStackDistance(args.tstack) # Default : 5.35
    opts.SetMaxSaltBridgeDistance(args.saltbridge) # Default : 5.0A
    opts.SetMaxCationPiDistance(args.cationpi) # Default : 5.5A
    opts.SetMinContactFraction(args.clashcontact) # Default : 0.8A
    oechem.OEPerceiveInteractionHints(asite, opts)

    return asite



def GetResidueName(residue):
    return {"ResName_Prot": residue.GetName(), "ResID_Prot": residue.GetResidueNumber(), "ChainID_Prot": residue.GetChainID()}

def GetSubsite(subsite, resid):
    if subsite == 1: 
        subsites = {"SS01":[30,71,108,115,118], "SS02":[32,228], "SS03":[72,73,107], "SS04":[230,231], "SS05":[34,198], "SS06":[110], "SS07":[232,233,235,325], "SS08":[35,69,70,76,126,128], "SS09":[11,13,14,229,335], "SS10":[224,226,329,332]}
        for sub in subsites:
            if resid in subsites[sub]:
                return sub
    return None

def GetResidBACE2(resid):
    '''if resid < 165:
        return resid + 16
    elif resid < 318:
        return resid + 13
    else: 
        return resid + 12''' # non-normalized 


def get_interactions(interactiontype):
    if interactiontype == "stacking":
        interactiontype = oechem.OEAndInteractionHint(oechem.OEIsStackingInteractionHint(),oechem.OEIsInterInteractionHint())
    elif interactiontype == "hbond":
        interactiontype = oechem.OEAndInteractionHint(oechem.OEIsHBondInteractionHint(),oechem.OEIsInterInteractionHint())
    elif interactiontype == "sbridge":
        interactiontype = oechem.OEAndInteractionHint(oechem.OEIsSaltBridgeInteractionHint(),oechem.OEIsInterInteractionHint())
    elif interactiontype == "halogen":
        interactiontype = oechem.OEAndInteractionHint(oechem.OEIsHalogenBondInteractionHint(),oechem.OEIsInterInteractionHint())
    elif interactiontype == "cation-pi":
        interactiontype = oechem.OEAndInteractionHint(oechem.OEIsCationPiInteractionHint(),oechem.OEIsInterInteractionHint())
    elif interactiontype == "clashcontact":
        interactiontype = oechem.OEAndInteractionHint(oechem.OEIsClashInteractionHint(),oechem.OEIsInterInteractionHint())
    elif interactiontype == "contact":
        interactiontype = oechem.OEAndInteractionHint(oechem.OEIsContactInteractionHint(),oechem.OEIsInterInteractionHint())
    return interactiontype



def get_protein_interactions(asite,list_type_contacts,subsites):
    prot = asite.GetMolecule(oechem.OEProteinInteractionHintComponent())
    fields = ["ResName_Prot","ResID_Prot","ChainID_Prot","AtomName_Prot","AtomsName_Lig","ContactType","Subsite"]
    listinteractions = []

    ftsv = open("interactions.tsv", "w")
    fjson = open("interactions.json", "w")
    
    writer_tsv = csv.DictWriter(ftsv, fields, delimiter='\t')
    writer_tsv.writeheader()
    
    for ContactType in list_type_contacts:
        atoms = list()
        for atom in prot.GetAtoms():
            if asite.HasInteraction(oechem.OEAndInteractionHint(oechem.OEHasInteractionHint(atom), get_interactions(ContactType))):
                atoms.append(atom)
        
        for i in atoms:
            ligatom = set()
            residue = oechem.OEAtomGetResidue(i)
            
            if "HOH" not in residue.GetName() and "TIP" not in residue.GetName():
                for inter in asite.GetInteractions(oechem.OEHasInteractionHint(i)):
                    ligfrag = inter.GetFragment(oechem.OELigandInteractionHintComponent())
                    if ligfrag is None:
                        continue
                    for latom in ligfrag.GetAtoms():
                        ligatom.add(str(latom).strip().replace(" ",""))
                interaction = GetResidueName(residue)
                interaction["AtomName_Prot"] = i.GetName()
                interaction["AtomsName_Lig"] = list(ligatom)
                interaction["ContactType"] = ContactType
                if subsites == 1: 
                    subsite = GetSubsite(subsites, interaction["ResID_Prot"])
                    interaction["Subsite"] = subsite
                else:
                    interaction["Subsite"] = None
                listinteractions.append(interaction)
    len_listinteractions = len(listinteractions)
    i = 1
    print("[", file=fjson)
    for dic in listinteractions:
        writer_tsv.writerow(dic)
        json_object = json.dumps(dic, indent = 4)
        if i < len_listinteractions: 
            print(json_object+",", file=fjson)
        else:
            print(json_object+"]", file=fjson)
        i += 1
    ftsv.close()
    fjson.close()
    return None

InterfaceData = '''
!BRIEF printinteractions [-complex] <input>

!CATEGORY "input/output options :"

  !PARAMETER -complex
    !ALIAS -c
    !TYPE string
    !KEYLESS 1
    !REQUIRED true
    !VISIBILITY simple
    !BRIEF Input filename of the protein complex
  !END

!END
'''

if __name__ == "__main__":
    sys.exit(main(sys.argv))
