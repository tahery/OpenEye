#!/usr/bin/env python3
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
# Depicts active site maps
#############################################################################

import sys
from openeye import oechem
from openeye import oedepict
from openeye import oegrapheme
import argparse

def main(argv=[__name__]):

    itf = oechem.OEInterface()
    oechem.OEConfigure(itf, InterfaceData)
    oedepict.OEConfigureImageWidth(itf, 900.0)
    oedepict.OEConfigureImageHeight(itf, 600.0)
    oechem.OEConfigureSplitMolComplexOptions(itf, oechem.OESplitMolComplexSetup_LigName |
                                             oechem.OESplitMolComplexSetup_CovLig)

    #if not oechem.OEParseCommandLine(itf, argv):
    #    return 1
    
    argParser = argparse.ArgumentParser()
    argParser.add_argument("--complex", help="Complex Input", dest="iname", required=True)
    argParser.add_argument("--out", help="Output in ONLY SVG", dest="oname", default="interactions.svg")
    argParser.add_argument("--pname", help="Protein Input", dest="pname", default="")
    argParser.add_argument("--lname", help="Ligand Input", dest="lname", default="")
    argParser.add_argument("--hbond", help="Maximum Distance of Hydrogen Bond (default=3.2A)", type=float, default=3.2)
    argParser.add_argument("--hbondni", help="Maximum Distance of Non-Ideal Hydrogen Bond (default=3.7A)", type=float, default=3.8)
    argParser.add_argument("--hbondca", help="Maximum Distance of Charge-Aided Hydrogen Bond (default=3.5A)", type=float, default=3.5)
    argParser.add_argument("--halogenbond", help="Maximum Distance of Halogen Bond (default=3.2A)", type=float, default=3.2)
    argParser.add_argument("--pistack", help="Maximum Distance of Pi-Stacking (default=5.0A)", type=float, default=5.0)
    argParser.add_argument("--tstack", help="Maximum Distance of T-Stacking (default=5.35A)", type=float, default=5.35)
    argParser.add_argument("--saltbridge", help="Maximum Distance of Salt Bridge (default=5.0A)", type=float, default=5.0)
    argParser.add_argument("--cationpi", help="Maximum Distance of Cation-Pi (default=5.5A)", type=float, default=5.5)
    argParser.add_argument("--contact", help="Maximum Distance of Contacts (default=1.2A)", type=float, default=1.2)
    argParser.add_argument("--clashcontact", help="Minimum Distance of Clash-Contacts (default=0.8A)", type=float, default=0.8)
    args = argParser.parse_args()

    oname = args.oname

    ext = oechem.OEGetFileExtension(oname)
    if ext != "svg":
        oechem.OEThrow.Fatal("Only available for SVG image type!")

    ofs = oechem.oeofstream()
    if not ofs.open(oname):
        oechem.OEThrow.Fatal("Cannot open output file!")

    # initialize protein and ligand

    protein = oechem.OEGraphMol()
    ligand = oechem.OEGraphMol()
    if not get_protein_and_ligands(protein, ligand, itf, args):
        oechem.OEThrow.Fatal("Cannot initialize protein and/or ligand!")

    # depict active site maps

    width, height = oedepict.OEGetImageWidth(itf), oedepict.OEGetImageHeight(itf)
    image = oedepict.OEImage(width, height)

    depict_activesite_maps(image, protein, ligand, args)

    oedepict.OEDrawCurvedBorder(image, oedepict.OELightGreyPen, 10.0)

    oedepict.OEWriteImage(oname, image)

    return 0


def perceive_interation_hints_user_def_params(protein, ligand, args):
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
    opts.SetMaxNonIdealHBondDistance(args.hbondni) # Default : 3.8A
    opts.SetMaxChargeAidedHBondDistance(args.hbondca) # Default : 3.5A
    opts.SetMaxContactFraction(args.contact) # Default : 1.2A
    opts.SetMaxHalogenBondDistance(args.halogenbond) # Default : 3.2A   
    opts.SetMaxPiStackDistance(args.pistack) # Default : 5.0
    opts.SetMaxTStackDistance(args.tstack) # Default : 5.35
    opts.SetMaxSaltBridgeDistance(args.saltbridge) # Default : 5.0A
    opts.SetMaxCationPiDistance(args.cationpi) # Default : 5.5A
    opts.SetMinContactFraction(args.clashcontact) # Default : 0.8A
    print("Hydrogen Bond", opts.GetMaxHBondDistance()) # Default : 3.2A
    print("Non-Ideal Hydrogen Bond", opts.GetMaxNonIdealHBondDistance()) # Default : 3.7A
    print("Max Charge Aided Hydrogen Bond", opts.GetMaxChargeAidedHBondDistance()) # Default : 3.5A
    print("Min Contact", opts.GetMinContactFraction()) # Default : 0.8A
    print("Max Contact", opts.GetMaxContactFraction()) # Default : 1.2A
    print("Max Halogen Bond", opts.GetMaxHalogenBondDistance()) # Default : 3.2A   
    print("Max Pi Stack", opts.GetMaxPiStackDistance()) # Default : 5.0
    print("Max T Stack", opts.GetMaxTStackDistance()) # Default : 5.35
    print("Max Salt Bridge", opts.GetMaxSaltBridgeDistance()) # Default : 5.0A
    print("Max Cation Pi", opts.GetMaxCationPiDistance()) # Default : 5.5A
    oechem.OEPerceiveInteractionHints(asite, opts)

    return asite


def depict_activesite_maps(image, protein, ligand, args):
    """
    :type image: oedepict.OEImageBase
    :type protein: oechem.OEMolBase
    :type ligand: oechem.OEMolBase
    """

    # perceive interactions

    #asite = oechem.OEInteractionHintContainer(protein, ligand)
    asite = perceive_interation_hints_user_def_params(protein,ligand,args)
    if not asite.IsValid():
        oechem.OEThrow.Fatal("Cannot initialize active site!")
    asite.SetTitle(ligand.GetTitle())

    #oechem.OEPerceiveInteractionHints(asite)

    # depiction

    oegrapheme.OEPrepareActiveSiteDepiction(asite)
    oegrapheme.OERenderActiveSiteMaps(image, asite)


def split_complex(protein, ligand, sopts, complexmol):

    water = oechem.OEGraphMol()
    other = oechem.OEGraphMol()

    pfilter = sopts.GetProteinFilter()
    wfilter = sopts.GetWaterFilter()
    sopts.SetProteinFilter(oechem.OEOrRoleSet(pfilter, wfilter))
    filter = oechem.OEMolComplexFilterCategory_Nothing
    sopts.SetWaterFilter(oechem.OEMolComplexFilterFactory(filter))

    oechem.OESplitMolComplex(ligand, protein, water, other, complexmol, sopts)

    return ligand.NumAtoms() != 0 and protein.NumAtoms() != 0


def get_protein_and_ligands(protein, ligand, itf, args):

    if (args.iname):

        # separate ligand and protein in complex

        iname = args.iname

        ifs = oechem.oemolistream()
        if not ifs.open(iname):
            oechem.OEThrow.Fatal("Cannot open input complex file!")

        complexmol = oechem.OEGraphMol()
        if not oechem.OEReadMolecule(ifs, complexmol):
            oechem.OEThrow.Fatal("Unable to read complex from %s" % iname)

        if not oechem.OEHasResidues(complexmol):
            oechem.OEPerceiveResidues(complexmol, oechem.OEPreserveResInfo_All)

        sopts = oechem.OESplitMolComplexOptions()
        oechem.OESetupSplitMolComplexOptions(sopts, itf)

        if not split_complex(protein, ligand, sopts, complexmol):
            oechem.OEThrow.Fatal("Cannot separate complex!")
    else:

        # read ligand and protein from separate files

        pname = args.pname

        ifs = oechem.oemolistream()
        if not ifs.open(pname):
            oechem.OEThrow.Fatal("Cannot open input protein file!")

        if not oechem.OEReadMolecule(ifs, protein):
            oechem.OEThrow.Fatal("Unable to read protein from %s" % pname)

        lname = args.lname

        ifs = oechem.oemolistream()
        if not ifs.open(lname):
            oechem.OEThrow.Fatal("Cannot open input ligand file!")

        if not oechem.OEReadMolecule(ifs, ligand):
            oechem.OEThrow.Fatal("Unable to read ligand from %s" % lname)

    return ligand.NumAtoms() != 0 and protein.NumAtoms() != 0


#############################################################################
# INTERFACE
#############################################################################

InterfaceData = '''
!CATEGORY "input/output options :" 1

  !PARAMETER -complex 1
    !ALIAS -c
    !TYPE string
    !REQUIRED false
    !VISIBILITY simple
    !BRIEF Input filename of the protein-ligand complex
  !END

  !PARAMETER -protein 2
    !ALIAS -p
    !TYPE string
    !REQUIRED false
    !VISIBILITY simple
    !BRIEF Input filename of the protein
  !END

  !PARAMETER -ligand 3
    !ALIAS -l
    !TYPE string
    !REQUIRED false
    !VISIBILITY simple
    !BRIEF Input filename of the ligand
  !END

  !PARAMETER -out 4
    !ALIAS -o
    !TYPE string
    !REQUIRED true
    !VISIBILITY simple
    !BRIEF Output filename of the generated image
  !END

!END
'''

if __name__ == "__main__":
    sys.exit(main(sys.argv))
