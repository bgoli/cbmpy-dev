"""

This algorithm selects candidate metabolites and creates individual sbml models with sink reactions for these candidates.


"""

import pyscescbm as cbm
import sys

byproducts = []

def getByProd(modelFile,coupledSubstrates, oldSubstrates, reactionList, exclusives):
    newSubstrates = []
    
    # Find product associated reactions
    for product in coupledSubstrates:
        follow = True
        for met in exclusives:
            if met in product:
                follow = False
                break
        if follow is True:
            reactions = modelFile.getFluxesAssociatedWithSpecies(product)
            for reaction in reactions:
                # Verify that this is an unique assay
                if reaction not in reactionList:
                    allProducts = modelFile.getReaction(reaction[0]).getProductIds()
                    # Select producing reaction
                    if product in allProducts:
                        reactionList.append(reaction[0])
                        # Define precursors and byproducts
                        subs = modelFile.getReaction(reaction[0]).getSubstrateIds()
                        for substrate in subs:
                                newSubstrates.append(substrate)
                        for extra in allProducts:
                            if not extra in oldSubstrates and extra != product:
                                save = True
                                for met in exclusives:
                                    if met in extra:
                                        save = False
                                        break
                                if save is True:
                                    byproducts.append(extra)
                                    print  product + ' with as byproduct\t' + str(modelFile.getSpeciesIds(extra))
    
    newSubstrates = list(set(newSubstrates)-set(oldSubstrates))
    oldSubstrates = oldSubstrates + newSubstrates

    return reactionList, byproducts, newSubstrates, oldSubstrates


def createSink(modelFile, to_sink, fileName, originalName, objMinFactor):
    
    modelFileCommon = modelFile
    sinkR = modelFile.getReactionIds('R_EX')
    byproductsFinal = []
    sinkM = {}
    sinkIDs = []
    
    # Define metabolites in cytoplasm with exchange reaction
    for exR in sinkR:
        reaction = modelFile.getReaction(exR)
        sinkM[reaction.getSpeciesIds()[0][:-2]+'_c'] = exR
    
    # Verify if exchange reaction excist for products of interest
    # Otherwise create sink reaction
    for byproduct in to_sink:
        if byproduct in sinkM:
            #Verify that the reaction can carry a flux
            testFlux1 = cbm.FluxVariabilityAnalysis(modelFileCommon, selected_reactions = [sinkM[byproduct]], optPercentage = objMinFactor*100)
            if abs(testFlux1[0][0][4]) > 0.0001:
                sinkIDs.append([sinkM[byproduct],originalName])
                byproductsFinal.append(byproduct)

        else:
                modelFileTemp = cbm.CBRead.readSBML3FBC('Models/'+originalName)
                s = modelFileTemp.getSpecies(byproduct)
                speciesname = s.getName()
                sinkid = 'R_EX_' + byproduct[2:]
                sinkname = speciesname + '_sink'
                
                # Make model for single sink
                singleFileName = fileName + '_' + byproduct[2:] +'.xml'
                modelFileTemp.createReaction(sinkid, name = sinkname, reversible=True, create_default_bounds=True)
                modelFileTemp.createReactionReagent(sinkid, byproduct, -1)
                modelFileTemp.setReactionBounds(sinkid,0,999999)
                modelFileTemp.getReaction(sinkid).setAnnotation('SUBSYSTEM','Exchange reactions')
                modelFileTemp.getReaction(sinkid).setAnnotation('NAME',sinkname)

                cbm.CBSolver.analyzeModel(modelFileTemp, oldlpgen=False)
                
                # Verify that the reaction can carry a flux
                testFlux2 = cbm.FluxVariabilityAnalysis(modelFileTemp, selected_reactions = [sinkid], optPercentage = objMinFactor*100)

                if abs(testFlux2[0][0][4]) > 0.0001:
                    cbm.CBWrite.writeSBML3FBC(modelFileTemp, singleFileName, directory='Models/SynthSinks')
                    sinkIDs.append([sinkid,singleFileName])
                    byproductsFinal.append(byproduct)

    return sinkIDs, byproductsFinal

def findProduct(modelFile, biomass_reaction, exclusives, objMinFactor):
    
    reactionList = []
    oldSubstrates = []
    
    originalName = modelFile
    fileName = modelFile[:-4]
    # Write orignal model file in common directory for later use
    modelFile = cbm.CBRead.readSBML3FBC('Models/'+modelFile)
    cbm.CBWrite.writeSBML3FBC(modelFile, originalName, directory='Models/SynthSinks')
    
    # Define co-factors that should not be followed
    bioSubstrates = modelFile.getReaction(biomass_reaction).getSubstrateIds()
    for metabolite in modelFile.getSpeciesIds():
        if not 'C' in modelFile.getSpecies(metabolite).getChemFormula() and metabolite not in bioSubstrates:
            exclusives.append(metabolite[:-1])

    # Take biomass compounds as first products to be analysed
    reactionList, byproducts, newSubstrates, oldSubstrates = getByProd(modelFile, bioSubstrates, oldSubstrates, reactionList, exclusives)

    # Repeat procedure for precursor metabolites found
    while len(newSubstrates) > 0:
        reactionList, byproducts, newSubstrates, oldSubstrates = getByProd(modelFile, newSubstrates, oldSubstrates, reactionList, exclusives)

    # Create a seperate sbml model for each candidate metabolite with its own sink reaction
    sinkIDs, byproductsFinal = createSink(modelFile,list(set(byproducts)), fileName, originalName, objMinFactor)

    return sinkIDs, byproductsFinal



