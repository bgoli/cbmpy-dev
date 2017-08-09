"""
    
This script verifies the optimisation results returned by the Optknock implementation, in terms 
of flux variance analysis on the metabolite of interest with the suggested gene/reaction knockout.

    
"""


import pyscescbm as cbm
import sys

def testKnockouts(modelFile,fuel,biomass,USE_GENE, delGen,delReact):

    model = cbm.CBRead.readSBML3FBC('Models/'+modelFile)
    model.splitEqualityFluxBounds()
    model.createGeneAssociationsFromAnnotations()


    # Apply Knockouts
    if USE_GENE == True:
        for g_ in delGen:
            model.setGeneInactive(g_, update_reactions=True)
    else:
        for r_ in delReact:
            model.setReactionBounds(r_, 0.0, 0.0)


    cbm.CBSolver.analyzeModel(model)

    ObjValue = model.getReaction(biomass).getValue()

    model.setReactionLowerBound(biomass, ObjValue)
    model.setReactionUpperBound(biomass, ObjValue)

    cbm.CBSolver.analyzeModel(model, oldlpgen=False)

    test = cbm.FluxVariabilityAnalysis(model, selected_reactions = [fuel])

    return test[0][0][2],test[0][0][3],ObjValue