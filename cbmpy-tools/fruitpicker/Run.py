"""
    MetaboliteFinder & pyOpt 1.0
    ==============
    
    A Python implementation of the OptKnock algorithm, extended with a metabolite candidate finder.
    Copyright (C) 2016 Coco van Boxtel and Brett G. Olivier, VU University Amsterdam, Amsterdam, The Netherlands
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>
    
    Authors: Coco van Boxtel and Brett G. Olivier
    Contact email: c.van.boxtel@vu.com
    Contact email: bgoli@users.sourceforge.net
    
"""

import os, sys
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
import pyOpt, productSelection, testIndependently
import pyscescbm as cbm

# ------
# Extra User parameters

# Define value that indicates 'unbounded' in sbml file, to be translated to cplex.inf
# If boundaries are already defined infinite, make sure this value is high
infinityValue = 99999
    
# Reduce number of knockouts
# Keep alpha value small enough to avoid potential dilution of the bilevel objective
# Can only be used if bilevelObjective is maximization
USE_KNOCKOUT_WEIGHTING = False
KNOCKOUT_WEIGHTING_ALPHA = 0.0002
    
# Allow a near optimal solution (e.g. 5% less than maximum)
# This will increase the solution speed.
SOLUTION_FROM_OPTIMUM = 0.0

# Define co-factors
exclusives = ['M_atp_','M_adp_','M_nadph_','M_nadh_','M_nad_', 'M_nadp_','M_ACP_','M_amp_','M_co2_','M_coa_','M_accoa_']

# ------


# Prompt
print 'Give name of file, located in the folder "Models"'
while True:
    modelFile = raw_input('> ')
    try:
        test = cbm.CBRead.readSBML3FBC('Models/'+modelFile)
        break
    except:
        print '\nThis file does not excist in the "Model" directory. Please enter the name of a SBML3 file including extension.'

print 'What is the reaction id of the biomass forming reaction?'
while True:
    bio_reaction = raw_input('> ')
    try:
        test.getReaction(bio_reaction).getName()
        break
    except:
        print '\nThis reaction ID does not excist in the given model.'


print '\nWhat fraction of biomass production do you want to maintain?'
while True:
    try:
        objMinFactor = float(raw_input('> '))
        if 0.0 <= objMinFactor <= 1.0:
            break
        else:
            print 'Fractions are between 0.0 and 1.0.'
    except:
        print 'Please give a number'
        continue

print '\nWhat is the maximum amount of knockouts you would like to have?'
while True:
    try:
        maxDelete = float(raw_input('> '))
        if 0.0 <= maxDelete:
            break
        else:
            print 'Fractions are between 0.0 and 1.0.'
    except:
        print 'Please give a number'
        continue

print '\nDo you want to know gene knockouts?'
doGenes = raw_input('> ')
if doGenes == 'Yes' or doGenes =='yes':
    USE_GENE = True
elif doGenes == 'No' or doGenes =='no':
    USE_GENE = False
else:
    USE_GENE = False
    print '\n!!!Reaction knockouts are assumed!!!'

print '\nDo you automatically want to search for product candidates?'
SingleRun = False
autoProd = raw_input('> ')
if autoProd == 'Yes' or autoProd =='yes':
    bilevelObjectives, byproducts = productSelection.findProduct(modelFile, bio_reaction, exclusives, objMinFactor)
else:
    SingleRun = True
    print '\nWhich reaction do you want to optimize for then?'
    while True:
        bilevelObjectives = raw_input('> ')
        try:
            test.getReaction(bilevelObjectives).getName()
            bilevelObjectives = [bilevelObjectives]
            break
        except:
            print '\nThis reaction ID does not excist in the given model.'

# Open new file for results
resultFile = file(os.path.join('Results', modelFile.replace('.xml', '') + '.result.txt'),'w')
resultFile.write('=================\n')
resultFile.write('pyOpt Solutions\n=================\n')
resultFile.write('\nMaximum deletions: {}'.format(maxDelete))
if len(bilevelObjectives) > 1:
    resultFile.write('\nAll metabolites tested: {}'.format(byproducts))
if USE_KNOCKOUT_WEIGHTING:
    resultFile.write('\nKnockout weighting (minimization): ON')
    resultFile.write('\nWARNING! The outer optimum value is slightly lower, due to use of knockout weighting')
resultFile.write('\n\n=================\n\n')

# Run pyOpt for compounds of interest
for reaction in bilevelObjectives:
    if SingleRun == True:
        [bilevelObjective, val1, bio_reaction, val2, delGen,delReact] = pyOpt.runOptKnock(modelFile, (reaction,1), bio_reaction, objMinFactor, maxDelete, infinityValue, USE_GENE, USE_KNOCKOUT_WEIGHTING, KNOCKOUT_WEIGHTING_ALPHA, SOLUTION_FROM_OPTIMUM)
    else:
        [bilevelObjective, val1, bio_reaction, val2, delGen,delReact] = pyOpt.runOptKnock('SynthSinks/'+reaction[1], (reaction[0],1), bio_reaction, objMinFactor, maxDelete, infinityValue, USE_GENE, USE_KNOCKOUT_WEIGHTING, KNOCKOUT_WEIGHTING_ALPHA, SOLUTION_FROM_OPTIMUM)

    if abs(round(val1,5)) != 0.0:
        # test outcome independent from OptKnock implementation
        if SingleRun == True:
            value1,value2,value3 = testIndependently.testKnockouts(modelFile, reaction, bio_reaction, USE_GENE,delGen,delReact)
        else:
            value1,value2,value3 = testIndependently.testKnockouts('SynthSinks/'+reaction[1], reaction[0], bio_reaction, USE_GENE,delGen,delReact)

        # write results in file
        resultFile.write('\nProduct reaction ({}): min {}, max {}'.format(bilevelObjective[0],round(value1,5),round(value2,5)))
        resultFile.write('\nBiomass reaction ({}): {}'.format(bio_reaction,round(value3,5)))
        if USE_GENE:
            resultFile.write('\nDeleted genes ({}): {}'.format(len(delGen), delGen))
        else:
            resultFile.write('\nDeleted reactions ({}): {}'.format(len(delReact), delReact))
        resultFile.write('\n_________________\n')

resultFile.close()







