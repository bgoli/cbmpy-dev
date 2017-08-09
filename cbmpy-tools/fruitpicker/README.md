Fruitpicker
===========

This script was developed by Coco van Boxtel (Vrije Universiteit Amsterdam) in collaboration with Dr. Brett G. Olivier (Vrije Universiteit Amsterdam), 2016

In order to be able to run the script, the following packages need to be installed and fully working:

 - Python 2	<www.python.org>
 - PySCeS CBMPy <www.cbmpy.sourceforge.net> 
 - IBM ILOG CPLEX Optimisation Studio <www.ibm.com>

Procedure:

1. Get your genome scale model and save it with its uptake reactions constrained to the growth condition of interest
2. Copy the model into the folder `Models`
3. run `Run.py`
4. You will be prompted to give info on filename (including extension `.xml`), the biomass reaction, growth rate maintenance (fraction of maximum biomass production rate that should be maintained (0 to 1)) and optionally the export reaction of a product of interest (classic Optknock implementation). 
5. Get the results from the folder `Results`, once the run is finished. 
6. Get the extended models (with synthetic export reactions) from the folder `Models`.

Note: 

Whenever the production rate of your compound of interest is 0.0, this means that the script is not able to interpret the gene annotations of your model correctly. The script is however still able to run for reaction knockouts. 

(c) Coco and Brett, Amsterdam 2017
