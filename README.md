# CAM-study-scripts
This repository contains scripts used to compare the energetics of CAM and C3 systems in fully developed leaves

## Dependencies:
libsbml version 5.12.1  
[cobrapy version 0.5.4](https://github.com/opencobra/cobrapy/tree/73ef5623ad)  

## Instructions:
1) Install python version 2, libsbml and cobrapy
2) Go to folder containing scripts
3) Open a console(mac), terminal(linux) or command prompt(windows)  
4) Run scripts (two options available)  
a) using the command `python scriptname.py` where scriptname.py is the filename of the script  
<p align='center'>
  OR
</p>
<p>
  &nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp b) Enter python terminal using the command "python" and copy paste script contents to the python terminal
</p>

  
The latter method is recommended as it maintains the data in the python environment allowing one to analyse metabolic fluxes directly and make your own figures.

## Index
* Script1.py - estimating non-growth associated maintenance (NGAM) in a fully developed C3 leaf  
* Script2.py - modelling C3 and CAM systems: PPFD = 200 umol/m2/s  
* Script3.py - modelling CAM in Opuntia ficus indica  
* Script4.py - modelling C3 and CAM systems with fixed metabolic output  
* Script5.py - comparing productivity of CAM and C3 systems  
* Script6.py - comparing productivity of CAM and C3 systems with varying NGAM  
* Script7.py - estimating NGAM in a fully developed C3 leaf in 18h6h day  
* Script8.py - estimating NGAM in a fully developed C3 leaf in 6h18h day  
* Script9.py - comparing productivity of CAM and C3 systems under different photoperiod conditions  
* Script10.py - estimating NGAM in a fully developed C3 leaf when sugar:amino acids in the phloem = 0.5  
* Script11.py - estimating NGAM in a fully developed C3 leaf when sugar:amino acids in the phloem = 0.7  
* Script12.py - estimating NGAM in a fully developed C3 leaf when sugar:amino acids in the phloem = 0.99  
* Script13.py - comparing productivity of CAM and C3 systems with varying sugar:amino acid ratio in the phloem  
* Script14.py - estimating NGAM in a fully developed C3 leaf when day-night phloem export rate ratio is 3:2  
* Script15.py - estimating NGAM in a fully developed C3 leaf when day-night phloem export rate ratio is 1:1  
* Script16.py - comparing productivity of CAM and C3 systems with varying day-night phloem export rate ratio 
* Script17.py - comparing productivity of fructan-storing CAM with symplastic and apoplastic loading
* diel_core_model - a model representing day and night time metabolism, that have yet to be configured for CAM or C3 metabolism  
* C3_model - a model representing primary metabolism in C3 leaf with constraints matching that described in the paper  
* CAM_model - a model representing primary metabolism in CAM leaf with constraints matching that described in the paper  
* rxnpathwaydict - a file representing reactions and associated pathways which is required to print our pathway information in outputs  
* pathwayprioritizer - a file listing pathways and a value for each pathway 1,2,3 ... which can be used to prioritize reactions in the output file for easy interpretation  
* FractionalCharges - a file containing metabolite charges (to complement charge attribute in the model)  
* MetabolitesToTransfer - a file listing metabolites that will be allowed to accumulate over the diel cycle  
* VacuolarMetabolitesAtpH3DOT3 - a file listing metabolites that exists as fractional charges at pH=3.3
* OpuntiaBiomass - phloem composition reported in Opuntia ficus indica
