#########################################################################
#This script can be used to model CAM across a wide range of VcVo ratios#
#and compare the metabolic productivity at different NGAM values	#
#									#
#########################################################################

from studyFunctions import *

################################ MAIN ######################################

#import libraries
from libsbml import readSBML
from cobra import io,flux_analysis
import re

#import model from SBML file
cobra_model = io.sbml.create_cobra_model_from_sbml_file("diel_biomass.xml")

#update metabolite charges
fin = open("FractionalCharges.csv","r")

ChargeDict=dict()
for line in fin:
  met=line.replace("\n","").split("\t")[0]
  met = met.replace("-","_")
  charge = line.replace("\n","").split("\t")[1]
  ChargeDict[met]=charge

fin.close()

for met in cobra_model.metabolites:
  tempMet=met.id
  if(met.id[len(met.id)-1]=="2" or met.id[len(met.id)-1]=="1"):
    tempMet = met.id[0:len(met.id)-1]
  if(ChargeDict.keys().__contains__(tempMet)):
    met.charge = ChargeDict.get(tempMet)
  if met.charge is None:
    met.charge=0


#add reactions to represent day-night metabolite accumulation
from cobra.core import Metabolite, Reaction

#Adding transfer reactions
tmfile = open("MetabolitesToTransfer.txt","r")
tmset=set()
for line in tmfile:
  tmset.add(line.replace("\n",""))

for met in tmset:
  tempRxn = Reaction(met+"_dielTransfer")
  tempRxn.add_metabolites({cobra_model.metabolites.get_by_id(met+"1"):-1,cobra_model.metabolites.get_by_id(met+"2"):1})
  tempRxn.lower_bound=-1000
  if not ((met == "STARCH_p") or (met == "SUCROSE_v") or (met == "MAL_v") or (met == "aMAL_v") or (met == "NITRATE_v") or (met == "CIT_v") or (met == "aCIT_v") or (met == "PROTON_v")):
    tempRxn.lower_bound=0
  tempRxn.upper_bound=1000
  cobra_model.add_reaction(tempRxn)


fractionMets=dict()
for rxn in cobra_model.reactions:
  for met in rxn.metabolites.keys():
    a=re.search("^a{1,3}",met.id)
    anion=""
    if a:
      anion=a.group(0)
    b=re.search("^b{1,3}",met.id)
    basic=""
    if b:
      basic=b.group(0)
    prefix = anion
    if prefix == "":
      prefix = basic
    if (abs(rxn.metabolites.get(met)) % 1 > 0 and (not prefix == "") and met.compartment == "v1"):
      fractionMets[met]=prefix


temp=cobra_model.copy()
for met in fractionMets.keys():
  for rxn in met.reactions:
    if rxn.id.__contains__("_dielTransfer"):
      continue
    else:
      mainMet = met.id[len(fractionMets[met]):]
      coeff1 = temp.reactions.get_by_id(rxn.id).metabolites.get(temp.metabolites.get_by_id(mainMet))
      coeff2 = temp.reactions.get_by_id(rxn.id).metabolites.get(temp.metabolites.get_by_id(met.id))
      total = coeff1 + coeff2
      coeff1 = float(coeff1)/total
      coeff2 = float(coeff2)/total
      if cobra_model.reactions.has_id(met.id[0:len(met.id)-1]+"_dielTransfer"):
	temp.reactions.get_by_id(met.id[0:len(met.id)-1]+"_dielTransfer").remove_from_model()
	temp.reactions.get_by_id(mainMet[0:len(mainMet)-1]+"_dielTransfer").remove_from_model()
	Reac = Reaction(mainMet[0:len(mainMet)-1]+"_dielTransfer",name=mainMet+"_dielTransfer")
	Reac.add_metabolites({temp.metabolites.get_by_id(met.id[0:len(met.id)-1]+"1"):-coeff2,temp.metabolites.get_by_id(met.id[0:len(met.id)-1]+"2"):coeff2,temp.metabolites.get_by_id(mainMet[0:len(mainMet)-1]+"1"):-coeff1,temp.metabolites.get_by_id(mainMet[0:len(mainMet)-1]+"2"):coeff1})
	Reac.lower_bound=-1000
	Reac.upper_bound=1000
	temp.add_reaction(Reac)
	print Reac.reaction
      break
   

cobra_model = temp.copy()



########################### SETUP CAM LEAF #####################################################
CAM_model=cobra_model.copy()

#leaf light constraints
CAM_model.reactions.get_by_id("diel_biomass").objective_coefficient=1
CAM_model.reactions.get_by_id("Sucrose_tx1").lower_bound=0
CAM_model.reactions.get_by_id("Sucrose_tx1").upper_bound=0
CAM_model.reactions.get_by_id("GLC_tx1").lower_bound=0
CAM_model.reactions.get_by_id("GLC_tx1").upper_bound=0
CAM_model.reactions.get_by_id("CO2_tx1").lower_bound=0
CAM_model.reactions.get_by_id("CO2_tx1").upper_bound=0
CAM_model.reactions.get_by_id("NH4_tx1").lower_bound=0
CAM_model.reactions.get_by_id("NH4_tx1").upper_bound=0
#leaf dark constraints
CAM_model.reactions.get_by_id("Sucrose_tx2").lower_bound=0
CAM_model.reactions.get_by_id("Sucrose_tx2").upper_bound=0
CAM_model.reactions.get_by_id("GLC_tx2").lower_bound=0
CAM_model.reactions.get_by_id("GLC_tx2").upper_bound=0
CAM_model.reactions.get_by_id("Photon_tx2").lower_bound=0
CAM_model.reactions.get_by_id("Photon_tx2").upper_bound=0
CAM_model.reactions.get_by_id("NH4_tx2").lower_bound=0
CAM_model.reactions.get_by_id("NH4_tx2").upper_bound=0



#nitrate uptake constrain
CAM_model.reactions.get_by_id("Nitrate_ec1").lower_bound=0
CAM_model.reactions.get_by_id("Nitrate_ec1").upper_bound=0

#Vc/Vo not set here
#Rubisco_balance = Metabolite("rubisco_bal_p1", name = "Weights to balance RuBP carboxygenase oxygenase balance", compartment = "p1")
#CAM_model.reactions.get_by_id("RXN_961_p1").add_metabolites({Rubisco_balance:5.15})
#CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1").add_metabolites({Rubisco_balance:-1})
#CAM_model.reactions.get_by_id("RXN_961_p2").add_metabolites({Rubisco_balance:5.15})
#CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2").add_metabolites({Rubisco_balance:-1})


#Remove PTOX activity
CAM_model.reactions.get_by_id("Plastoquinol_Oxidase_p1").lower_bound=0
CAM_model.reactions.get_by_id("Plastoquinol_Oxidase_p1").upper_bound=0



#constraint ATPase:NADPHoxidase = 3:1 and also ATPase (day)=ATPase(night)
Maintenance_constraint = Metabolite("ATPase_NADPHoxidase_constraint_c1",name =  "ATPase_NADPHoxidase_constraint_c1", compartment = "c1")
Maintenance_constraint2 = Metabolite("ATPase_NADPHoxidase_constraint_c2",name =  "ATPase_NADPHoxidase_constraint_c2", compartment = "c2")
Maintenance_constraint3 = Metabolite("Light_dark_maintainence_constraint",name =  "Light_dark_maintainence_constraint", compartment = "c1")
CAM_model.reactions.get_by_id("ATPase_tx1").add_metabolites({Maintenance_constraint:1,Maintenance_constraint3:1})
CAM_model.reactions.get_by_id("ATPase_tx2").add_metabolites({Maintenance_constraint2:1,Maintenance_constraint3:-1})
#CAM_model.reactions.get_by_id("ATPase_tx1").add_metabolites({Maintenance_constraint:1})
#CAM_model.reactions.get_by_id("ATPase_tx2").add_metabolites({Maintenance_constraint2:1})
CAM_model.reactions.get_by_id("NADPHoxc_tx1").add_metabolites({Maintenance_constraint:-3})
CAM_model.reactions.get_by_id("NADPHoxc_tx2").add_metabolites({Maintenance_constraint2:-3})
CAM_model.reactions.get_by_id("NADPHoxm_tx1").add_metabolites({Maintenance_constraint:-3})
CAM_model.reactions.get_by_id("NADPHoxm_tx2").add_metabolites({Maintenance_constraint2:-3})
CAM_model.reactions.get_by_id("NADPHoxp_tx1").add_metabolites({Maintenance_constraint:-3})
CAM_model.reactions.get_by_id("NADPHoxp_tx2").add_metabolites({Maintenance_constraint2:-3})

#constrain sucrose, fructan and starch storage
CAM_model.reactions.get_by_id("SUCROSE_v_dielTransfer").lower_bound=0
CAM_model.reactions.get_by_id("SUCROSE_v_dielTransfer").upper_bound=0
CAM_model.reactions.get_by_id("FRUCTAN_v_dielTransfer").lower_bound=0
CAM_model.reactions.get_by_id("FRUCTAN_v_dielTransfer").upper_bound=0


#Setting chloroplastic NADPH dehydrogenase to 0  ((Yamamoto et al., 2011)
CAM_model.reactions.get_by_id("NADPH_Dehydrogenase_p1").lower_bound=0
CAM_model.reactions.get_by_id("NADPH_Dehydrogenase_p1").upper_bound=0
CAM_model.reactions.get_by_id("NADPH_Dehydrogenase_p2").lower_bound=0
CAM_model.reactions.get_by_id("NADPH_Dehydrogenase_p2").upper_bound=0


#Set biomass to zero
CAM_model.reactions.get_by_id("Biomass_tx1").lower_bound=0
CAM_model.reactions.get_by_id("Biomass_tx1").upper_bound=0
CAM_model.reactions.get_by_id("Biomass_tx2").lower_bound=0
CAM_model.reactions.get_by_id("Biomass_tx2").upper_bound=0

#constraining light to 200
CAM_model.reactions.get_by_id("Photon_tx1").lower_bound=200
CAM_model.reactions.get_by_id("Photon_tx1").upper_bound=200

#get solution to check if model works
#solution=flux_analysis.parsimonious.optimize_minimal_flux(CAM_model)

#set pH in night time vacuole to 3.3
temp = convertToClassicalModel(CAM_model,comp="v2",updateCharges = "FractionalCharges.csv")
CAM_model_pH = convertToFractionalCharges(temp,"VacuolarMetabolitesAtpH3DOT3.csv")

#perform parsimonious optimization
solution=flux_analysis.parsimonious.optimize_minimal_flux(CAM_model_pH)

#constrain for starch phosphorylating CAM
CAM_model_pH.reactions.get_by_id("MALTODEXGLUCOSID_RXN_p2").lower_bound=0
CAM_model_pH.reactions.get_by_id("MALTODEXGLUCOSID_RXN_p2").upper_bound=0
CAM_model_pH.reactions.get_by_id("RXN_1827_p2").lower_bound=0
CAM_model_pH.reactions.get_by_id("RXN_1827_p2").upper_bound=0

#perform parsimonious optimization
solution=flux_analysis.parsimonious.optimize_minimal_flux(CAM_model_pH)


########################### SETUP C3 LEAF #####################################################

C3_model=cobra_model.copy()
#leaf light constraints
C3_model.reactions.get_by_id("diel_biomass").objective_coefficient=1
C3_model.reactions.get_by_id("Sucrose_tx1").lower_bound=0
C3_model.reactions.get_by_id("Sucrose_tx1").upper_bound=0
C3_model.reactions.get_by_id("GLC_tx1").lower_bound=0
C3_model.reactions.get_by_id("GLC_tx1").upper_bound=0
C3_model.reactions.get_by_id("CO2_tx1").lower_bound=0
C3_model.reactions.get_by_id("NH4_tx1").lower_bound=0
C3_model.reactions.get_by_id("NH4_tx1").upper_bound=0
#leaf dark constraints
C3_model.reactions.get_by_id("Sucrose_tx2").lower_bound=0
C3_model.reactions.get_by_id("Sucrose_tx2").upper_bound=0
C3_model.reactions.get_by_id("GLC_tx2").lower_bound=0
C3_model.reactions.get_by_id("GLC_tx2").upper_bound=0
C3_model.reactions.get_by_id("Photon_tx2").lower_bound=0
C3_model.reactions.get_by_id("Photon_tx2").upper_bound=0
C3_model.reactions.get_by_id("NH4_tx2").lower_bound=0
C3_model.reactions.get_by_id("NH4_tx2").upper_bound=0
C3_model.reactions.get_by_id("CO2_tx2").upper_bound=0



from cobra.core import Metabolite, Reaction


#nitrate uptake constraint
Nitrate_balance = Metabolite("Nitrate_bal_c", name = "Weights to balance nitrate uptake", compartment = "c1")
C3_model.reactions.get_by_id("Nitrate_ec1").add_metabolites({Nitrate_balance:-2})
C3_model.reactions.get_by_id("Nitrate_ec2").add_metabolites({Nitrate_balance:3})


#constraint VcVo
Rubisco_balance = Metabolite("rubisco_bal_p1", name = "Weights to balance RuBP carboxygenase oxygenase balance", compartment = "p1")
C3_model.reactions.get_by_id("RXN_961_p1").add_metabolites({Rubisco_balance:3})#changed from Alison Smith's 3 to Leegood's 2.5 to back to back to 3 based on Ma et al 2014 ()
C3_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1").add_metabolites({Rubisco_balance:-1})
#C3_model.reactions.get_by_id("RXN_961_p2").add_metabolites({Rubisco_balance:3})
#C3_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2").add_metabolites({Rubisco_balance:-1})

#constraint ATPase:NADPHoxidase = 3:1 and also ATPase (day)=ATPase(night)
Maintenance_constraint = Metabolite("ATPase_NADPHoxidase_constraint_c1",name =  "ATPase_NADPHoxidase_constraint_c1", compartment = "c1")
Maintenance_constraint2 = Metabolite("ATPase_NADPHoxidase_constraint_c2",name =  "ATPase_NADPHoxidase_constraint_c2", compartment = "c2")
Maintenance_constraint3 = Metabolite("Light_dark_maintainence_constraint",name =  "Light_dark_maintainence_constraint", compartment = "c1")
C3_model.reactions.get_by_id("ATPase_tx1").add_metabolites({Maintenance_constraint:1,Maintenance_constraint3:1})
C3_model.reactions.get_by_id("ATPase_tx2").add_metabolites({Maintenance_constraint2:1,Maintenance_constraint3:-1})
#C3_model.reactions.get_by_id("ATPase_tx1").add_metabolites({Maintenance_constraint:1})
#C3_model.reactions.get_by_id("ATPase_tx2").add_metabolites({Maintenance_constraint2:1})
C3_model.reactions.get_by_id("NADPHoxc_tx1").add_metabolites({Maintenance_constraint:-3})
C3_model.reactions.get_by_id("NADPHoxc_tx2").add_metabolites({Maintenance_constraint2:-3})
C3_model.reactions.get_by_id("NADPHoxm_tx1").add_metabolites({Maintenance_constraint:-3})
C3_model.reactions.get_by_id("NADPHoxm_tx2").add_metabolites({Maintenance_constraint2:-3})
C3_model.reactions.get_by_id("NADPHoxp_tx1").add_metabolites({Maintenance_constraint:-3})
C3_model.reactions.get_by_id("NADPHoxp_tx2").add_metabolites({Maintenance_constraint2:-3})

#Remove PTOX activity
C3_model.reactions.get_by_id("Plastoquinol_Oxidase_p1").lower_bound=0
C3_model.reactions.get_by_id("Plastoquinol_Oxidase_p1").upper_bound=0

#Plastid enolase was not detected in Arabidopsis mesophyll tissue
C3_model.reactions.get_by_id("2PGADEHYDRAT_RXN_p1").lower_bound=0
C3_model.reactions.get_by_id("2PGADEHYDRAT_RXN_p1").upper_bound=0
C3_model.reactions.get_by_id("2PGADEHYDRAT_RXN_p2").lower_bound=0
C3_model.reactions.get_by_id("2PGADEHYDRAT_RXN_p2").upper_bound=0

#constrain sucrose,fructan and starch storage
C3_model.reactions.get_by_id("SUCROSE_v_dielTransfer").lower_bound=0
C3_model.reactions.get_by_id("SUCROSE_v_dielTransfer").upper_bound=0
C3_model.reactions.get_by_id("FRUCTAN_v_dielTransfer").lower_bound=0
C3_model.reactions.get_by_id("FRUCTAN_v_dielTransfer").upper_bound=0

#Setting chloroplastic NADPH dehydrogenase to 0  ((Yamamoto et al., 2011)
C3_model.reactions.get_by_id("NADPH_Dehydrogenase_p1").lower_bound=0
C3_model.reactions.get_by_id("NADPH_Dehydrogenase_p1").upper_bound=0
C3_model.reactions.get_by_id("NADPH_Dehydrogenase_p2").lower_bound=0
C3_model.reactions.get_by_id("NADPH_Dehydrogenase_p2").upper_bound=0

#Set biomass to zero
C3_model.reactions.get_by_id("Biomass_tx1").lower_bound=0
C3_model.reactions.get_by_id("Biomass_tx1").upper_bound=0
C3_model.reactions.get_by_id("Biomass_tx2").lower_bound=0
C3_model.reactions.get_by_id("Biomass_tx2").upper_bound=0

#constraining light to 200
C3_model.reactions.get_by_id("Photon_tx1").lower_bound=200
C3_model.reactions.get_by_id("Photon_tx1").upper_bound=200

#constraint C3 as starch hydrolyzing
C3_model.reactions.get_by_id("RXN_1826_p2").lower_bound=0
C3_model.reactions.get_by_id("RXN_1826_p2").upper_bound=0
solution=flux_analysis.parsimonious.optimize_minimal_flux(C3_model)

#perform parsimonious optimization for maintenance ATPase = 0, 4, 8, 12, 16
C3_sols = dict()
for i in range(0,17,4):
	C3_model.reactions.get_by_id("ATPase_tx2").lower_bound = i
	C3_model.reactions.get_by_id("ATPase_tx2").upper_bound = i
	solution=flux_analysis.parsimonious.optimize_minimal_flux(C3_model)
	C3_sols[i] = solution.f*4*12*60*60




################################### Vc/Vo CURVE FOR CAM WITH DIFFERENT NGAM VALUES #############


test_model=CAM_model_pH.copy()
test_model.reactions.get_by_id("MALIC_NADP_RXN_c1").lower_bound=0
test_model.reactions.get_by_id("MALIC_NADP_RXN_c1").upper_bound=0
test_model.reactions.get_by_id("MALIC_NADP_RXN_p1").lower_bound=0
test_model.reactions.get_by_id("MALIC_NADP_RXN_p1").upper_bound=0
test_model.reactions.get_by_id("MALIC_NAD_RXN_m1").lower_bound=0
test_model.reactions.get_by_id("MALIC_NAD_RXN_m1").upper_bound=0
test_model.reactions.get_by_id("STARCH_p_dielTransfer").lower_bound = -1000
test_model.reactions.get_by_id("STARCH_p_dielTransfer").upper_bound = 1000
test_model.reactions.get_by_id("SUCROSE_v_dielTransfer").lower_bound = 0
test_model.reactions.get_by_id("SUCROSE_v_dielTransfer").upper_bound = 0
test_model.reactions.get_by_id("FRUCTAN_v_dielTransfer").lower_bound=0
test_model.reactions.get_by_id("FRUCTAN_v_dielTransfer").upper_bound=0
solution=flux_analysis.parsimonious.optimize_minimal_flux(test_model)
backup1 = test_model.copy()


ansA = dict()

for j in range(0,17,4):
  test_model = backup1.copy()
  test_model.reactions.get_by_id("ATPase_tx2").lower_bound = j
  test_model.reactions.get_by_id("ATPase_tx2").upper_bound = j
  solution=flux_analysis.parsimonious.optimize_minimal_flux(test_model)
  backup2 = test_model.copy()
  ansB = dict()
  solutiondict=dict()
  ansB[10000]=test_model.solution.f*4*12*60*60
  solutiondict[10000] = test_model.solution.x_dict
  for i in [1000,160,80,40,20,10,5,3,2]:
    b = float(i)
    print b
    Rubisco_balance = Metabolite("rubisco_bal_p1", name = "Weights to balance RuBP carboxygenase oxygenase balance", compartment = "p1")
    
    temp = backup2.copy()
    temp.reactions.get_by_id("STARCH_p_dielTransfer").lower_bound = -1000
    temp.reactions.get_by_id("STARCH_p_dielTransfer").upper_bound = 1000
    temp.reactions.get_by_id("SUCROSE_v_dielTransfer").lower_bound = 0
    temp.reactions.get_by_id("SUCROSE_v_dielTransfer").upper_bound = 0
    temp.reactions.get_by_id("FRUCTAN_v_dielTransfer").lower_bound=0
    temp.reactions.get_by_id("FRUCTAN_v_dielTransfer").upper_bound=0
    temp.reactions.get_by_id("RXN_961_p1").add_metabolites({Rubisco_balance:-b})
    temp.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1").add_metabolites({Rubisco_balance:1})
    solution=flux_analysis.parsimonious.optimize_minimal_flux(temp)
    
    test_model=temp.copy()
    #temp_model.reactions.get_by_id("diel_biomass").lower_bound=0.64
    #temp_model.reactions.get_by_id("diel_biomass").upper_bound=0.64
    #solution=flux_analysis.parsimonious.optimize_minimal_flux(temp_model)
    print str(temp.solution.x_dict.get("RXN_961_p1"))+"\t"+str(temp.solution.x_dict.get("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1"))
    print "Storing results from "+str(b)
    ansB[b]=solution.f*4*12*60*60
    solutiondict[b]=test_model.solution.x_dict
    print ansB[b]
    #break
  writeAllFluxes(CAM_model_pH,solutiondict,"CAM_PEPCK_starch_maintenance"+str(j)+".csv")
  ansA[j] = ansB

######################################Plot metabolic productivities#####################

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 18}) #sets a global fontsize
plt.rcParams['xtick.major.size'] = 5 # adjusts tick line length and width
plt.rcParams['xtick.major.width'] = 1
plt.rcParams['ytick.major.size'] = 5
plt.rcParams['ytick.major.width'] = 1
plt.rcParams['axes.linewidth']=3 # makes axes line thicker
#plt.figure(figsize=(5,15))

xvalues=list()
xvalues2=list()
yvalues=list()
yvalues2=list()
yvalues3=list()
yvalues4=list()
yvalues5=list()
yvalues6=list()
yvalues7=list()
yvalues8=list()
for k in sorted(ansA[0].keys()):
  if k==10000:
    xvalues.append(r'$\infty$')
  xvalues.append(int(k))
  xvalues2.append(np.log(k))
  yvalues.append(ansA[0].get(k)/1000)
  yvalues2.append(ansA[4].get(k)/1000)
  yvalues3.append(ansA[8].get(k)/1000)
  yvalues4.append(ansA[12].get(k)/1000)
  yvalues5.append(ansA[16].get(k)/1000)

curve = plt.subplot(121)
curve.plot(xvalues2,yvalues,marker="o",label="ATPase = 0",color="green")
curve.plot(xvalues2,yvalues2,marker="o",label="ATPase = 4",color="red")
curve.plot(xvalues2,yvalues3,marker="o",label="ATPase = 8",color="blue")
curve.plot(xvalues2,yvalues4,marker="o",label="ATPase = 12",color="purple")
curve.plot(xvalues2,yvalues5,marker="o",label="ATPase = 16",color="orange")
curve.axvspan(np.log(float(1)/0.33),np.log(float(1)/0.14),alpha=0.3,color="green")
curve.plot(xvalues2[1],C3_sols[0]/1000,marker="s",markersize=10,color="green")
curve.axhline(y=C3_sols[0]/1000,ls="--",alpha=0.3,color="green")
curve.plot(xvalues2[1],C3_sols[4]/1000,marker="s",markersize=10,color="red")
curve.axhline(y=C3_sols[4]/1000,ls="--",alpha=0.3,color="red")
curve.plot(xvalues2[1],C3_sols[8]/1000,marker="s",markersize=10,color="blue")
curve.axhline(y=C3_sols[8]/1000,ls="--",alpha=0.3,color="blue")
curve.plot(xvalues2[1],C3_sols[12]/1000,marker="s",markersize=10,color="purple")
curve.axhline(y=C3_sols[12]/1000,ls="--",alpha=0.3,color="purple")
curve.plot(xvalues2[1],C3_sols[16]/1000,marker="s",markersize=10,color="orange")
curve.axhline(y=C3_sols[16]/1000,ls="--",alpha=0.3,color="orange")
curve.set_xlabel(r'$Vc/Vo$',fontsize=20)
curve.set_xticks(xvalues2)
curve.set_xticklabels(xvalues)
curve.set_ylabel("Phloem output ("+r'$mmol.m^{-2}.d^{-1}$'+")",fontsize=20)


#plt.legend(loc="best")
plt.legend(bbox_to_anchor = (1.65,1.03))
plt.show()

