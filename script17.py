#########################################################################
#This script can be used to model fructan-storing CAM substype	with 	#
#symplastic or apoplastic phloem loading acorss a range of VcVo ratios 	#
#and compare their metabolic productivity	 			#
#									#
#########################################################################

from studyFunctions import *

################################ MAIN ######################################

#import libraries
from libsbml import readSBML
from cobra import io,flux_analysis
import re

#import model from SBML file
cobra_model = io.sbml.create_cobra_model_from_sbml_file("diel_model.xml")


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
CAM_model.reactions.get_by_id("ATPase_tx2").lower_bound=8.5
CAM_model.reactions.get_by_id("ATPase_tx2").upper_bound=8.5

#get solution to check if model works
#solution=flux_analysis.parsimonious.optimize_minimal_flux(CAM_model)

#set pH in night time vacuole to 3.3
temp = convertToClassicalModel(CAM_model,comp="v2",updateCharges = "FractionalCharges.csv")
CAM_model_pH = convertToFractionalCharges(temp,"VacuolarMetabolitesAtpH3DOT3.csv")

#perform parsimonious optimization
solution=flux_analysis.parsimonious.optimize_minimal_flux(CAM_model_pH)


#perform parsimonious optimization of starch hydrolyzing CAM
temp_model=CAM_model_pH.copy()
temp_model.reactions.get_by_id("RXN_1826_p2").lower_bound=0
temp_model.reactions.get_by_id("RXN_1826_p2").upper_bound=0
solution=flux_analysis.parsimonious.optimize_minimal_flux(temp_model)


#perform parsimonious optimization of starch phosphorylating CAM
temp_model2=CAM_model_pH.copy()
temp_model2.reactions.get_by_id("MALTODEXGLUCOSID_RXN_p2").lower_bound=0
temp_model2.reactions.get_by_id("MALTODEXGLUCOSID_RXN_p2").upper_bound=0
temp_model2.reactions.get_by_id("RXN_1827_p2").lower_bound=0
temp_model2.reactions.get_by_id("RXN_1827_p2").upper_bound=0
solution=flux_analysis.parsimonious.optimize_minimal_flux(temp_model2)


########################### SETUP C3 LEAF #####################################################

C3_model=cobra_model.copy()
#Leaves - light
C3_model.reactions.get_by_id("diel_biomass").objective_coefficient=1
C3_model.reactions.get_by_id("Sucrose_tx1").lower_bound=0
C3_model.reactions.get_by_id("Sucrose_tx1").upper_bound=0
C3_model.reactions.get_by_id("GLC_tx1").lower_bound=0
C3_model.reactions.get_by_id("GLC_tx1").upper_bound=0
C3_model.reactions.get_by_id("CO2_tx1").lower_bound=0
C3_model.reactions.get_by_id("NH4_tx1").lower_bound=0
C3_model.reactions.get_by_id("NH4_tx1").upper_bound=0
#Leaves - dark
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

#nitrate uptake constrain
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
C3_model.reactions.get_by_id("ATPase_tx2").lower_bound=8.5
C3_model.reactions.get_by_id("ATPase_tx2").upper_bound=8.5

#get solution to check if model works
solution=flux_analysis.parsimonious.optimize_minimal_flux(C3_model)

#perform parsimonious optimization of starch hydrolyzing C3
temp_model3 = C3_model.copy()
temp_model3.reactions.get_by_id("RXN_1826_p2").lower_bound=0
temp_model3.reactions.get_by_id("RXN_1826_p2").upper_bound=0
solution=flux_analysis.parsimonious.optimize_minimal_flux(temp_model3)


#perform parsimonious optimization of starch phosphorylating C3
temp_model4=C3_model.copy()
temp_model4.reactions.get_by_id("MALTODEXGLUCOSID_RXN_p2").lower_bound=0
temp_model4.reactions.get_by_id("MALTODEXGLUCOSID_RXN_p2").upper_bound=0
temp_model4.reactions.get_by_id("RXN_1827_p2").lower_bound=0
temp_model4.reactions.get_by_id("RXN_1827_p2").upper_bound=0
solution=flux_analysis.parsimonious.optimize_minimal_flux(temp_model4)





############# Vc/Vo CRUVE FOR FRUCTAN-STORING CAM SUBTYPE WITH SYMPLASTIC AND APOPLATIC PHLOEM LOADING #############

###############fructan storing - active Sucrose vacuolar loading - apoplstic Phloem loading
test_model=CAM_model_pH.copy()
test_model.reactions.get_by_id("MALIC_NADP_RXN_c1").lower_bound=0
test_model.reactions.get_by_id("MALIC_NADP_RXN_c1").upper_bound=0
test_model.reactions.get_by_id("MALIC_NADP_RXN_p1").lower_bound=0
test_model.reactions.get_by_id("MALIC_NADP_RXN_p1").upper_bound=0
test_model.reactions.get_by_id("MALIC_NAD_RXN_m1").lower_bound=0
test_model.reactions.get_by_id("MALIC_NAD_RXN_m1").upper_bound=0
test_model.reactions.get_by_id("STARCH_p_dielTransfer").lower_bound = 0
test_model.reactions.get_by_id("STARCH_p_dielTransfer").upper_bound = 0
test_model.reactions.get_by_id("SUCROSE_v_dielTransfer").lower_bound = 0
test_model.reactions.get_by_id("SUCROSE_v_dielTransfer").upper_bound = 0
test_model.reactions.get_by_id("FRUCTAN_v_dielTransfer").lower_bound= -1000
test_model.reactions.get_by_id("FRUCTAN_v_dielTransfer").upper_bound= 1000
solution=flux_analysis.parsimonious.optimize_minimal_flux(test_model)

ansA = dict()
ansA[10000]=test_model.solution.f

solutiondict=dict()
solutiondict[10000]=test_model.solution.x_dict


for i in [1000,160,80,40,20,10,5,3,2]:
  b = float(i)
  print b
  Rubisco_balance = Metabolite("rubisco_bal_p1", name = "Weights to balance RuBP carboxygenase oxygenase balance", compartment = "p1")
  
  temp = CAM_model_pH.copy()
  temp.reactions.get_by_id("STARCH_p_dielTransfer").lower_bound = 0
  temp.reactions.get_by_id("STARCH_p_dielTransfer").upper_bound = 0
  temp.reactions.get_by_id("SUCROSE_v_dielTransfer").lower_bound = 0
  temp.reactions.get_by_id("SUCROSE_v_dielTransfer").upper_bound = 0
  temp.reactions.get_by_id("FRUCTAN_v_dielTransfer").lower_bound= -1000
  temp.reactions.get_by_id("FRUCTAN_v_dielTransfer").upper_bound= 1000
  temp.reactions.get_by_id("RXN_961_p1").add_metabolites({Rubisco_balance:-b})
  temp.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1").add_metabolites({Rubisco_balance:1})
  solution=flux_analysis.parsimonious.optimize_minimal_flux(temp)
  
  test_model=temp.copy()
  #temp_model.reactions.get_by_id("diel_biomass").lower_bound=0.64
  #temp_model.reactions.get_by_id("diel_biomass").upper_bound=0.64
  #solution=flux_analysis.parsimonious.optimize_minimal_flux(temp_model)
  
  print "Storing results from "+str(b)
  ansA[b]=test_model.solution.f
  solutiondict[b]=test_model.solution.x_dict
  #break



writeAllFluxes(CAM_model_pH,solutiondict,"PEPCK_fructan_activeSuc_activePhlo.csv")

###############fructan storing - passive Sucrose vacuolar loading - apoplastic Phloem loading
test_model=CAM_model_pH.copy()
test_model.reactions.get_by_id("MALIC_NADP_RXN_c1").lower_bound=0
test_model.reactions.get_by_id("MALIC_NADP_RXN_c1").upper_bound=0
test_model.reactions.get_by_id("MALIC_NADP_RXN_p1").lower_bound=0
test_model.reactions.get_by_id("MALIC_NADP_RXN_p1").upper_bound=0
test_model.reactions.get_by_id("MALIC_NAD_RXN_m1").lower_bound=0
test_model.reactions.get_by_id("MALIC_NAD_RXN_m1").upper_bound=0
test_model.reactions.get_by_id("STARCH_p_dielTransfer").lower_bound = 0
test_model.reactions.get_by_id("STARCH_p_dielTransfer").upper_bound = 0
test_model.reactions.get_by_id("SUCROSE_v_dielTransfer").lower_bound = 0
test_model.reactions.get_by_id("SUCROSE_v_dielTransfer").upper_bound = 0
test_model.reactions.get_by_id("FRUCTAN_v_dielTransfer").lower_bound= -1000
test_model.reactions.get_by_id("FRUCTAN_v_dielTransfer").upper_bound= 1000

test_model.reactions.get_by_id("SUCROSE_PROTON_vc1").pop(test_model.metabolites.get_by_id("PROTON_v1"))
test_model.reactions.get_by_id("SUCROSE_PROTON_vc2").pop(test_model.metabolites.get_by_id("PROTON_v2"))
test_model.reactions.get_by_id("SUCROSE_PROTON_vc1").pop(test_model.metabolites.get_by_id("PROTON_c1"))
test_model.reactions.get_by_id("SUCROSE_PROTON_vc2").pop(test_model.metabolites.get_by_id("PROTON_c2"))

solution=flux_analysis.parsimonious.optimize_minimal_flux(test_model)

ansB = dict()
ansB[10000]=test_model.solution.f

solutiondict=dict()
solutiondict[10000]=test_model.solution.x_dict


for i in [1000,160,80,40,20,10,5,3,2]:
  b = float(i)
  print b
  Rubisco_balance = Metabolite("rubisco_bal_p1", name = "Weights to balance RuBP carboxygenase oxygenase balance", compartment = "p1")
  
  temp = CAM_model_pH.copy()
  temp.reactions.get_by_id("STARCH_p_dielTransfer").lower_bound = 0
  temp.reactions.get_by_id("STARCH_p_dielTransfer").upper_bound = 0
  temp.reactions.get_by_id("SUCROSE_v_dielTransfer").lower_bound = 0
  temp.reactions.get_by_id("SUCROSE_v_dielTransfer").upper_bound = 0
  temp.reactions.get_by_id("FRUCTAN_v_dielTransfer").lower_bound= -1000
  temp.reactions.get_by_id("FRUCTAN_v_dielTransfer").upper_bound= 1000
  temp.reactions.get_by_id("RXN_961_p1").add_metabolites({Rubisco_balance:-b})
  temp.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1").add_metabolites({Rubisco_balance:1})
  
  
  temp.reactions.get_by_id("SUCROSE_PROTON_vc1").pop(temp.metabolites.get_by_id("PROTON_v1"))
  temp.reactions.get_by_id("SUCROSE_PROTON_vc2").pop(temp.metabolites.get_by_id("PROTON_v2"))
  temp.reactions.get_by_id("SUCROSE_PROTON_vc1").pop(temp.metabolites.get_by_id("PROTON_c1"))
  temp.reactions.get_by_id("SUCROSE_PROTON_vc2").pop(temp.metabolites.get_by_id("PROTON_c2"))
  
  solution=flux_analysis.parsimonious.optimize_minimal_flux(temp)
  
  test_model=temp.copy()
  #temp_model.reactions.get_by_id("diel_biomass").lower_bound=0.64
  #temp_model.reactions.get_by_id("diel_biomass").upper_bound=0.64
  #solution=flux_analysis.parsimonious.optimize_minimal_flux(temp_model)
  
  print "Storing results from "+str(b)
  ansB[b]=test_model.solution.f
  solutiondict[b]=test_model.solution.x_dict
  #break



writeAllFluxes(CAM_model_pH,solutiondict,"PEPCK_fructan_passiveSuc_activePhlo.csv")



###############fructan storing - active Sucrose vacuolar loading - symplastic Phloem loading
test_model=CAM_model_pH.copy()
test_model.reactions.get_by_id("MALIC_NADP_RXN_c1").lower_bound=0
test_model.reactions.get_by_id("MALIC_NADP_RXN_c1").upper_bound=0
test_model.reactions.get_by_id("MALIC_NADP_RXN_p1").lower_bound=0
test_model.reactions.get_by_id("MALIC_NADP_RXN_p1").upper_bound=0
test_model.reactions.get_by_id("MALIC_NAD_RXN_m1").lower_bound=0
test_model.reactions.get_by_id("MALIC_NAD_RXN_m1").upper_bound=0
test_model.reactions.get_by_id("STARCH_p_dielTransfer").lower_bound = 0
test_model.reactions.get_by_id("STARCH_p_dielTransfer").upper_bound = 0
test_model.reactions.get_by_id("SUCROSE_v_dielTransfer").lower_bound = 0
test_model.reactions.get_by_id("SUCROSE_v_dielTransfer").upper_bound = 0
test_model.reactions.get_by_id("FRUCTAN_v_dielTransfer").lower_bound= -1000
test_model.reactions.get_by_id("FRUCTAN_v_dielTransfer").upper_bound= 1000

test_model.reactions.get_by_id("Phloem_output_tx1").pop(test_model.metabolites.get_by_id("PROTON_e1"))
test_model.reactions.get_by_id("Phloem_output_tx2").pop(test_model.metabolites.get_by_id("PROTON_e2"))
test_model.reactions.get_by_id("Phloem_output_tx1").pop(test_model.metabolites.get_by_id("PROTON_c1"))
test_model.reactions.get_by_id("Phloem_output_tx2").pop(test_model.metabolites.get_by_id("PROTON_c2"))

solution=flux_analysis.parsimonious.optimize_minimal_flux(test_model)

ansC = dict()
ansC[10000]=test_model.solution.f

solutiondict=dict()
solutiondict[10000]=test_model.solution.x_dict


for i in [1000,160,80,40,20,10,5,3,2]:
  b = float(i)
  print b
  Rubisco_balance = Metabolite("rubisco_bal_p1", name = "Weights to balance RuBP carboxygenase oxygenase balance", compartment = "p1")
  
  temp = CAM_model_pH.copy()
  temp.reactions.get_by_id("STARCH_p_dielTransfer").lower_bound = 0
  temp.reactions.get_by_id("STARCH_p_dielTransfer").upper_bound = 0
  temp.reactions.get_by_id("SUCROSE_v_dielTransfer").lower_bound = 0
  temp.reactions.get_by_id("SUCROSE_v_dielTransfer").upper_bound = 0
  temp.reactions.get_by_id("FRUCTAN_v_dielTransfer").lower_bound= -1000
  temp.reactions.get_by_id("FRUCTAN_v_dielTransfer").upper_bound= 1000
  temp.reactions.get_by_id("RXN_961_p1").add_metabolites({Rubisco_balance:-b})
  temp.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1").add_metabolites({Rubisco_balance:1})
  
  temp.reactions.get_by_id("Phloem_output_tx1").pop(temp.metabolites.get_by_id("PROTON_e1"))
  temp.reactions.get_by_id("Phloem_output_tx2").pop(temp.metabolites.get_by_id("PROTON_e2"))
  temp.reactions.get_by_id("Phloem_output_tx1").pop(temp.metabolites.get_by_id("PROTON_c1"))
  temp.reactions.get_by_id("Phloem_output_tx2").pop(temp.metabolites.get_by_id("PROTON_c2"))
  
  solution=flux_analysis.parsimonious.optimize_minimal_flux(temp)
  
  test_model=temp.copy()
  #temp_model.reactions.get_by_id("diel_biomass").lower_bound=0.64
  #temp_model.reactions.get_by_id("diel_biomass").upper_bound=0.64
  #solution=flux_analysis.parsimonious.optimize_minimal_flux(temp_model)
  
  print "Storing results from "+str(b)
  ansC[b]=test_model.solution.f
  solutiondict[b]=test_model.solution.x_dict
  #break



writeAllFluxes(CAM_model_pH,solutiondict,"PEPCK_fructan_activeSuc_passivePhlo.csv")


###############fructan storing - passive Sucrose vacuolar loading - symplastic Phloem loading
test_model=CAM_model_pH.copy()
test_model.reactions.get_by_id("MALIC_NADP_RXN_c1").lower_bound=0
test_model.reactions.get_by_id("MALIC_NADP_RXN_c1").upper_bound=0
test_model.reactions.get_by_id("MALIC_NADP_RXN_p1").lower_bound=0
test_model.reactions.get_by_id("MALIC_NADP_RXN_p1").upper_bound=0
test_model.reactions.get_by_id("MALIC_NAD_RXN_m1").lower_bound=0
test_model.reactions.get_by_id("MALIC_NAD_RXN_m1").upper_bound=0
test_model.reactions.get_by_id("STARCH_p_dielTransfer").lower_bound = 0
test_model.reactions.get_by_id("STARCH_p_dielTransfer").upper_bound = 0
test_model.reactions.get_by_id("SUCROSE_v_dielTransfer").lower_bound = 0
test_model.reactions.get_by_id("SUCROSE_v_dielTransfer").upper_bound = 0
test_model.reactions.get_by_id("FRUCTAN_v_dielTransfer").lower_bound= -1000
test_model.reactions.get_by_id("FRUCTAN_v_dielTransfer").upper_bound= 1000

test_model.reactions.get_by_id("SUCROSE_PROTON_vc1").pop(test_model.metabolites.get_by_id("PROTON_v1"))
test_model.reactions.get_by_id("SUCROSE_PROTON_vc2").pop(test_model.metabolites.get_by_id("PROTON_v2"))
test_model.reactions.get_by_id("SUCROSE_PROTON_vc1").pop(test_model.metabolites.get_by_id("PROTON_c1"))
test_model.reactions.get_by_id("SUCROSE_PROTON_vc2").pop(test_model.metabolites.get_by_id("PROTON_c2"))
test_model.reactions.get_by_id("Phloem_output_tx1").pop(test_model.metabolites.get_by_id("PROTON_e1"))
test_model.reactions.get_by_id("Phloem_output_tx2").pop(test_model.metabolites.get_by_id("PROTON_e2"))
test_model.reactions.get_by_id("Phloem_output_tx1").pop(test_model.metabolites.get_by_id("PROTON_c1"))
test_model.reactions.get_by_id("Phloem_output_tx2").pop(test_model.metabolites.get_by_id("PROTON_c2"))

solution=flux_analysis.parsimonious.optimize_minimal_flux(test_model)

ansD = dict()
ansD[10000]=test_model.solution.f

solutiondict=dict()
solutiondict[10000]=test_model.solution.x_dict


for i in [1000,160,80,40,20,10,5,3,2]:
  b = float(i)
  print b
  Rubisco_balance = Metabolite("rubisco_bal_p1", name = "Weights to balance RuBP carboxygenase oxygenase balance", compartment = "p1")
  
  temp = CAM_model_pH.copy()
  temp.reactions.get_by_id("STARCH_p_dielTransfer").lower_bound = 0
  temp.reactions.get_by_id("STARCH_p_dielTransfer").upper_bound = 0
  temp.reactions.get_by_id("SUCROSE_v_dielTransfer").lower_bound = 0
  temp.reactions.get_by_id("SUCROSE_v_dielTransfer").upper_bound = 0
  temp.reactions.get_by_id("FRUCTAN_v_dielTransfer").lower_bound= -1000
  temp.reactions.get_by_id("FRUCTAN_v_dielTransfer").upper_bound= 1000
  temp.reactions.get_by_id("RXN_961_p1").add_metabolites({Rubisco_balance:-b})
  temp.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1").add_metabolites({Rubisco_balance:1})
  
  temp.reactions.get_by_id("SUCROSE_PROTON_vc1").pop(temp.metabolites.get_by_id("PROTON_v1"))
  temp.reactions.get_by_id("SUCROSE_PROTON_vc2").pop(temp.metabolites.get_by_id("PROTON_v2"))
  temp.reactions.get_by_id("SUCROSE_PROTON_vc1").pop(temp.metabolites.get_by_id("PROTON_c1"))
  temp.reactions.get_by_id("SUCROSE_PROTON_vc2").pop(temp.metabolites.get_by_id("PROTON_c2"))
  temp.reactions.get_by_id("Phloem_output_tx1").pop(temp.metabolites.get_by_id("PROTON_e1"))
  temp.reactions.get_by_id("Phloem_output_tx2").pop(temp.metabolites.get_by_id("PROTON_e2"))
  temp.reactions.get_by_id("Phloem_output_tx1").pop(temp.metabolites.get_by_id("PROTON_c1"))
  temp.reactions.get_by_id("Phloem_output_tx2").pop(temp.metabolites.get_by_id("PROTON_c2"))
  
  solution=flux_analysis.parsimonious.optimize_minimal_flux(temp)
  
  test_model=temp.copy()
  #temp_model.reactions.get_by_id("diel_biomass").lower_bound=0.64
  #temp_model.reactions.get_by_id("diel_biomass").upper_bound=0.64
  #solution=flux_analysis.parsimonious.optimize_minimal_flux(temp_model)
  
  print "Storing results from "+str(b)
  ansD[b]=test_model.solution.f
  solutiondict[b]=test_model.solution.x_dict
  #break



writeAllFluxes(CAM_model_pH,solutiondict,"PEPCK_fructan_passiveSuc_passivePhlo.csv")


###########################Plot metabolic productivities expressed as % of starch-storing C3#####################

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
for k in sorted(ansA.keys()):
  if k==10000:
    xvalues.append(r'$\infty$')
  xvalues.append(int(k))
  xvalues2.append(np.log(k))
  yvalues.append(((ansA.get(k))/temp_model3.solution.f)*100)
  yvalues2.append(((ansB.get(k))/temp_model3.solution.f)*100)
  yvalues3.append(((ansC.get(k))/temp_model3.solution.f)*100)
  yvalues4.append(((ansD.get(k))/temp_model3.solution.f)*100)

curve = plt.subplot(111)
curve.plot(xvalues2,yvalues,marker="o",label="active vacuolar sucrose uptake - apoplastic phloem loading ",color="blue")
curve.plot(xvalues2,yvalues2,marker="o",label="passive vacuolar sucrose uptake - appoplastic phloem loading ",color="red")
curve.plot(xvalues2,yvalues3,marker="o",label="active vacuolar sucrose uptake - symplastic phloem loading ",color="green")
curve.plot(xvalues2,yvalues4,marker="o",label="passive vacuolar sucrose uptake - symplastic phloem loading ",color="yellow")
curve.axvspan(np.log(float(1)/0.33),np.log(float(1)/0.14),alpha=0.3,color="green")
curve.axhline(y=100,alpha=0.3,color="grey")
#curve.set_xlim(1,100)
curve.set_xlabel(r'$Vc/Vo$',fontsize=20)
curve.set_xticks(xvalues2)
curve.set_xticklabels(xvalues)
curve.set_ylabel("Phloem output as % of C"+r'$_3$'+" output",fontsize=20)


#plt.xscale('log')
plt.legend(loc="best")
plt.show()



