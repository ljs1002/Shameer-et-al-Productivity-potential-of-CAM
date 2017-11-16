#########################################################################
#This script can be used to model mature C3 and CAM leaves with fixed	#
#producitivty								#
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

fin = open("FractionalCharges.csv","r")


#update metabolite charges
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

#constraint VcVo
Rubisco_balance = Metabolite("rubisco_bal_p1", name = "Weights to balance RuBP carboxygenase oxygenase balance", compartment = "p1")
CAM_model.reactions.get_by_id("RXN_961_p1").add_metabolites({Rubisco_balance:5.15})
CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1").add_metabolites({Rubisco_balance:-1})
CAM_model.reactions.get_by_id("RXN_961_p2").add_metabolites({Rubisco_balance:5.15})
CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2").add_metabolites({Rubisco_balance:-1})


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

#constraining light to 200 and maintenance ATPase flux to 8.5 based on result from script1
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

#perform parsimonious optimization of starch hydrolyzing CAM - maltose released from plastid
temp_model_1=CAM_model_pH.copy()
temp_model_1.reactions.get_by_id("RXN_1826_p2").lower_bound=0
temp_model_1.reactions.get_by_id("RXN_1826_p2").upper_bound=0
temp_model_1.reactions.get_by_id("MALTODEXGLUCOSID_RXN_p2").lower_bound=0
temp_model_1.reactions.get_by_id("MALTODEXGLUCOSID_RXN_p2").upper_bound=0
temp_model_1.reactions.get_by_id("diel_biomass").lower_bound=0.178
temp_model_1.reactions.get_by_id("diel_biomass").upper_bound=0.178
solution=flux_analysis.parsimonious.optimize_minimal_flux(temp_model_1)

#perform parsimonious optimization of starch hydrolyzing CAM - glucose released from plastid
temp_model_2=CAM_model_pH.copy()
temp_model_2.reactions.get_by_id("RXN_1826_p2").lower_bound=0
temp_model_2.reactions.get_by_id("RXN_1826_p2").upper_bound=0
temp_model_2.reactions.get_by_id("RXN_1827_p2").lower_bound=0
temp_model_2.reactions.get_by_id("RXN_1827_p2").upper_bound=0
temp_model_2.reactions.get_by_id("diel_biomass").lower_bound=0.178
temp_model_2.reactions.get_by_id("diel_biomass").upper_bound=0.178
solution=flux_analysis.parsimonious.optimize_minimal_flux(temp_model_2)


#perform parsimonious optimization of starch phosphorylating CAM
temp_model2=CAM_model_pH.copy()
temp_model2.reactions.get_by_id("MALTODEXGLUCOSID_RXN_p2").lower_bound=0
temp_model2.reactions.get_by_id("MALTODEXGLUCOSID_RXN_p2").upper_bound=0
temp_model2.reactions.get_by_id("RXN_1827_p2").lower_bound=0
temp_model2.reactions.get_by_id("RXN_1827_p2").upper_bound=0
temp_model2.reactions.get_by_id("diel_biomass").lower_bound=0.178
temp_model2.reactions.get_by_id("diel_biomass").upper_bound=0.178
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


#nitrate uptake constraint
Nitrate_balance = Metabolite("Nitrate_bal_c", name = "Weights to balance nitrate uptake", compartment = "c1")
C3_model.reactions.get_by_id("Nitrate_ec1").add_metabolites({Nitrate_balance:-2})
C3_model.reactions.get_by_id("Nitrate_ec2").add_metabolites({Nitrate_balance:3})


#constraint VcVo
Rubisco_balance = Metabolite("rubisco_bal_p1", name = "Weights to balance RuBP carboxygenase oxygenase balance", compartment = "p1")
C3_model.reactions.get_by_id("RXN_961_p1").add_metabolites({Rubisco_balance:3})#changed from Alison Smith's 3 to Leegood's 2.5 to back to back to 3 based on Ma et al 2014 ()
C3_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1").add_metabolites({Rubisco_balance:-1})

#constraint ATPase:NADPHoxidase = 3:1 and also ATPase (day)=ATPase(night)
Maintenance_constraint = Metabolite("ATPase_NADPHoxidase_constraint_c1",name =  "ATPase_NADPHoxidase_constraint_c1", compartment = "c1")
Maintenance_constraint2 = Metabolite("ATPase_NADPHoxidase_constraint_c2",name =  "ATPase_NADPHoxidase_constraint_c2", compartment = "c2")
Maintenance_constraint3 = Metabolite("Light_dark_maintainence_constraint",name =  "Light_dark_maintainence_constraint", compartment = "c1")
C3_model.reactions.get_by_id("ATPase_tx1").add_metabolites({Maintenance_constraint:1,Maintenance_constraint3:1})
C3_model.reactions.get_by_id("ATPase_tx2").add_metabolites({Maintenance_constraint2:1,Maintenance_constraint3:-1})
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

#perform parsimonious optimization
solution=flux_analysis.parsimonious.optimize_minimal_flux(C3_model)

#perform parsimonious optimization of starch hydrolyzing C3 - maltose released from plastid
temp_model3_1 = C3_model.copy()
temp_model3_1.reactions.get_by_id("RXN_1826_p2").lower_bound=0
temp_model3_1.reactions.get_by_id("RXN_1826_p2").upper_bound=0
temp_model3_1.reactions.get_by_id("MALTODEXGLUCOSID_RXN_p2").lower_bound=0
temp_model3_1.reactions.get_by_id("MALTODEXGLUCOSID_RXN_p2").upper_bound=0
temp_model3_1.reactions.get_by_id("diel_biomass").lower_bound=0.178
temp_model3_1.reactions.get_by_id("diel_biomass").upper_bound=0.178
solution=flux_analysis.parsimonious.optimize_minimal_flux(temp_model3_1)

#perform parsimonious optimization of starch hydrolyzing C3 - glucose released from plastid
temp_model3_2 = C3_model.copy()
temp_model3_2.reactions.get_by_id("RXN_1826_p2").lower_bound=0
temp_model3_2.reactions.get_by_id("RXN_1826_p2").upper_bound=0
temp_model3_2.reactions.get_by_id("RXN_1827_p2").lower_bound=0
temp_model3_2.reactions.get_by_id("RXN_1827_p2").upper_bound=0
temp_model3_2.reactions.get_by_id("diel_biomass").lower_bound=0.178
temp_model3_2.reactions.get_by_id("diel_biomass").upper_bound=0.178
solution=flux_analysis.parsimonious.optimize_minimal_flux(temp_model3_2)

#perform parsimonious optimization of starch phosphorylating C3
temp_model4=C3_model.copy()
temp_model4.reactions.get_by_id("MALTODEXGLUCOSID_RXN_p2").lower_bound=0
temp_model4.reactions.get_by_id("MALTODEXGLUCOSID_RXN_p2").upper_bound=0
temp_model4.reactions.get_by_id("RXN_1827_p2").lower_bound=0
temp_model4.reactions.get_by_id("RXN_1827_p2").upper_bound=0
temp_model4.reactions.get_by_id("diel_biomass").lower_bound=0.178
temp_model4.reactions.get_by_id("diel_biomass").upper_bound=0.178
solution=flux_analysis.parsimonious.optimize_minimal_flux(temp_model4)




####################### WRITING FLUXES #####################################

fin = file("rxnpathwaydict.csv","r")
pathclass = dict()
pathname=dict()
for line in fin:
  line=line.replace("\n","")
  lineparts=line.split(",")
  pathname[lineparts[0]]=lineparts[1]


fin = file("pathwayprioritizer.csv","r")
pathindex = dict()

for line in fin:
  line=line.replace("\n","")
  lineparts=line.split("\t")
  pathindex[lineparts[0]]=lineparts[1]



fout = open("CAM_C3_fluxes_fixedOutput.csv","w")
fout.write("Reaction ID\tReaction Name\tEquation\tCompartment\tPathway\tPathway index\tE.C number\tC3 optimal\tC3 starch-hydrolyzing & maltose-exporting\tC3 starch-hydrolyzing & glucose-exporting\tC3 no hydrolysis\tCAM optimal\tCAM starch-hydrolyzing & maltose-exporting\tCAM starch-hydrolyzing & glucose-exporting\tCAM no hydrolysis\n")
for rxn in C3_model.reactions:
  compSet=set()
  for met in rxn.metabolites.keys():
    compSet.add(met.compartment)
  comp="";
  for C in compSet:
    comp=comp+C
  path=""
  index=""
  EC=""
  
  RXN=rxn.id
  if(rxn.id[len(rxn.id)-1]=="1" or rxn.id[len(rxn.id)-1]=="2"):
    RXN = rxn.id[0:len(rxn.id)-1]
  
  if(pathname.keys().__contains__(RXN)):
    path=pathname.get(RXN)
  
  if pathindex.keys().__contains__(path):
    index=pathindex.get(path)
  
  if(not str(rxn.notes)=="{}"):
    EC=rxn.notes.get("PROTEIN CLASS")[0]
  
  fout.write(rxn.id+"\t"+rxn.name+"\t"+rxn.reaction+"\t"+comp+"\t"+path+"\t"+index+"\t"+EC+"\t"+str(C3_model.solution.x_dict.get(rxn.id))+"\t"+str(temp_model3_1.solution.x_dict.get(rxn.id))+"\t"+str(temp_model3_2.solution.x_dict.get(rxn.id))+"\t"+str(temp_model4.solution.x_dict.get(rxn.id))+"\t"+str(CAM_model_pH.solution.x_dict.get(rxn.id))+"\t"+str(temp_model_1.solution.x_dict.get(rxn.id))+"\t"+str(temp_model_2.solution.x_dict.get(rxn.id))+"\t"+str(temp_model2.solution.x_dict.get(rxn.id))+"\n")


for rxn in CAM_model_pH.reactions:
  if len(C3_model.reactions.query(rxn.id)) == 0:
    compSet=set()
    for met in rxn.metabolites.keys():
      compSet.add(met.compartment)
    comp="";
    for C in compSet:
      comp=comp+C
    path=""
    index=""
    EC=""
    
    RXN=rxn.id
    if(rxn.id[len(rxn.id)-1]=="1" or rxn.id[len(rxn.id)-1]=="2"):
      RXN = rxn.id[0:len(rxn.id)-1]
    
    if(pathname.keys().__contains__(RXN)):
      path=pathname.get(RXN)
    
    if pathindex.keys().__contains__(path):
      index=pathindex.get(path)
    
    if(not str(rxn.notes)=="{}"):
      EC=rxn.notes.get("PROTEIN CLASS")[0]
    
    fout.write(rxn.id+"\t"+rxn.name+"\t"+rxn.reaction+"\t"+comp+"\t"+path+"\t"+index+"\t"+EC+"\t0\t0\t0\t0\t"+str(CAM_model_pH.solution.x_dict.get(rxn.id))+"\t"+str(temp_model_1.solution.x_dict.get(rxn.id))+"\t"+str(temp_model_2.solution.x_dict.get(rxn.id))+"\t"+str(temp_model2.solution.x_dict.get(rxn.id))+"\n")
    

fout.close()


################################ TO PRINT ATP AND NADPH REACTIONS ###########################
#replace suffix 1 with 2 to print dark fluxes
fout = open("NADPH_fluxes","w")

for p in ("c","p","m"):
  met=C3_model.metabolites.get_by_id("NADPH_"+p+"1")
  for rxn in met.reactions:
    #if (rxn.metabolites.get(met)>0) or (rxn.reversibility and rxn.metabolites.get(met)<0):
      sto=rxn.metabolites.get(met)
      sto1=0#rxn.metabolites.get(met1)
      fout.write(rxn.id+"#"+rxn.reaction+"#"+str(C3_model.solution.x_dict.get(rxn.id)*(sto+sto1))+"#"+str(temp_model3_1.solution.x_dict.get(rxn.id)*(sto+sto1))+"#"+str(temp_model3_2.solution.x_dict.get(rxn.id)*(sto+sto1))+"#"+str(temp_model4.solution.x_dict.get(rxn.id)*(sto+sto1))+"#"+str(CAM_model_pH.solution.x_dict.get(rxn.id)*(sto+sto1))+"#"+str(temp_model_1.solution.x_dict.get(rxn.id)*(sto+sto1))+"#"+str(temp_model_2.solution.x_dict.get(rxn.id)*(sto+sto1))+"#"+str(temp_model2.solution.x_dict.get(rxn.id)*(sto+sto1))+"#"+met.compartment)

fout.close()

fout = open("NADH_fluxes","w")

for p in ("c","p","m","x"):
  met=C3_model.metabolites.get_by_id("NADH_"+p+"1")
  for rxn in met.reactions:
      sto=rxn.metabolites.get(met)
      sto1=0
      fout.write(rxn.id+"#"+rxn.reaction+"#"+str(C3_model.solution.x_dict.get(rxn.id)*(sto+sto1))+"#"+str(temp_model3_1.solution.x_dict.get(rxn.id)*(sto+sto1))+"#"+str(temp_model3_2.solution.x_dict.get(rxn.id)*(sto+sto1))+"#"+str(temp_model4.solution.x_dict.get(rxn.id)*(sto+sto1))+"#"+str(CAM_model_pH.solution.x_dict.get(rxn.id)*(sto+sto1))+"#"+str(temp_model_1.solution.x_dict.get(rxn.id)*(sto+sto1))+"#"+str(temp_model_2.solution.x_dict.get(rxn.id)*(sto+sto1))+"#"+str(temp_model2.solution.x_dict.get(rxn.id)*(sto+sto1))+"#"+met.compartment)

fout.close()

fout = open("ATP_fluxes","w")

for p in ("c","p","m","x"):
  met=C3_model.metabolites.get_by_id("ATP_"+p+"1")
  met1=C3_model.metabolites.get_by_id("aATP_"+p+"1")
  for rxn in met.reactions:
    #if (rxn.metabolites.get(met)>0) or (rxn.reversibility and rxn.metabolites.get(met)<0):
      sto=rxn.metabolites.get(met)
      sto1=rxn.metabolites.get(met1)
      fout.write(rxn.id+"#"+rxn.reaction+"#"+str(C3_model.solution.x_dict.get(rxn.id)*(sto+sto1))+"#"+str(temp_model3_1.solution.x_dict.get(rxn.id)*(sto+sto1))+"#"+str(temp_model3_2.solution.x_dict.get(rxn.id)*(sto+sto1))+"#"+str(temp_model4.solution.x_dict.get(rxn.id)*(sto+sto1))+"#"+str(CAM_model_pH.solution.x_dict.get(rxn.id)*(sto+sto1))+"#"+str(temp_model_1.solution.x_dict.get(rxn.id)*(sto+sto1))+"#"+str(temp_model_2.solution.x_dict.get(rxn.id)*(sto+sto1))+"#"+str(temp_model2.solution.x_dict.get(rxn.id)*(sto+sto1))+"#"+met.compartment)

fout.close()


################################ TO PRINT PROTON REACTIONS ###########################
#replace suffix 1 with 2 to print dark fluxes
fout = open("PROTON_fluxes","w")

for p in ("c","p","m","x","e"):
  met=C3_model.metabolites.get_by_id("PROTON_"+p+"1")
  for rxn in met.reactions:
    #if (rxn.metabolites.get(met)>0) or (rxn.reversibility and rxn.metabolites.get(met)<0):
      sto=rxn.metabolites.get(met)
      sto1=0#rxn.metabolites.get(met1)
      fout.write(rxn.id+"#"+rxn.reaction+"#"+str(C3_model.solution.x_dict.get(rxn.id)*(sto+sto1))+"#"+str(temp_model3_1.solution.x_dict.get(rxn.id)*(sto+sto1))+"#"+str(temp_model3_2.solution.x_dict.get(rxn.id)*(sto+sto1))+"#"+str(temp_model4.solution.x_dict.get(rxn.id)*(sto+sto1))+"#"+str(CAM_model_pH.solution.x_dict.get(rxn.id)*(sto+sto1))+"#"+str(temp_model_1.solution.x_dict.get(rxn.id)*(sto+sto1))+"#"+str(temp_model_2.solution.x_dict.get(rxn.id)*(sto+sto1))+"#"+str(temp_model2.solution.x_dict.get(rxn.id)*(sto+sto1))+"#"+met.compartment)

fout.close()
