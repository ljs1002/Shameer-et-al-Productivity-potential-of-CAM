#########################################################################
#This script can be used to model chlorenchyma metabolism Opuntia ficus	#
#indica	 								#
#									#
#########################################################################


########################## FUNCTIONS ####################################

#This function can be used to convert a model with reversible reactions to a model with only irreversible reactions
def rev2irrev(cobra_model):
  exp_model=cobra_model.copy()
  
  for RXN in cobra_model.reactions:
    rxn=exp_model.reactions.get_by_id(RXN.id)
    if (rxn.lower_bound < 0 and rxn.upper_bound > 0):
      rxn_reverse = rxn.copy()
      rxn_reverse.id = "%s_reverse" %(rxn.id)
      rxn.lower_bound = 0
      rxn_reverse.upper_bound = 0
      exp_model.add_reaction(rxn_reverse)
  
  return exp_model


#This function updates phloem composition of a model from a python dictionary
def updatePhloemComposition(model,PhCompF):
  newPhloem = dict()
  fin = open(PhCompF,"r")
  for line in fin:
    line=line.replace("\n","")
    lineparts = line.split("\t")
    newPhloem[lineparts[0]]=float(lineparts[1])
  for i in [1,2]:
    temp = list()
    for met in model.reactions.get_by_id("Phloem_output_tx"+str(i)).metabolites:
      if(not met.id=="X_Phloem_contribution_t"+str(i)):
	temp.append(met)
    for met in temp:
      model.reactions.get_by_id("Phloem_output_tx"+str(i)).pop(met)
    for met in newPhloem.keys():
      model.reactions.get_by_id("Phloem_output_tx"+str(i)).add_metabolites({model.metabolites.get_by_id(met+str(i)):-1*newPhloem.get(met)})
    return model
  

#this function converts a model which fractionally charged metabolites to a classical model
def convertToClassicalModel(cobra_model2,comp="*",updateCharges=""):
  
  if not updateCharges == "":
    #override charges from file
    fin = open(updateCharges,"r")
    
    ChargeDict=dict()
    for line in fin:
      met=line.replace("\n","").split("\t")[0]
      met = met.replace("-","_")
      charge = line.replace("\n","").split("\t")[1]
      ChargeDict[met]=charge
    
    fin.close()
    
    temp=cobra_model2.copy()
    for met in temp.metabolites:
      if(met.compartment=="b"):
	cobra_model2.metabolites.get_by_id(met.id).remove_from_model()
    
    for met in cobra_model2.metabolites:
      tempMet=met.id
      if(met.id[len(met.id)-1]=="2" or met.id[len(met.id)-1]=="1"):
	tempMet = met.id[0:len(met.id)-1]
      if(ChargeDict.keys().__contains__(tempMet)):
	met.charge = ChargeDict.get(tempMet)
      if met.charge is None:
	met.charge=0
    
  
  import re
  #List metabolites that were added
  fractionMets=set()
  for rxn in cobra_model2.reactions:
    for met in rxn.metabolites.keys():
      a=re.search("^a{1,3}",met.id)
      anion=""
      if a:
	anion=a.group(0)
      b=re.search("^b{1,3}",met.id)
      basic=""
      if b:
	basic=b.group(0)
      if (abs(rxn.metabolites.get(met)) % 1 > 0 and (not anion == "" or not basic== "")):
	if not comp == "*":
	  if met.compartment == comp:
	    fractionMets.add(met.id)
	else:
	  fractionMets.add(met.id)
  
  
  
  cobra_model = cobra_model2.copy()
  
  
  for met in fractionMets:
    for rxn in cobra_model.metabolites.get_by_id(met).reactions:
      a=re.search("^a{1,3}",met)
      anion=""
      if a:
	anion=a.group(0)
      b=re.search("^b{1,3}",met)
      basic=""
      if b:
	basic=b.group(0)
      prefix = anion
      if prefix == "":
	prefix = basic
      if(rxn.id=="Protein_Processing_c1"):
	print rxn.reaction
      coeff1 = rxn.pop(cobra_model.metabolites.get_by_id(met[len(prefix):]))
      #print rxn.reaction
      coeff2 = rxn.pop(cobra_model.metabolites.get_by_id(met))
      #print rxn.reaction
      Charge=(coeff1*float(cobra_model.metabolites.get_by_id(met[len(prefix):]).charge))+(coeff2*float(cobra_model.metabolites.get_by_id(met).charge))-((coeff1+coeff2)*float(cobra_model.metabolites.get_by_id(met[len(prefix):]).charge))
      if rxn.id.__contains__("_dielTransfer"):
	rxn.add_metabolites({cobra_model.metabolites.get_by_id(met[len(prefix):]):coeff1+coeff2,cobra_model.metabolites.get_by_id("PROTON_"+cobra_model.metabolites.get_by_id(met).compartment[0:len(cobra_model.metabolites.get_by_id(met).compartment)-1]+"1"):round(Charge,4)})
      else:
	rxn.add_metabolites({cobra_model.metabolites.get_by_id(met[len(prefix):]):coeff1+coeff2,cobra_model.metabolites.get_by_id("PROTON_"+cobra_model.metabolites.get_by_id(met).compartment):round(Charge,4)})
      if(rxn.id=="Protein_Processing_c1"):
	print Charge
	print rxn.reaction
  
  
  uncharged = cobra_model.copy()
  return uncharged




#this function converts a classical model  a model with fractional charges, a file were fractional charges are recorded is required
def convertToFractionalCharges(uncharged,filename,compH={"v1":5.5,"v2":3.3}):
  
  fin = open(filename,"r")
  
  i=0
  FractionDict=dict()
  FractionCharges=dict()
  for line in fin:
    if i==0:
      i=i+1
      continue
    lineparts = line.replace("\n","").split("\t")
    masterMet = lineparts[0]
    numberMet = lineparts[1]
    #print numberMet
    chargeFraction = dict()
    for j in range(1,int(numberMet)+1):
      chargeFraction[lineparts[(3*int(j))-1]]=lineparts[(3*int(j))+1]
      #print lineparts[(3*int(j))-1]
      #print lineparts[(3*int(j))]
      FractionCharges[lineparts[(3*int(j))-1]]=lineparts[(3*int(j))]
    FractionDict[masterMet]=chargeFraction
    #print chargeFraction
  
  
  fin.close()
  
  
  charged = uncharged.copy()
  #uncharged=backup.copy()
  for met in FractionDict.keys():
    #print met
    for rxn in uncharged.metabolites.get_by_id(met).reactions:
      #print met
      coeff = charged.reactions.get_by_id(rxn.id).pop(charged.metabolites.get_by_id(met))
      for tmet in FractionDict.get(met).keys():
	if not uncharged.metabolites.has_id(tmet):
	  metab = charged.metabolites.get_by_id(met).copy()
	  metab.id = tmet
	  metab.charge = int(FractionCharges.get(tmet))
	else:
	  metab = charged.metabolites.get_by_id(tmet)
	if not FractionDict.get(met).get(tmet) == "0":
	  charged.reactions.get_by_id(rxn.id).add_metabolites({metab:coeff*float(FractionDict.get(met).get(tmet))})
	#rxn.reaction
      Charge=0
      for tmet in charged.reactions.get_by_id(rxn.id).metabolites.keys():
	Charge = Charge + (int(tmet.charge) * charged.reactions.get_by_id(rxn.id).metabolites.get(tmet))
      
      
      if not Charge == 0:
	#print rxn.id+"\t"+str(Charge)
	#break
	if rxn.id.__contains__("_dielTransfer"):
	  compList=list()
	  for Rmet in rxn.metabolites:
	    compList.append(Rmet.compartment)
	  compSet = set(compList)
	  compList = list(compSet)
	  for comp in compList:
	    if charged.reactions.get_by_id(rxn.id).metabolites.has_key(charged.metabolites.get_by_id("PROTON_"+comp)):
	      #print charged.reactions.get_by_id(rxn.id).reaction
	      Charge = Charge - charged.reactions.get_by_id(rxn.id).pop(charged.metabolites.get_by_id("PROTON_"+comp))
	      #print charged.reactions.get_by_id(rxn.id).reaction
		  
	  if compH.get(compList[0]) < compH.get(compList[1]):
	    protSource = compList[1]
	    protSink = compList[0]
	  else:
	    protSource = compList[0]
	    protSink = compList[1]
	  if Charge > 0:
	    t = protSink
	    protSink = protSource
	    protSource = t
	  
	  if rxn.lower_bound < 0 and rxn.upper_bound > 0:
	    rxn_reverse = charged.reactions.get_by_id(rxn.id).copy()
	    rxn_reverse.id = rxn_reverse.id+"_R"
	    
	    if Charge > 0:
	      charged.reactions.get_by_id(rxn.id).add_metabolites({charged.metabolites.get_by_id("PROTON_"+protSource):-Charge})
	      charged.reactions.get_by_id(rxn.id).lower_bound = 0
	      rxn_reverse.add_metabolites({charged.metabolites.get_by_id("PROTON_"+protSink):-Charge})
	      rxn_reverse.upper_bound=0
	      charged.add_reactions({rxn_reverse})
	    else:
	      charged.reactions.get_by_id(rxn.id).add_metabolites({charged.metabolites.get_by_id("PROTON_"+protSink):-Charge})
	      charged.reactions.get_by_id(rxn.id).lower_bound = 0
	      rxn_reverse.add_metabolites({charged.metabolites.get_by_id("PROTON_"+protSource):-Charge})
	      rxn_reverse.upper_bound=0
	  elif rxn.lower_bound < 0 and rxn.upper_bound == 0:
	    if (Charge > 0):
	      charged.reactions.get_by_id(rxn.id).add_metabolites({charged.metabolites.get_by_id("PROTON_"+protSink):-Charge})
	    else:
	      charged.reactions.get_by_id(rxn.id).add_metabolites({charged.metabolites.get_by_id("PROTON_"+protSource):-Charge})
	  elif rxn.lower_bound == 0 and rxn.upper_bound > 0:
	    if (Charge > 0):
	      charged.reactions.get_by_id(rxn.id).add_metabolites({charged.metabolites.get_by_id("PROTON_"+protSource):-Charge})
	    else:
	      charged.reactions.get_by_id(rxn.id).add_metabolites({charged.metabolites.get_by_id("PROTON_"+protSink):-Charge})
	  #print charged.reactions.get_by_id(rxn.id).reaction
	  #print rxn_reverse.reaction
	else:
	  charged.reactions.get_by_id(rxn.id).add_metabolites({charged.metabolites.get_by_id("PROTON_"+charged.metabolites.get_by_id(met).compartment):-Charge})
	  #print rxn.reaction 
  return charged






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

#transform phloem composition to that reported in Opuntia ficus indica for validation
CAM_model_pH = updatePhloemComposition(CAM_model_pH,"OpuntiaBiomass.csv")
CAM_model_pH.reactions.get_by_id("Photon_tx1").lower_bound=0
CAM_model_pH.reactions.get_by_id("Photon_tx1").upper_bound=1000

#Opuntia CO2 assimilation vs PPFD curve - Nobel and Hartsock 1983
CO2_data_opuntia = {200:59.2,300:148.3,400:247.9,500:296.5,600:323.8,700:336.0,800:338.5}

#Constraint CO2 uptake
CAM_model_pH.reactions.get_by_id("CO2_tx2").lower_bound=float(296.5*1000)/(12*60*60)
CAM_model_pH.reactions.get_by_id("CO2_tx2").upper_bound=float(296.5*1000)/(12*60*60)

#Adding linear maintenance constraint ATPase =(0.0423454545455 * photon) + 0.06
met1 = Metabolite("Maintenance_estimator1",name="metabolite to represent ATPase maintence flux",compartment="c1")
rxn = Reaction("LR_intercet")
rxn.add_metabolites({met1:1})
rxn.lower_bound=0.06
rxn.upper_bound=0.06
CAM_model_pH.add_reaction(rxn)

met2 = Metabolite("Maintenance_estimator2",name="metaboltie to represent ATPase maintenance as a function of photon flux", compartment="c1")
rxn = Reaction("LR_slope")
rxn.add_metabolites({met2:-1,met1:1})
rxn.lower_bounds=0
rxn.upper_bound=1000
CAM_model_pH.add_reaction(rxn)

CAM_model_pH.reactions.get_by_id("ATPase_tx2").add_metabolites({met1:-1})
CAM_model_pH.reactions.get_by_id("Photon_tx1").add_metabolites({met2:0.0423454545455})

#perform parsimonious optimization
solution=flux_analysis.parsimonious.optimize_minimal_flux(CAM_model_pH)

#check starch and malate accumulation
CAM_model_pH.reactions.get_by_id("MAL_v_dielTransfer_R").x * (float(12*60*60)/1000)
CAM_model_pH.reactions.get_by_id("STARCH_p_dielTransfer").x* (float(12*60*60)/1000)