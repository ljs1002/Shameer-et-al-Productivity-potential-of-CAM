#########################################################################
#This script can be used to model CAM in different phloem compositions	#
#across a wide range of VcVo ratios and compare the metabolic 		#
#productivity at different NGAM values					#
#									#
#########################################################################

########################## FUNCTIONS ####################################

#write output from varying VcVo ratio
def writeAllFluxes(model,solutiondict,filename):
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
  
  
  
  fout = open(filename,"w")
  fout.write("Reaction ID\tReaction Name\tEquation\tCompartment\tPathway\tPathway index\tE.C number\tVc/Vo=inf\tVc/Vo=1000\tVc/Vo=160\tVc/Vo=80\tVc/Vo=40\tVc/Vo=20\tVc/Vo=10\tVc/Vo=5\tVc/Vo=3\tVc/Vo=2\n")
  for rxn in model.reactions:
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
    
    fout.write(rxn.id+"\t"+rxn.name+"\t"+rxn.reaction+"\t"+comp+"\t"+path+"\t"+index+"\t"+EC)
    for k in (10000,1000,160,80,40,20,10,5,3,2):
      fout.write("\t"+str(solutiondict.get(k).get(rxn.id)))
    fout.write("\n")
    
  fout.close()
      


#a function to update the phloem composition
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


#nitrate uptake constraint
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

#starch phosphorylating
CAM_model_pH.reactions.get_by_id("MALTODEXGLUCOSID_RXN_p2").lower_bound=0
CAM_model_pH.reactions.get_by_id("MALTODEXGLUCOSID_RXN_p2").upper_bound=0
CAM_model_pH.reactions.get_by_id("RXN_1827_p2").lower_bound=0
CAM_model_pH.reactions.get_by_id("RXN_1827_p2").upper_bound=0

#31
CAM_model_pH_31 = CAM_model_pH.copy()
CAM_model_pH_31.reactions.get_by_id("ATPase_tx2").lower_bound=8.5
CAM_model_pH_31.reactions.get_by_id("ATPase_tx2").upper_bound=8.5
solution=flux_analysis.parsimonious.optimize_minimal_flux(CAM_model_pH_31)

#32
CAM_model_pH_32 = CAM_model_pH.copy()
CAM_model_pH_32.reactions.get_by_id("ATPase_tx2").lower_bound=11.5
CAM_model_pH_32.reactions.get_by_id("ATPase_tx2").upper_bound=11.5

rxn = CAM_model_pH_32.reactions.get_by_id("diel_biomass")
rxn.add_metabolites({CAM_model_pH_32.metabolites.get_by_id("X_Phloem_contribution_t2"):-1})

solution=flux_analysis.parsimonious.optimize_minimal_flux(CAM_model_pH_32)

#11
CAM_model_pH_11 = CAM_model_pH.copy()
CAM_model_pH_11.reactions.get_by_id("ATPase_tx2").lower_bound=13.0
CAM_model_pH_11.reactions.get_by_id("ATPase_tx2").upper_bound=13.0

rxn = CAM_model_pH_11.reactions.get_by_id("diel_biomass")
rxn.add_metabolites({CAM_model_pH_11.metabolites.get_by_id("X_Phloem_contribution_t1"):2})

solution=flux_analysis.parsimonious.optimize_minimal_flux(CAM_model_pH_11)



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

#Set G6P transporter and Phosphorylase to 0
C3_model.reactions.get_by_id("G6P_Pi_pc1").lower_bound=0
C3_model.reactions.get_by_id("G6P_Pi_pc1").upper_bound=0
C3_model.reactions.get_by_id("G6P_Pi_pc2").lower_bound=0
C3_model.reactions.get_by_id("G6P_Pi_pc2").upper_bound=0

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

#starch hydrolyzing
C3_model.reactions.get_by_id("RXN_1826_p2").lower_bound=0
C3_model.reactions.get_by_id("RXN_1826_p2").upper_bound=0

#31
C3_model_31 = C3_model.copy()
C3_model_31.reactions.get_by_id("ATPase_tx2").lower_bound=8.5
C3_model_31.reactions.get_by_id("ATPase_tx2").upper_bound=8.5
solution=flux_analysis.parsimonious.optimize_minimal_flux(C3_model_31)

#32
C3_model_32 = C3_model.copy()
C3_model_32.reactions.get_by_id("ATPase_tx2").lower_bound=11.5
C3_model_32.reactions.get_by_id("ATPase_tx2").upper_bound=11.5

rxn = C3_model_32.reactions.get_by_id("diel_biomass")
rxn.add_metabolites({C3_model_32.metabolites.get_by_id("X_Phloem_contribution_t2"):-1})

solution=flux_analysis.parsimonious.optimize_minimal_flux(C3_model_32)

#11
C3_model_11 = C3_model.copy()
C3_model_11.reactions.get_by_id("ATPase_tx2").lower_bound=13.0
C3_model_11.reactions.get_by_id("ATPase_tx2").upper_bound=13.0

rxn = C3_model_11.reactions.get_by_id("diel_biomass")
rxn.add_metabolites({C3_model_11.metabolites.get_by_id("X_Phloem_contribution_t1"):2})

solution=flux_analysis.parsimonious.optimize_minimal_flux(C3_model_11)




################################### Vc/Vo CURVE FOR CAM WITH DIFFERENT PHLOEM COMPOSITIONS #############

#3:1 day-night phloem export
test_model=CAM_model_pH_31.copy()
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


ansA = dict()
ansA[10000]=(test_model.reactions.get_by_id("Phloem_output_tx1").x+test_model.reactions.get_by_id("Phloem_output_tx2").x)*12*60*60/1000
solutiondict=dict()
solutiondict[10000]=test_model.solution.x_dict
#print ans[0]
for i in [1000,160,80,40,20,10,5,3,2]:
  b = float(i)
  print b
  Rubisco_balance = Metabolite("rubisco_bal_p1", name = "Weights to balance RuBP carboxygenase oxygenase balance", compartment = "p1")
  
  temp = CAM_model_pH_31.copy()
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
  ansA[b]=(test_model.reactions.get_by_id("Phloem_output_tx1").x+test_model.reactions.get_by_id("Phloem_output_tx2").x)*12*60*60/1000
  solutiondict[b]=test_model.solution.x_dict
  #break


writeAllFluxes(CAM_model_pH_31,solutiondict,"PEPCK_starch_31.csv")

#########################################Vc/Vo curve for different CAM subtypes###

#3:2 day-night phloem export
test_model=CAM_model_pH_32.copy()
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


ansA_1 = dict()
ansA_1[10000]=(test_model.reactions.get_by_id("Phloem_output_tx1").x+test_model.reactions.get_by_id("Phloem_output_tx2").x)*12*60*60/1000
solutiondict=dict()
solutiondict[10000]=test_model.solution.x_dict
#print ans[0]
for i in [1000,160,80,40,20,10,5,3,2]:
  b = float(i)
  print b
  Rubisco_balance = Metabolite("rubisco_bal_p1", name = "Weights to balance RuBP carboxygenase oxygenase balance", compartment = "p1")
  
  temp = CAM_model_pH_32.copy()
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
  ansA_1[b]=(test_model.reactions.get_by_id("Phloem_output_tx1").x+test_model.reactions.get_by_id("Phloem_output_tx2").x)*12*60*60/1000
  solutiondict[b]=test_model.solution.x_dict
  #break


writeAllFluxes(CAM_model_pH_32,solutiondict,"PEPCK_starch_32.csv")

#########################################Vc/Vo curve for different CAM subtypes###

#1:1 day-night phloem export
test_model=CAM_model_pH_11.copy()
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


ansA_2 = dict()
ansA_2[10000]=(test_model.reactions.get_by_id("Phloem_output_tx1").x+test_model.reactions.get_by_id("Phloem_output_tx2").x)*12*60*60/1000
solutiondict=dict()
solutiondict[10000]=test_model.solution.x_dict
#print ans[0]
for i in [1000,160,80,40,20,10,5,3,2]:
  b = float(i)
  print b
  Rubisco_balance = Metabolite("rubisco_bal_p1", name = "Weights to balance RuBP carboxygenase oxygenase balance", compartment = "p1")
  
  temp = CAM_model_pH_11.copy()
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
  ansA_2[b]=(test_model.reactions.get_by_id("Phloem_output_tx1").x+test_model.reactions.get_by_id("Phloem_output_tx2").x)*12*60*60/1000
  solutiondict[b]=test_model.solution.x_dict
  #break


writeAllFluxes(CAM_model_pH_11,solutiondict,"PEPCK_starch_11.csv")

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
for k in sorted(ansA.keys()):
  if k==10000:
    xvalues.append(r'$\infty$')
  xvalues.append(int(k))
  xvalues2.append(np.log(k))
  yvalues.append(ansA.get(k))
  yvalues2.append(ansA_1.get(k))
  yvalues3.append(ansA_2.get(k))

curve = plt.subplot(111)
curve.plot(xvalues2,yvalues,marker="o",label="3:1 CAM",color="green")
curve.plot(xvalues2,yvalues2,marker="o",label="3:2 CAM",color="red")
curve.plot(xvalues2,yvalues3,marker="o",label="1:1 CAM",color="blue")
curve.axvspan(np.log(float(1)/0.33),np.log(float(1)/0.14),alpha=0.3,color="green")

curve.plot(xvalues2[1],(C3_model_31.reactions.get_by_id("Phloem_output_tx1").x+C3_model_31.reactions.get_by_id("Phloem_output_tx2").x)*12*60*60/1000,marker="s",markersize=10,color="green")
curve.axhline(y=(C3_model_31.reactions.get_by_id("Phloem_output_tx1").x+C3_model_31.reactions.get_by_id("Phloem_output_tx2").x)*12*60*60/1000,ls="--",alpha=0.3,color="green")
curve.plot(xvalues2[1],(C3_model_32.reactions.get_by_id("Phloem_output_tx1").x+C3_model_32.reactions.get_by_id("Phloem_output_tx2").x)*12*60*60/1000,marker="s",markersize=10,color="red")
curve.axhline(y=(C3_model_32.reactions.get_by_id("Phloem_output_tx1").x+C3_model_32.reactions.get_by_id("Phloem_output_tx2").x)*12*60*60/1000,ls="--",alpha=0.3,color="red")
curve.plot(xvalues2[1],(C3_model_11.reactions.get_by_id("Phloem_output_tx1").x+C3_model_11.reactions.get_by_id("Phloem_output_tx2").x)*12*60*60/1000,marker="s",markersize=10,color="blue")
curve.axhline(y=(C3_model_11.reactions.get_by_id("Phloem_output_tx1").x+C3_model_11.reactions.get_by_id("Phloem_output_tx2").x)*12*60*60/1000,ls="--",alpha=0.3,color="blue")
curve.set_xlabel(r'$Vc/Vo$',fontsize=20)
curve.set_xticks(xvalues2)
curve.set_xticklabels(xvalues)
curve.set_ylabel("Phloem output ("+r'$mmol.m^{-2}.d^{-1}$'+")",fontsize=20)


#plt.xscale('log')
plt.legend(loc="best")
plt.show()

