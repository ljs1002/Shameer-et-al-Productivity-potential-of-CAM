#########################################################################
#This script can be used to estimate NGAM when 90% of phloem is sucrose	#
#growth condition based on the assumption that carbon conversion 	#
#efficiency (CCE) of night time metaoblism is 50%			#
#									#
#########################################################################


########################## FUNCTIONS ####################################

#This function predicts the carbon conversion efficiency of night time metabolism
def predictCCE(C3_model):
  import re
  for met in C3_model.metabolites:
    if not met.formula:
      met.formula=""
  Cin = 0
  Cout = 0
  for rxn in C3_model.reactions.query("dielTransfer"):
    if round(rxn.x,5)>0:
      for met in rxn.products:
        print met
        if met.formula.__contains__("C"):
          #print str(Cin)+"---"+met.id
          print str(rxn.x)+"\t"+str(rxn.metabolites.get(met))+"\t"+str(int(re.split(r"[C,H]",met.formula)[1]))
          Cin = Cin + (rxn.x * rxn.metabolites.get(met) * int(re.split(r"[C,H]",met.formula)[1]))
  
  for rxn in C3_model.reactions.query("dielTransfer"):
    if round(rxn.x,5)<0:
      for met in rxn.reactants:
        if met.formula.__contains__("C"):
          print str(Cout)+"---"+met.id
          print str(rxn.x)+"\t"+str(rxn.metabolites.get(met))+"\t"+str(int(re.split(r"[C,H]",met.formula)[1]))
          Cout = Cout + (rxn.x * rxn.metabolites.get(met) * int(re.split(r"[C,H]",met.formula)[1]))
  
  rxn = C3_model.reactions.get_by_id("Phloem_output_tx2")
  for met in rxn.reactants:
    if met.formula.__contains__("C"):
      print str(Cout)+"---"+met.id
      print str(rxn.x)+"\t"+str(-1 * rxn.metabolites.get(met))+"\t"+str(int(re.split(r"[C,H]",met.formula)[1]))
      Cout = Cout + (rxn.x * -1 * rxn.metabolites.get(met) * int(re.split(r"[C,H]",met.formula)[1]))
  
  Cout = Cout + (-1*C3_model.reactions.get_by_id("CO2_tx2").x)
  
  if(not round(Cin,5) == round(Cout,5)):
    print "Error, Cin = "+str(Cin)+" and Cout = "+str(Cout)
    return 0
  else:
    print "Cin = "+str(Cin)+"\tCO2 ="+str(C3_model.reactions.get_by_id("CO2_tx2").x)+"\t"+str(1 + ((C3_model.reactions.get_by_id("CO2_tx2").x)/Cin))
    return 1 + ((C3_model.reactions.get_by_id("CO2_tx2").x)/Cin)




################################ MAIN ######################################


#import libraries
from libsbml import readSBML
from cobra import io,flux_analysis
import re

#import model from SBML file
cobra_model = io.sbml.create_cobra_model_from_sbml_file("diel_model.xml")

#update metabolite charges
fin = open("home/sanu/Documents/FractionalCharges.csv","r")

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


#leaf light constraints
cobra_model.reactions.get_by_id("diel_biomass").objective_coefficient=1
cobra_model.reactions.get_by_id("Sucrose_tx1").lower_bound=0
cobra_model.reactions.get_by_id("Sucrose_tx1").upper_bound=0
cobra_model.reactions.get_by_id("GLC_tx1").lower_bound=0
cobra_model.reactions.get_by_id("GLC_tx1").upper_bound=0
cobra_model.reactions.get_by_id("CO2_tx1").lower_bound=0
cobra_model.reactions.get_by_id("NH4_tx1").lower_bound=0
cobra_model.reactions.get_by_id("NH4_tx1").upper_bound=0
#leaf dark constraints
cobra_model.reactions.get_by_id("Sucrose_tx2").lower_bound=0
cobra_model.reactions.get_by_id("Sucrose_tx2").upper_bound=0
cobra_model.reactions.get_by_id("GLC_tx2").lower_bound=0
cobra_model.reactions.get_by_id("GLC_tx2").upper_bound=0
cobra_model.reactions.get_by_id("Photon_tx2").lower_bound=0
cobra_model.reactions.get_by_id("Photon_tx2").upper_bound=0
cobra_model.reactions.get_by_id("NH4_tx2").lower_bound=0
cobra_model.reactions.get_by_id("NH4_tx2").upper_bound=0
cobra_model.reactions.get_by_id("CO2_tx2").upper_bound=0



#Set G6P transporter to 0 for C3
cobra_model.reactions.get_by_id("G6P_Pi_pc1").lower_bound=0
cobra_model.reactions.get_by_id("G6P_Pi_pc1").upper_bound=0
cobra_model.reactions.get_by_id("G6P_Pi_pc2").lower_bound=0
cobra_model.reactions.get_by_id("G6P_Pi_pc2").upper_bound=0

#Turn off PTOX
cobra_model.reactions.get_by_id("Plastoquinol_Oxidase_p1").lower_bound=0
cobra_model.reactions.get_by_id("Plastoquinol_Oxidase_p1").upper_bound=0

#nitrate uptake constrain
Nitrate_balance = Metabolite("Nitrate_bal_c", name = "Weights to balance nitrate uptake", compartment = "c1")
cobra_model.reactions.get_by_id("Nitrate_ec1").add_metabolites({Nitrate_balance:-2})
cobra_model.reactions.get_by_id("Nitrate_ec2").add_metabolites({Nitrate_balance:3})


#constraint VcVo
Rubisco_balance = Metabolite("rubisco_bal_p1", name = "Weights to balance RuBP carboxygenase oxygenase balance", compartment = "p1")
cobra_model.reactions.get_by_id("RXN_961_p1").add_metabolites({Rubisco_balance:3})
cobra_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1").add_metabolites({Rubisco_balance:-1})
#cobra_model2.reactions.get_by_id("RXN_961_p2").add_metabolites({Rubisco_balance:3})
#cobra_model2.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2").add_metabolites({Rubisco_balance:-1})

#constraint ATPase:NADPHoxidase = 3:1 and also ATPase (day)=ATPase(night)
Maintenance_constraint = Metabolite("ATPase_NADPHoxidase_constraint_c1",name =  "ATPase_NADPHoxidase_constraint_c1", compartment = "c1")
Maintenance_constraint2 = Metabolite("ATPase_NADPHoxidase_constraint_c2",name =  "ATPase_NADPHoxidase_constraint_c2", compartment = "c2")
Maintenance_constraint3 = Metabolite("Light_dark_maintainence_constraint",name =  "Light_dark_maintainence_constraint", compartment = "c1")
cobra_model.reactions.get_by_id("ATPase_tx1").add_metabolites({Maintenance_constraint:1,Maintenance_constraint3:1})
cobra_model.reactions.get_by_id("ATPase_tx2").add_metabolites({Maintenance_constraint2:1,Maintenance_constraint3:-1})
cobra_model.reactions.get_by_id("NADPHoxc_tx1").add_metabolites({Maintenance_constraint:-3})
cobra_model.reactions.get_by_id("NADPHoxc_tx2").add_metabolites({Maintenance_constraint2:-3})
cobra_model.reactions.get_by_id("NADPHoxm_tx1").add_metabolites({Maintenance_constraint:-3})
cobra_model.reactions.get_by_id("NADPHoxm_tx2").add_metabolites({Maintenance_constraint2:-3})
cobra_model.reactions.get_by_id("NADPHoxp_tx1").add_metabolites({Maintenance_constraint:-3})
cobra_model.reactions.get_by_id("NADPHoxp_tx2").add_metabolites({Maintenance_constraint2:-3})

##constrain sucrose, starch and fructan storage
cobra_model.reactions.get_by_id("SUCROSE_v_dielTransfer").upper_bound=0
cobra_model.reactions.get_by_id("SUCROSE_v_dielTransfer").upper_bound=0
cobra_model.reactions.get_by_id("FRUCTAN_v_dielTransfer").lower_bound=0
cobra_model.reactions.get_by_id("FRUCTAN_v_dielTransfer").upper_bound=0


#Plastid enolase was not detected in Arabidopsis mesophyll tissue
cobra_model.reactions.get_by_id("2PGADEHYDRAT_RXN_p1").lower_bound=0
cobra_model.reactions.get_by_id("2PGADEHYDRAT_RXN_p1").upper_bound=0


#Setting chloroplastic NADPH dehydrogenase to 0  ((Yamamoto et al., 2011)
cobra_model.reactions.get_by_id("NADPH_Dehydrogenase_p1").lower_bound=0
cobra_model.reactions.get_by_id("NADPH_Dehydrogenase_p1").upper_bound=0
cobra_model.reactions.get_by_id("NADPH_Dehydrogenase_p2").lower_bound=0
cobra_model.reactions.get_by_id("NADPH_Dehydrogenase_p2").upper_bound=0

#Set biomass to zero
cobra_model.reactions.get_by_id("Biomass_tx1").lower_bound=0
cobra_model.reactions.get_by_id("Biomass_tx1").upper_bound=0
cobra_model.reactions.get_by_id("Biomass_tx2").lower_bound=0
cobra_model.reactions.get_by_id("Biomass_tx2").upper_bound=0

#constrain sucrose to be 90% of the phloem
for i in range(1,3):
  cobra_model.reactions.get_by_id("Phloem_output_tx"+str(i)).add_metabolites({cobra_model.metabolites.get_by_id("sSUCROSE_b"+str(i)):-1.32,cobra_model.metabolites.get_by_id("PROTON_e"+str(i)):-1.32,cobra_model.metabolites.get_by_id("PROTON_c"+str(i)):1.32})



Orig = cobra_model.copy()	

#Find ATPase at which CCE ~ 50%
ATPaseValuesL=dict()
ATPaseValuesL2=dict()
DielBiomassValuesL=dict()
ATPaseValues=dict()
ATPaseValues2=dict()
DielBiomassValues=dict()
for a in range(100,1001,100):
  #constraining light to 500
  #a=500
  cobra_model2=Orig.copy()
  cobra_model2.reactions.get_by_id("Photon_tx1").lower_bound=a
  cobra_model2.reactions.get_by_id("Photon_tx1").upper_bound=a
  #get solution to check if model works
  solution=cobra_model2.optimize()
  print solution.f
  
  if solution.f is None:
    print "Solution infeasible at %s" % (a)
    continue
  
  solution=flux_analysis.parsimonious.optimize_minimal_flux(cobra_model2)
  cobra_model2.reactions.get_by_id("diel_biomass").objective_coefficient=0
  cobra_model2.reactions.get_by_id("ATPase_tx1").objective_coefficient=1
  cobra_model2.optimize()
  
  backup = cobra_model2.copy()
  
  j=0
  o=0
  
  for x in range(0,int(backup.solution.f)*10):
    j=j+1
    y=float(x)/10
    cobra_model=backup.copy()
    cobra_model.reactions.get_by_id("ATPase_tx1").lower_bound = y
    cobra_model.reactions.get_by_id("ATPase_tx1").upper_bound = y
    cobra_model.reactions.get_by_id("diel_biomass").objective_coefficient=1
    cobra_model.reactions.get_by_id("ATPase_tx1").objective_coefficient=0
    solution = flux_analysis.parsimonious.optimize_minimal_flux(cobra_model)
    
    CCE=predictCCE(cobra_model)
    print str(a)+"---\t"+str(y)+": max objective= "+str(solution.f)+"\t"+str(CCE)+"===="+str(cobra_model.reactions.get_by_id("CO2_tx2").x)
    
    if CCE < 0.50:
      ATPaseValuesL[a]=y
      ATPaseValuesL2[y]=solution.f
      print y
      DielBiomassValuesL[a]=solution.f
      break
  

#################################### LINE FITTING ################################
from scipy import stats
x=list()
y=list()
for k in ATPaseValuesL.keys():
  y.append(ATPaseValuesL.get(k))
  x.append(DielBiomassValuesL.get(k))

slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
if(p_value > 0.01):
  print "CHECK QUALITY"
else:
  print "ATPase =("+ str(slope)+" * diel_biomass) + "+str(intercept)



from scipy import stats
x=list()
y=list()
for k in ATPaseValuesL.keys():
  y.append(ATPaseValuesL.get(k))
  x.append(k)

slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
if(p_value > 0.01):
  print "CHECK QUALITY"
else:
  print "ATPase =("+ str(slope)+" * photon) + "+str(intercept)


