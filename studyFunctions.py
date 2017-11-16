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
        #print met
        if met.formula.__contains__("C"):
          #print str(Cin)+"---"+met.id
          #print str(rxn.x)+"\t"+str(rxn.metabolites.get(met))+"\t"+str(int(re.split(r"[C,H]",met.formula)[1]))
          Cin = Cin + (rxn.x * rxn.metabolites.get(met) * int(re.split(r"[C,H]",met.formula)[1]))
  
  for rxn in C3_model.reactions.query("dielTransfer"):
    if round(rxn.x,5)<0:
      for met in rxn.reactants:
        if met.formula.__contains__("C"):
          #print str(Cout)+"---"+met.id
          #print str(rxn.x)+"\t"+str(rxn.metabolites.get(met))+"\t"+str(int(re.split(r"[C,H]",met.formula)[1]))
          Cout = Cout + (rxn.x * rxn.metabolites.get(met) * int(re.split(r"[C,H]",met.formula)[1]))
  
  rxn = C3_model.reactions.get_by_id("Phloem_output_tx2")
  for met in rxn.reactants:
    if met.formula.__contains__("C"):
      #print str(Cout)+"---"+met.id
      #print str(rxn.x)+"\t"+str(-1 * rxn.metabolites.get(met))+"\t"+str(int(re.split(r"[C,H]",met.formula)[1]))
      Cout = Cout + (rxn.x * -1 * rxn.metabolites.get(met) * int(re.split(r"[C,H]",met.formula)[1]))
  
  Cout = Cout + (-1*C3_model.reactions.get_by_id("CO2_tx2").x)
  
  if(not round(Cin,5) == round(Cout,5)):
    print "Error, Cin = "+str(Cin)+" and Cout = "+str(Cout)
    return 0
  else:
    print "Cin = "+str(Cin)+"\tCO2 ="+str(C3_model.reactions.get_by_id("CO2_tx2").x)+"\t"+str(1 + ((C3_model.reactions.get_by_id("CO2_tx2").x)/Cin))
    return 1 + ((C3_model.reactions.get_by_id("CO2_tx2").x)/Cin)



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
      


#convert moles of Phloem to moles of C
def convertOut2Cout(model,out):
  import re
  rxn = model.reactions.get_by_id("Phloem_output_tx1")
  Cout = 0
  for met in rxn.metabolites.keys():
    if met.formula.__contains__("C"):
      Cout = Cout + (-1* rxn.metabolites.get(met) * int(re.split(r"[C,H]",met.formula)[1]))
  
  Cout = out * Cout
  return Cout


#a function to calculate Sucrose to AA ratio
def predictSucAAratio(model):
  n=0
  s=0
  for met in model.reactions.get_by_id("Phloem_output_tx1").metabolites:
    if met.id == "X_Phloem_contribution_t1" or met.id=="GLC_c1" or met.id=="FRU_c1":
      continue
    elif met.id == "sSUCROSE_b1":
      s = model.reactions.get_by_id("Phloem_output_tx1").metabolites.get(met)
    else:
      n = n+model.reactions.get_by_id("Phloem_output_tx1").metabolites.get(met)
  return float(s)/n


#a function to calculate % of sucrose in phloem
def predictSucPercent(model):
  n=0
  s=0
  for met in model.reactions.get_by_id("Phloem_output_tx1").metabolites:
    if met.id == "X_Phloem_contribution_t1":
      continue
    elif met.id == "sSUCROSE_b1":
      s = model.reactions.get_by_id("Phloem_output_tx1").metabolites.get(met)
      n = n+model.reactions.get_by_id("Phloem_output_tx1").metabolites.get(met)
    else:
      n = n+model.reactions.get_by_id("Phloem_output_tx1").metabolites.get(met)
  return float(s)/n*100


#a function to calculate Sucrose to AA ratio
def predictSugAAratio(model):
  n=0
  s=0
  for met in model.reactions.get_by_id("Phloem_output_tx1").metabolites:
    if met.id == "X_Phloem_contribution_t1":
      continue
    elif met.id == "sSUCROSE_b1" or met.id=="GLC_c1" or met.id=="FRU_c1":
      s = s + model.reactions.get_by_id("Phloem_output_tx1").metabolites.get(met)
    else:
      n = n+model.reactions.get_by_id("Phloem_output_tx1").metabolites.get(met)
  return float(s)/n


#a function to calculate % of sucrose in phloem
def predictAAPercent(model):
  n=0
  a=0
  for met in model.reactions.get_by_id("Phloem_output_tx1").metabolites:
    if met.id == "X_Phloem_contribution_t1":
      continue
    elif met.id == "sSUCROSE_b1" or met.id == "GLC_c1" or met.id == "FRU_c1":
      n = n+model.reactions.get_by_id("Phloem_output_tx1").metabolites.get(met)
    else:
      a = a+model.reactions.get_by_id("Phloem_output_tx1").metabolites.get(met)
      n = n+model.reactions.get_by_id("Phloem_output_tx1").metabolites.get(met)
  return float(a)/n*100





