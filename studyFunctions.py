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



