
rm /tmp/partial*

#rm ~/work/ppg14/trees/photon10_aa/caloana.root
#hadd -k -j 4 ~/work/ppg14/trees/photon10_aa/caloana.root photon10/condorout/OutDir*/caloana.root
#cd ~/work/ppg14/ppg12/FunWithxgboost/
#root -b -q "apply_BDT.C(\"config_split.yaml\", \"data\", \"/sphenix/u/bseidlitz/work/ppg14/trees/photon10_aa/caloana.root\") "
#cd - 

#rm ~/work/ppg14/trees/photon20_aa/caloana.root
#hadd -k -j 4 ~/work/ppg14/trees/photon20_aa/caloana.root test/condorout/OutDir*/caloana.root
#cd ~/work/ppg14/ppg12/FunWithxgboost/
#root -b -q "apply_BDT.C(\"config_split.yaml\", \"data\", \"/sphenix/u/bseidlitz/work/ppg14/trees/photon20_aa/caloana.root\") "
#cd - 

 rm ~/work/ppg14/trees/jet20_aa/caloana.root
 hadd -k -j 8 ~/work/ppg14/trees/jet20_aa/caloana.root jet20/condorout/OutDir*/caloana.root
 cd ~/work/ppg14/ppg12/FunWithxgboost/
 root -b -q "apply_BDT.C(\"config_split.yaml\", \"data\", \"/sphenix/u/bseidlitz/work/ppg14/trees/jet20_aa/caloana.root\") "
 cd - 

rm /tmp/partial*
