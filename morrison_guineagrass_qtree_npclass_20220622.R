
setwd("~/Documents/GradStudents/Morrison_Colin")

library(phytools)
library(picante)



feat = read.csv("~/Documents/GradStudents/Morrison_Colin/ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-a6ee92f8-download_qza_table_data/quantification_table/quantification_table-00000.csv")

dim(feat)
# [1] 8691  154

# meta = read.table("~/Documents/GradStudents/Greig_Keri/Metadata_HerbariumProject.txt", header = T, sep = "	", comment.char = "")
# meta = read.table("~/Documents/GradStudents/Greig_Keri/Metadata_HerbariumProject_xblanks.txt", header = T, sep = "	", comment.char = "")
tree = read.tree("~/Documents/GradStudents/Morrison_Colin/morr4_qemistree.nwk")
siri = read.table("~/Documents/GradStudents/Morrison_Colin/morr4_classified_feature_data.tsv", header = T, sep = "	", comment.char = "", , quote = "")
names(siri)[1] = "id"

dim(siri)
# [1] 5196   15

npclass = read.table(file = "~/Documents/GradStudents/Morrison_Colin/ProteoSAFe-NPCLASSIFIER-97cea901-view_results/NPCLASSIFIER-97cea901-view_results-main.tsv", header = T, sep = "\t", comment.char = "")

dim(npclass)
# [1] 3804    5

tree
# Phylogenetic tree with 5196 tips and 5195 internal nodes.

# Tip labels:
  # de15dd049c8d0fe63e0b5cb3369c95ac, 3388471770f94696e97be3d623ad72a1, 1cb86e64548f9d71fac44864d1346b4a, ba4d5ce67020581d977171af54a64b95, bb089437cf8e57b5dc91560472601f35, bfee5209e928688ecfa9b38f2eac91a8, ...

# Rooted; includes branch lengths.

nsamps = ncol(feat)-4
nsamps
# [1] 150


heat = as.data.frame(matrix(0,nrow = nrow(siri), ncol = (nsamps)))
row.names(heat) = siri$id
names(heat) = names(feat)[4:(ncol(feat)-1)]

dim(heat)
# [1] 2156  158

multips = grep(",", siri$table_number)
for(i in 1:nrow(siri)){
	labelid = as.character(siri$id[i])
	# featid = as.character(siri$X.featureID[i])
	featid = unlist(strsplit(as.character(siri$X.featureID[i]), split = ","))
	if(length(featid) == 1){
		table = siri$table_number[i]
		if(table == 1){fbmntable = feat}
		# if(table == 1){fbmntable = feat.kanupa44}
		# if(table == 2){fbmntable = feat.titiri42}
		# if(table == 3){fbmntable = feat.tintay25}
		for(j in 4:(ncol(fbmntable)-1)){
				samp = names(fbmntable)[j]
				heat[i,which(names(heat) == samp)] = fbmntable[which(fbmntable$row.ID == featid),j]
		}
	
	}
	if(length(featid) > 1){
		table = unlist(strsplit(as.character(siri$table_number[i]), split = ","))
		for(n in 1:length(featid)){
			# featidn = featid[n]
			if(table[n] == 1){fbmntable = feat}
			# if(table[n] == 1){fbmntable = feat.kanupa44}
			# if(table[n] == 2){fbmntable = feat.titiri42}
			# if(table[n] == 3){fbmntable = feat.tintay25}
			for(j in 4:(ncol(fbmntable)-1)){
				samp = names(fbmntable)[j]
				heat[i,which(names(heat) == samp)] = heat[i,which(names(heat) == samp)] + fbmntable[which(fbmntable$row.ID == featid[n]),j]
			}					
		}
	}	
}

dim(heat)
# [1] 5196  150

# length(which(meta$ATTRIBUTE_SampleType == "BLANK"))

# grep("LANK", names(heat))
# which(meta$ATTRIBUTE_SampleType == "BLANK")
# meta$filename[which(meta$ATTRIBUTE_SampleType == "BLANK")]

grep("lank", names(heat))
# [1] 2 4 5

# blanklist = c(meta$filename[which(meta$ATTRIBUTE_SampleType == "BLANK")])
# blankspots = rep(0,length(which(meta$ATTRIBUTE_SampleType == "BLANK")))
# for(i in 1:length(blanklist)){
	# blankspots[i] = grep(blanklist[i],names(heat))
# }
# blankspots
# heat.real = heat[which(rowSums(heat[,blankspots]) == 0),]


blanks = grep("lank", names(heat))
blanks
# [1] 2 4 5

heat.real = heat[-which(heat[,blanks] >0),]


dim(heat.real)
# [1] 3416  150

colSums(heat.real)
 # X8560_5_blank.mzXML.Peak.area   X8560_4_soil.mzXML.Peak.area   X8560_2_root.mzXML.Peak.area 
                             # 0                       12984875                     1223481515 
# X8560_3_litter.mzXML.Peak.area   X8560_1_leaf.mzXML.Peak.area 
                     # 531576447                     3119285276 

# colSums(heat.real)[blankspots]
# heat.real = heat.real[,-blankspots]

heat.real = heat.real[,-blanks]

dim(heat.real)
#[1] 3416  147

siri.heat = siri[which(siri$id %in% row.names(heat.real)),]

dim(siri.heat)
# [1] 3416   15



head(npclass)
                                                           smiles          class_results superclass_results
1 CC(=CCOC1=C(C2=C(C(=C1)O)C(=O)C(=C(O2)C3=CC4=C(C=C3)OCO4)O)OC)C              Flavonols         Flavonoids
2                       C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)CO)O)O)N  Purine nucleos(t)ides        Nucleosides
3                     CC(C(=O)NC(C)C(=O)O)NC(=O)C(CCC1(CNC1=O)O)N Dipeptides,Tripeptides     Small peptides
4                C=C1CSC(N=C1C(=O)O)C(C(=O)O)NC(=O)CCCCC(C(=O)O)N             Aminoacids     Small peptides
5              CCCCCCCCC=CC=CCCCCCCCCCCCCCCC=CC=CCCCC(=O)OC(CO)CO        Diacylglycerols      Glycerolipids
6              CCCCCCCCC=CC=CCCCCCCCCCCCCCCC=CC=CCCCC(=O)OCC(CO)O        Diacylglycerols      Glycerolipids
                  pathway_results isglycoside
1 Shikimates and Phenylpropanoids       False
2                   Carbohydrates       False
3        Amino acids and Peptides       False
4        Amino acids and Peptides       False
5                     Fatty acids       False
6                     Fatty acids       False



siri.itol = siri.heat

for(i in 1:nrow(siri.itol)){
	smiles_i = siri.itol$smiles[i]
	if(smiles_i %in% npclass$smiles){
		siri.itol$class_results[i] = npclass$class_results[which(npclass$smiles == smiles_i)]
		siri.itol$superclass_results[i] = npclass$superclass_results[which(npclass$smiles == smiles_i)]
		siri.itol$pathway_results[i] = npclass$pathway_results[which(npclass$smiles == smiles_i)]
		siri.itol$isglycoside[i] = npclass$isglycoside[which(npclass$smiles == smiles_i)]
	}
}






head(siri.itol)

levels(as.factor(siri.itol$pathway_results))

 # [1] ""                                                                  
 # [2] "Alkaloids"                                                         
 # [3] "Alkaloids,Amino acids and Peptides"                                
 # [4] "Alkaloids,Amino acids and Peptides,Shikimates and Phenylpropanoids"
 # [5] "Alkaloids,Shikimates and Phenylpropanoids"                         
 # [6] "Alkaloids,Terpenoids"                                              
 # [7] "Amino acids and Peptides"                                          
 # [8] "Amino acids and Peptides,Carbohydrates"                            
 # [9] "Amino acids and Peptides,Fatty acids"                              
# [10] "Amino acids and Peptides,Polyketides"                              
# [11] "Amino acids and Peptides,Shikimates and Phenylpropanoids"          
# [12] "Carbohydrates"                                                     
# [13] "Fatty acids"                                                       
# [14] "Fatty acids,Terpenoids"                                            
# [15] "Polyketides"                                                       
# [16] "Polyketides,Shikimates and Phenylpropanoids"                       
# [17] "Shikimates and Phenylpropanoids"                                   
# [18] "Shikimates and Phenylpropanoids,Terpenoids"                        
# [19] "Terpenoids"   



save(heat.real, siri.heat, heat, heat.real, tree,
# meta, 
feat, tree, siri, 
siri.itol,
file = "Morrison_GuineaGrass_metab_20220622.RData")


write.table(siri.itol, file = "Morrison_GuineaGrass_metab_20220622.tsv", sep = "\t", quote = F, row.names = F)










date = 20220210
outfile = paste("CSCSsppMorrisonGuineaGrass-", date, ".RData",sep="")

sampsByCompounds = as.data.frame(t(as.data.frame(heat.real)))
names(sampsByCompounds) = row.names(heat.real)

dim(sampsByCompounds)
# [1]  150 1504


# # # sppByCompounds = as.data.frame(matrix(0,ncol = ncol(sampsByCompounds), nrow = length(specieslist)))
# sppByCompounds = as.data.frame(matrix(0,ncol = nrow(heat), nrow = length(specieslist)))
# names(sppByCompounds) = names(sampsByCompounds)
# row.names(sppByCompounds) = specieslist

# for(i in 1:length(specieslist)){
	# sppByCompounds[i,] = colMeans(sampsByCompounds[grep(specieslist[i],row.names(sampsByCompounds)),])
# }
sppByCompounds = sampsByCompounds




## if calculating CSCS for species pairs rather than sample pairs, set to TRUE
# species = FALSE
# species = TRUE

# outfile = paste("CSCSsamples", date, ".RData",sep="")
# outfile = paste("CSCSspecies", date, ".RData",sep="")



sampsByCompounds = sppByCompounds





# bci.comm.mat = as.data.frame(matrix(1, nrow=2, ncol = ncol(sampsByCompoundsSiri)))
# names(bci.comm.mat) = names(sampsByCompoundsSiri)
comm.mat = as.data.frame(matrix(1, nrow=2, ncol = ncol(sampsByCompounds)))
names(comm.mat) = names(sampsByCompounds)
tree.bci = prune.sample(samp = comm.mat, phylo = tree)

maxdist = max(cophenetic(tree.bci))

# maxdist = pd.calc(cm = tree, tip.subset = c(tree$tip.label[1], tree$tip.label[length(tree$tip.label)]))[1]
# pairwise.comps = cophenetic(tree.bci)
pairwise.comps = 1-(cophenetic(tree.bci)/maxdist)

dim(pairwise.comps)
# [1] 1504 1504

max(pairwise.comps)
# [1] 1



# net.comps = c(levels(network$CLUSTERID1), levels(network$CLUSTERID2))
net.comps = nrow(pairwise.comps)
nspp = nrow(sampsByCompounds)
ncomps = ncol(sampsByCompounds)
pairwise.spp = as.data.frame(matrix(0,nrow = nspp, ncol = nspp))
names(pairwise.spp) = row.names(sampsByCompounds)
row.names(pairwise.spp) = row.names(sampsByCompounds)
sampsCompsStand = sampsByCompounds
for(i in 1:nrow(sampsByCompounds)){	
	sampsCompsStand[i,] = sampsByCompounds[i,]/sum(sampsByCompounds[i,])
}
diags = pairwise.spp
for (k in 1:nspp){
	sppX = as.character(row.names(sampsCompsStand)[k])
	cat("Comparing ", sppX, " to itself", "\n", sep = "")
	sppXonly = sampsCompsStand[k,which(sampsCompsStand[k,]>0)]
	ncomps = length(sppXonly)
	pairwise.comps.temp = pairwise.comps[names(sppXonly),names(sppXonly)]
	diags[k,k] = sum(((outer(as.numeric(sppXonly), as.numeric(sppXonly)))*pairwise.comps.temp),na.rm = T)	
}
save(sampsByCompounds, pairwise.comps, diags, file = outfile)
for (i in 1:nspp){
	spp1 = as.character(row.names(sampsCompsStand)[i])
	for (j in i:nspp){
		spp2 = as.character(row.names(sampsCompsStand)[j])
		cat("Comparing ", spp1, " to ", spp2, "\n", sep = "")
		#identify which compounds are in each species
		spp1comps = sampsCompsStand[spp1,]
		spp2comps = sampsCompsStand[spp2,]
		spp_pair = rbind(spp1comps,spp2comps)
		paircomps = spp_pair[,which(colSums(spp_pair)>0)]
		#make a pairwise.comps matrix for only those compounds found in either species in the species pair
		ncomps = ncol(paircomps)
		pairwise.comps.temp = pairwise.comps[names(paircomps),names(paircomps)]
		pairwise.spp[i,j] = pairwise.spp[j,i] = sum(((outer(as.numeric(paircomps[1,]), as.numeric(paircomps[2,])))*pairwise.comps.temp), na.rm = T)/max(diags[i,i], diags[j,j])
	}
}
cscs = pairwise.spp
save(sppByCompounds, sampsByCompounds, pairwise.comps, sampsCompsStand, diags, cscs, file = outfile)
	
cat("Completed cacluation of CSCS for all sample pairs","\n")






setwd("~/Documents/GradStudents/Morrison_Colin")

# save(heat.real, siri.heat, heat, 
# # meta, 
# feat, tree, siri, file = "Morrison_GuineaGrass_metab_20220210.RData")

load("Morrison_GuineaGrass_metab_20220210.RData")
ls()
# [1] "feat"      "heat"      "heat.real" "siri"      "siri.heat" "tree"    

dim(siri.heat)
# [1] 50 15

dim(heat.real)
# [1] 50  4

tree

# Phylogenetic tree with 118 tips and 117 internal nodes.

# Tip labels:
  # cc839d6329d1c8fb577114b5a8c713f8, 3a8499a4892a3ee234046fa58cde127a, 9f954da642b23dd7cae786b6ec3d331d, 21b9d719e5b035a3f3f3d7fc415f0c04, f6d0aa1b6594c5e85c112013dd94ba7c, 9fb7f6f88db6d39f86a6132a8285c9e6, ...

# Rooted; includes branch lengths.

npclass = read.table(file = "~/Documents/GradStudents/Morrison_Colin/ProteoSAFe-NPCLASSIFIER-97cea901-view_results/NPCLASSIFIER-97cea901-view_results-main.tsv", header = T, sep = "\t", comment.char = "")


dim(npclass)
# [1] 3804    5


head(npclass)
                                                           smiles          class_results superclass_results
1 CC(=CCOC1=C(C2=C(C(=C1)O)C(=O)C(=C(O2)C3=CC4=C(C=C3)OCO4)O)OC)C              Flavonols         Flavonoids
2                       C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)CO)O)O)N  Purine nucleos(t)ides        Nucleosides
3                     CC(C(=O)NC(C)C(=O)O)NC(=O)C(CCC1(CNC1=O)O)N Dipeptides,Tripeptides     Small peptides
4                C=C1CSC(N=C1C(=O)O)C(C(=O)O)NC(=O)CCCCC(C(=O)O)N             Aminoacids     Small peptides
5              CCCCCCCCC=CC=CCCCCCCCCCCCCCCC=CC=CCCCC(=O)OC(CO)CO        Diacylglycerols      Glycerolipids
6              CCCCCCCCC=CC=CCCCCCCCCCCCCCCC=CC=CCCCC(=O)OCC(CO)O        Diacylglycerols      Glycerolipids
                  pathway_results isglycoside
1 Shikimates and Phenylpropanoids       False
2                   Carbohydrates       False
3        Amino acids and Peptides       False
4        Amino acids and Peptides       False
5                     Fatty acids       False
6                     Fatty acids       False



siri.itol = siri.heat

for(i in 1:nrow(siri.itol)){
	smiles_i = siri.itol$smiles[i]
	if(smiles_i %in% npclass$smiles){
		siri.itol$class_results[i] = npclass$class_results[which(npclass$smiles == smiles_i)]
		siri.itol$superclass_results[i] = npclass$superclass_results[which(npclass$smiles == smiles_i)]
		siri.itol$pathway_results[i] = npclass$pathway_results[which(npclass$smiles == smiles_i)]
		siri.itol$isglycoside[i] = npclass$isglycoside[which(npclass$smiles == smiles_i)]
	}
}






head(siri.itol)

levels(as.factor(siri.itol$pathway_results))



































