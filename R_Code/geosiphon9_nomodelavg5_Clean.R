require(phytools)
#read MCC tree
tr<-read.nexus(file="...nuSSU_Parsed_rc_renamed_reduced_ginsi_trimmed_to_200m_mcctree_10pcntburnin_medhts.tre")

outs<-c("Endogone_pisformis_DQ322628_Endogone_Endogonaceae_Endogonales","Phycomyces_blakesleeanus_AY635837_Phycomyces_Phycomycetaceae_Mucorales","Rhizopus_oryzae_NG_062621_Rhizopus_Rhizopodaceae_Mucorales","Umbelopsis_ramanniana_DQ322627_Umbelopsis_Umbelopsidaceae_Umbelopsidales")

tr<-drop.tip(tr,outs)

#make states
df<-data.frame(matrix(nrow=length(tr$tip.label),ncol=7))
colnames(df)<-c("Taxon","State","Acc","VT","Genus","Family","Order")
df$Taxon<-tr$tip.label
rownames(df)<-df$Taxon

#straighten out taxonomy for Glomeromycotina
for(x in 1:nrow(df)){
	df$VT[x]<-strsplit(df$Taxon[x],"_")[[1]][1]
	df$Acc[x]<-strsplit(df$Taxon[x],"_")[[1]][2]
	df$Genus[x]<-strsplit(df$Taxon[x],"_")[[1]][3]
	df$Family[x]<-strsplit(df$Taxon[x],"_")[[1]][4]
	df$Order[x]<-strsplit(df$Taxon[x],"_")[[1]][5]
}
#add trait states
df$State<-"AM"
df$State[df$Genus %in% "Geosiphon"]<-"Cyanobacterial_Endosymbiont"

write.csv(df,file="...Glomeromycotina_States.csv")


require(phytools)
#read in trait data
tax<-read.csv(file="...Glomeromycotina_States.csv",stringsAsFactors=TRUE,row.names=1)
dat<-tax$State
names(dat)<-tax$Taxon
head(dat)


#read in trees
tr<-read.nexus(file="...nuSSU_Parsed_rc_renamed_reduced_ginsi_trimmed_to_200m_mcctree_10pcntburnin_medhts.tre")
tr<-read.nexus(file="...nuSSU_Parsed_rc_renamed_reduced_ginsi_trimmed_to_200m_mcctree_10pcntburnin_meanhts.tre")

outs<-c("Endogone_pisformis_DQ322628_Endogone_Endogonaceae_Endogonales","Phycomyces_blakesleeanus_AY635837_Phycomyces_Phycomycetaceae_Mucorales","Rhizopus_oryzae_NG_062621_Rhizopus_Rhizopodaceae_Mucorales","Umbelopsis_ramanniana_DQ322627_Umbelopsis_Umbelopsidaceae_Umbelopsidales")

tr<-drop.tip(tr,outs)
tr<-ladderize(tr,FALSE)

write.nexus(tr,file="...nuSSU_Parsed_rc_renamed_reduced_ginsi_trimmed_to_200m_mcctree_10pcntburnin_medhts_prunedladderized.tre")
tr<-read.nexus(file="...nuSSU_Parsed_rc_renamed_reduced_ginsi_trimmed_to_200m_mcctree_10pcntburnin_medhts_prunedladderized.tre")

require(phytools)
#read MCC tree
tr<-read.nexus(file="...nuSSU_Parsed_rc_renamed_reduced_ginsi_trimmed_to_200m_mcctree_10pcntburnin_medhts.tre")

#make states
df<-data.frame(matrix(nrow=length(tr$tip.label),ncol=8))
colnames(df)<-c("Taxon","State","Acc","VT","Genus","Family","Order","Ingroup")
df$Taxon<-tr$tip.label
rownames(df)<-df$Taxon

#straighten out taxonomy for Glomeromycotina
for(x in 5:nrow(df)){
	df$VT[x]<-strsplit(df$Taxon[x],"_")[[1]][1]
	df$Acc[x]<-strsplit(df$Taxon[x],"_")[[1]][2]
	df$Genus[x]<-strsplit(df$Taxon[x],"_")[[1]][3]
	df$Family[x]<-strsplit(df$Taxon[x],"_")[[1]][4]
	df$Order[x]<-strsplit(df$Taxon[x],"_")[[1]][5]
	df$Ingroup[x]<-"Yes"
}

#straighten out taxonomy for Outs
for(x in c(1,2,4)){
	df$VT[x]<-paste(strsplit(df$Taxon[x],"_")[[1]][1],strsplit(df$Taxon[x],"_")[[1]][2],sep="_")
	df$Acc[x]<-strsplit(df$Taxon[x],"_")[[1]][3]
	df$Genus[x]<-strsplit(df$Taxon[x],"_")[[1]][4]
	df$Family[x]<-strsplit(df$Taxon[x],"_")[[1]][5]
	df$Order[x]<-strsplit(df$Taxon[x],"_")[[1]][6]
	df$Ingroup[x]<-"No"
}
for(x in 3){
	df$VT[x]<-paste(strsplit(df$Taxon[x],"_")[[1]][1],strsplit(df$Taxon[x],"_")[[1]][2],sep="_")
	df$Acc[x]<-paste(strsplit(df$Taxon[x],"_")[[1]][3],strsplit(df$Taxon[x],"_")[[1]][4],sep="_")
	df$Genus[x]<-strsplit(df$Taxon[x],"_")[[1]][5]
	df$Family[x]<-strsplit(df$Taxon[x],"_")[[1]][6]
	df$Order[x]<-strsplit(df$Taxon[x],"_")[[1]][7]
	df$Ingroup[x]<-"No"
}


#add trait states
df$State<-"AM"
df$State[df$Genus %in% "Geosiphon"]<-"Cyanobacterial_Endosymbiont"
df$State[df$Genus %in% "Endogone"]<-"ECM"
df$State[df$Genus %in% c("Phycomyces","Rhizopus","Umbelopsis")]<-"Saprotroph"

write.csv(df,file="...Mucoromycota_States.csv")


#Stochastic Mapping
#Pi constrained to be non-ECM, Q="empirical"
require(phytools)
#read in distribution of trees
tr<-read.nexus(file="...nuSSU_Parsed_rc_renamed_reduced_ginsi_trimmed_to_200m_10pcntburnin.trees")

#select 100 from the 3601 trees here
sams<-sort(sample(1:3601,100,replace=FALSE))

tr.sams<-tr[sams]
#tr.sams
#  [1]    4   35   36   82  105  146  187  228  236  265  382  401  518  535  554
# [16]  559  605  630  637  663  714  784  785  805  826  842  854  896  913  931
# [31] 1084 1155 1181 1203 1267 1271 1272 1343 1361 1367 1375 1433 1505 1536 1545
# [46] 1553 1556 1693 1737 1774 1812 1817 1830 1853 1898 1901 1914 1921 1972 1973
# [61] 2017 2018 2019 2074 2276 2277 2326 2400 2413 2481 2517 2572 2577 2621 2629
# [76] 2651 2674 2690 2807 2874 2995 3004 3016 3030 3090 3100 3143 3146 3152 3220
# [91] 3225 3249 3285 3290 3292 3320 3416 3429 3485 3594

write.tree(tr.sams,file="...nuSSU_Parsed_rc_renamed_reduced_ginsi_trimmed_to_200m_10pcntburnin_100randomsamples.trees")

#read in trees
sams<-read.tree(file="...nuSSU_Parsed_rc_renamed_reduced_ginsi_trimmed_to_200m_10pcntburnin_100randomsamples.trees")
tr<-read.nexus(file="......nuSSU_Parsed_rc_renamed_reduced_ginsi_trimmed_to_200m_mcctree_10pcntburnin_medhts.tre")

#now add MCC tree to end as tree 101
boots.mod<-c(sams,tr)
class(boots.mod)<-"multiPhylo"
write.tree(boots.mod,file="...nuSSU_Parsed_rc_renamed_reduced_ginsi_trimmed_to_200m_10pcntburnin_100randomsamples_Tree101isMCC.trees")

#drop from test tree and all trees
require(phytools)
#read in trees
trs.clean<-read.tree(file="...nuSSU_Parsed_rc_renamed_reduced_ginsi_trimmed_to_200m_10pcntburnin_100randomsamples_Tree101isMCC.trees")

outs<-c("Endogone_pisformis_DQ322628_Endogone_Endogonaceae_Endogonales","Phycomyces_blakesleeanus_AY635837_Phycomyces_Phycomycetaceae_Mucorales","Rhizopus_oryzae_NG_062621_Rhizopus_Rhizopodaceae_Mucorales","Umbelopsis_ramanniana_DQ322627_Umbelopsis_Umbelopsidaceae_Umbelopsidales")

trs.cleaner<-NULL
for(x in 1:length(trs.clean)){
	trs.cleaner[[x]]<-drop.tip(trs.clean[[x]],outs)
}
class(trs.cleaner)<-"multiPhylo"
trs.cleaner[[1]]
write.tree(trs.cleaner,"...nuSSU_Parsed_rc_renamed_reduced_ginsi_trimmed_to_200m_10pcntburnin_100randomsamples_Tree101isMCC_noOuts.trees")

##############################################
#######PLOT MCC TREE - FIGURE S1
##############################################
require(phyloch)
require(phytools)
mcc<-read.beast(file="...nuSSU_Parsed_rc_renamed_reduced_ginsi_trimmed_to_200m_mcctree_10pcntburnin_medhts.tre")
png("...nuSSU_Parsed_rc_renamed_reduced_ginsi_trimmed_to_200m_mcctree_10pcntburnin_medhts_update.png",width=10,height=10,units="in",res=2400)

arch<-mcc$tip.label[grep("_Archaeosporales",mcc$tip.label)]
others<-mcc$tip.label[!mcc$tip.label %in% geo]
others<-others[!others %in% arch]
geo<-mcc$tip.label[grep("_Geosiphon_",mcc$tip.label)]
arch<-arch[!arch %in% geo]

#color tips
collabs<-mcc$tip.label
collabs[collabs %in% others]<-"darkgrey"
collabs[collabs %in% arch]<-"red"
collabs[collabs %in% geo]<-"deepskyblue"

plot(mcc,show.tip.label=TRUE,x.lim=c(-(650-max(nodeHeights(mcc))),max(nodeHeights(mcc))+55),edge.width=1,cex=0.15,tip.color=collabs)
HPDbars(mcc,col="darkolivegreen3",lwd=0.5)
node.support(mcc$posterior, cutoff = 0.95, mode = "dots",col="orange")
axisPhylo(cex.axis=0.6)

dev.off()
#plot(mcc,show.tip.label=TRUE,edge.width=3,add=TRUE)
#plot.phylo.upon(mcc,edge.width=3)


##############################################
######## Maximum Parsimony Analysis ##########
##############################################
require(castor)
require(ape)

#read in trait data
tax<-read.csv(file="...Mucoromycota_States.csv",stringsAsFactors=FALSE,row.names=1)

dat<-tax$State
names(dat)<-tax$Taxon
head(dat)

outs<-c("Endogone_pisformis_DQ322628_Endogone_Endogonaceae_Endogonales","Phycomyces_blakesleeanus_AY635837_Phycomyces_Phycomycetaceae_Mucorales","Rhizopus_oryzae_NG_062621_Rhizopus_Rhizopodaceae_Mucorales","Umbelopsis_ramanniana_DQ322627_Umbelopsis_Umbelopsidaceae_Umbelopsidales")

#remove outgroups
dat<-dat[!names(dat) %in% outs]

dat[dat=="AM"]<-1
dat[dat=="Cyanobacterial_Endosymbiont"]<-2

#read in trees
trs.cleaner<-read.tree(file="...nuSSU_Parsed_rc_renamed_reduced_ginsi_trimmed_to_200m_10pcntburnin_100randomsamples_Tree101isMCC_noOuts.trees")

costs1<-matrix(data=c(0,1,1,0),nrow=2,ncol=2)
costs2<-matrix(data=c(0,1,2,0),nrow=2,ncol=2)
costs3<-matrix(data=c(0,1,3,0),nrow=2,ncol=2)
costs4<-matrix(data=c(0,1,4,0),nrow=2,ncol=2)
costs5<-matrix(data=c(0,1,5,0),nrow=2,ncol=2)
costs6<-matrix(data=c(0,1,6,0),nrow=2,ncol=2)
costs7<-matrix(data=c(0,1,7,0),nrow=2,ncol=2)
costs8<-matrix(data=c(0,1,8,0),nrow=2,ncol=2)
costs9<-matrix(data=c(0,1,9,0),nrow=2,ncol=2)
costs10<-matrix(data=c(0,1,10,0),nrow=2,ncol=2)
costs11<-matrix(data=c(0,1,11,0),nrow=2,ncol=2)
costs12<-matrix(data=c(0,1,12,0),nrow=2,ncol=2)
costs13<-matrix(data=c(0,1,13,0),nrow=2,ncol=2)
costs14<-matrix(data=c(0,1,14,0),nrow=2,ncol=2)
costs15<-matrix(data=c(0,1,15,0),nrow=2,ncol=2)
costs16<-matrix(data=c(0,1,16,0),nrow=2,ncol=2)
costs17<-matrix(data=c(0,1,17,0),nrow=2,ncol=2)
costs18<-matrix(data=c(0,1,18,0),nrow=2,ncol=2)
costs19<-matrix(data=c(0,1,19,0),nrow=2,ncol=2)
costs20<-matrix(data=c(0,1,20,0),nrow=2,ncol=2)
costs50<-matrix(data=c(0,1,50,0),nrow=2,ncol=2)
costs100<-matrix(data=c(0,1,100,0),nrow=2,ncol=2)

#create a list of the matrices
cost.list<-list(costs1,costs2,costs3,costs4,costs5,costs6,costs7,costs8,costs9,costs10,costs11,costs12,costs13,costs14,costs15,costs16,costs17,costs18,costs19,costs20)


mpsample<-function(datadf,treeset,costmatrices){
	for(q in 1:length(costmatrices)){
		resdf<-NULL
		#make an empty dataframe for results
		cols<-c("Tree","All105","All10","NodeNo","Nodes1","Nodes2","NodesTie","TotalCost","RootState1","RootState2")
		resdf<-data.frame(matrix(nrow=101,ncol=length(cols)))
		colnames(resdf)<-cols
		for(x in 1:length(treeset)){
			resdf[x,"Tree"]<-x
			datn<-NULL
			dato<-NULL
			res<-NULL
			nodes<-NULL
			#order states the same as the tips
			datn<-datadf[order(match(datadf,treeset[[x]]$tip.label))]
			dato<-as.integer(datn)
			res<-asr_max_parsimony(tree= treeset[[x]],tip_states=dato,Nstates=2,transition_costs=costmatrices[[q]])
			#do all have node probabilities of 0, 0.5, 1?
			print(all(res$ancestral_likelihoods %in% c(0,1,0.5)))
			resdf[x,"All105"]<-all(res$ancestral_likelihoods %in% c(0,1,0.5))
			print(all(res$ancestral_likelihoods %in% c(0,1)))
			resdf[x,"All10"]<-all(res$ancestral_likelihoods %in% c(0,1))
			nodes<-NULL
			for(p in 1:nrow(res$ancestral_likelihoods)){
				if(p==1){
					if(all(res$ancestral_likelihoods[p,] %in% c(1,0))){
						nodes<-which.max(res$ancestral_likelihoods[p,])
					}else if(all(res$ancestral_likelihoods[p,] %in% c(0.5))){
						nodes<-"tie"
					}
				}
				if(p>1){
					if(all(res$ancestral_likelihoods[p,] %in% c(1,0))){
						nodes<-c(nodes,which.max(res$ancestral_likelihoods[p,]))
					}else if(all(res$ancestral_likelihoods[p,] %in% c(0.5))){
						nodes<-c(nodes,"tie")
					}	
				}
			}
			#nodes<-max.col(res$ancestral_likelihoods)
			print(table(nodes))
			resdf[x,"NodeNo"]<-length(nodes)	
			resdf[x,"Nodes1"]<-sum(nodes=="1")
			resdf[x,"Nodes2"]<-sum(nodes=="2")
			resdf[x,"NodesTie"]<-sum(nodes=="tie")	
			print(res$total_cost)
			resdf[x,"TotalCost"]<-res$total_cost	
			print(res$ancestral_likelihoods[1,])
			resdf[x,"RootState1"]<-res$ancestral_likelihoods[1,1]
			resdf[x,"RootState2"]<-res$ancestral_likelihoods[1,2]
		}
		write.csv(resdf,file=paste("~/geosiphon/Results_MP_Cost",q,".csv",sep=""))
	}
}

#run it.
mpsample(datadf=dat,treeset=trs.cleaner,costmatrices=cost.list)


#now, read files back in and extract info about root node and put in a dataframe
#make empty data frame first
rootst<-data.frame(matrix(ncol=20,nrow=101))
colnames(rootst)<-1:20

#probability it is cyanobacterial
for(x in 1:20){
	ind.df<-NULL
	ind.df<-read.csv(file=paste("...Results_MP_Cost",x,".csv",sep=""),stringsAsFactors=FALSE)
	rootst[1:101,x]<-ind.df["RootState2"]
}

#now, read files back in and extract info about total cost and put in a dataframe
#make empty data frame first
totalcost<-data.frame(matrix(ncol=20,nrow=101))
colnames(totalcost)<-1:20

#probability it is cyanobacterial
for(x in 1:20){
	ind.df<-NULL
	ind.df<-read.csv(file=paste("...Results_MP_Cost",x,".csv",sep=""),stringsAsFactors=FALSE)
	totalcost[1:101,x]<-ind.df["TotalCost"]
}

#Nodes1 - non-cyano
totalnodes1<-data.frame(matrix(ncol=20,nrow=101))
colnames(totalnodes1)<-1:20

#Number of nodes coded as non-cyano
for(x in 1:20){
	ind.df<-NULL
	ind.df<-read.csv(file=paste("...Results_MP_Cost",x,".csv",sep=""),stringsAsFactors=FALSE)
	totalnodes1[1:101,x]<-ind.df["Nodes1"]
}

#Nodes2 - cyano
totalnodes2<-data.frame(matrix(ncol=20,nrow=101))
colnames(totalnodes2)<-1:20

#Number of nodes coded as non-cyano
for(x in 1:20){
	ind.df<-NULL
	ind.df<-read.csv(file=paste("...Results_MP_Cost",x,".csv",sep=""),stringsAsFactors=FALSE)
	totalnodes2[1:101,x]<-ind.df["Nodes2"]
}

#NodesTie
totalnodestied<-data.frame(matrix(ncol=20,nrow=101))
colnames(totalnodestied)<-1:20

#Number of nodes tied
for(x in 1:20){
	ind.df<-NULL
	ind.df<-read.csv(file=paste("...Results_MP_Cost",x,".csv",sep=""),stringsAsFactors=FALSE)
	totalnodestied[1:101,x]<-ind.df["NodesTie"]
}


rootstt<-t(rootst)
colnames(rootstt)<-c(1:100,"MCC")
rownames(rootstt)<-paste("1:",1:20,sep="")

require(reshape)
df<-melt(rootstt)
df<-df[,c("X2","X1","value")]
colnames(df) <- c("Tree", "Cost", "value")

#reorder levels
df$Cost<-factor(df$Cost,levels=c(paste("1:",1:20,sep="")))
df$Tree<-factor(df$Tree,levels=c(1:101,"MCC"))

require(ggplot2)
rootimage<-ggplot(df, aes(x = Tree, y = Cost, fill = value)) + geom_tile(color = "black") + geom_text(aes(label = value), color = "white", size = 0.0001) + coord_fixed() + theme(text = element_text(size=6), axis.text.x = element_text(angle=90, hjust=1))
pdf(file="...Results_RootState.pdf")
rootimage
dev.off()


#####################
#####################Stochastic Mapping
#####################
#####################

conda activate r4-base
require("phytools")
require("MuMIn")

#read in trait data
tax<-read.csv(file="...Mucoromycota_States.csv",stringsAsFactors=TRUE,row.names=1)
dat<-tax$State
names(dat)<-tax$Taxon
head(dat)

outs<-c("Endogone_pisformis_DQ322628_Endogone_Endogonaceae_Endogonales","Phycomyces_blakesleeanus_AY635837_Phycomyces_Phycomycetaceae_Mucorales","Rhizopus_oryzae_NG_062621_Rhizopus_Rhizopodaceae_Mucorales","Umbelopsis_ramanniana_DQ322627_Umbelopsis_Umbelopsidaceae_Umbelopsidales")

#remove outgroups
dat<-dat[!names(dat) %in% outs]
dat<-droplevels(dat)

#read in trees
trs.cleaner<-read.tree(file="...nuSSU_Parsed_rc_renamed_reduced_ginsi_trimmed_to_200m_10pcntburnin_100randomsamples_Tree101isMCC_noOuts.trees")

#make empty data frame for results
cols<-c("Tree","ER_Lik","ER_AICc","ER_AICcw","ARD_Lik","ARD_AICc","ARD_AICcw","Best","Changes","AM_Gains","Cyano_Gains")
df<-data.frame(matrix(nrow=length(trs.cleaner),ncol=length(cols)))
colnames(df)<-cols
head(df)

#set numbers for aicc
n<-length(dat)
k.er<-1
k.ard<-2

for(z in 1:length(trs.cleaner)){
	cat(paste("Beginning tree ", z, "\n",sep=""))
	tr<-NULL
	tr<-trs.cleaner[[z]]

	#determine whether ER or ARD is best fit using fitMk and calculating AIC and AICcweights
	mk.er<-NULL
	mk.ard<-NULL
	mk.er<-fitMk(tree=tr,x=dat,model="ER",pi="fitzjohn")
	df[z,c("ER_Lik")]<-mk.er$logLik
	mk.ard<-fitMk(tree=tr,x=dat,model="ARD",pi="fitzjohn")
	df[z,c("ARD_Lik")]<-mk.ard$logLik

	#from geiger:::.aic
	aicc.vals<-NULL
	aiccw<-NULL
	df[z,"ER_AICc"]<-2*k.er*(n/(n-k.er-1))-(2*df[z,"ER_Lik"])
	df[z,"ARD_AICc"]<-2*k.ard*(n/(n-k.ard-1))-(2*df[z,"ARD_Lik"])
	aicc.vals<-c(df[z,"ER_AICc"],df[z,"ARD_AICc"])
	names(aicc.vals)<-c("ER","ARD")
	
	#get weights using MuMIn
	aiccw<-Weights(aicc.vals)
	df[z,"ER_AICcw"]<-aiccw[[1]]
	df[z,"ARD_AICcw"]<-aiccw[[2]]
	df[z,"Best"]<-names(aiccw)[aiccw %in% max(aiccw)]

	#make simmap with best-fit model
	tr.fit<-NULL
	qq<-NULL
	tr.fit<-make.simmap(tree=tr,x=dat,model=df[z,"Best"],pi="fitzjohn",Q="empirical",nsim=1000)
	qq<-colMeans(countSimmap(tr.fit)$Tr)
	df[z,"Changes"]<-qq["N"][[1]]
	df[z,"Cyano_Gains"]<-qq["AM,Cyanobacterial_Endosymbiont"][[1]]
	df[z,"AM_Gains"]<-qq["Cyanobacterial_Endosymbiont,AM"][[1]]
	save(tr.fit,dat,tr,file=paste("...SIMMAP_1K_RootFitzJohn_empirical_Glomeromycotina2state_ARDvsER_Tree_",z,".RData",sep=""))
	cat(paste(df[z,"Cyano_Gains"], " Cyano gains; ", df[z,"AM_Gains"], " AM gains\n",sep=""))
	cat(paste(z, "_finished\n",sep=""))
}

summary(df$Cyano_Gains)
summary(df$AM_Gains)
sum(df$AM_Gains>0.999)
table(df$Best)
write.csv(df,file="...SIMMAP_1K_RootFitzJohn_empirical_Glomeromycotina2state_ARDvsER_RESULTS.csv")
#df<-read.csv(file="...SIMMAP_1K_RootFitzJohn_empirical_Glomeromycotina2state_MODELAVERAGED_RESULTS.csv",stringsAsFactors=FALSE)

#make density maps for all and save
require(phytools)
for(z in 1:101){
	tr.fit<-NULL
	dat<-NULL
	tr<-NULL
	obj<-NULL
	load(paste("...SIMMAP_1K_RootFitzJohn_empirical_Glomeromycotina2state_ARDvsER_Tree_",z,".RData",sep=""))
	obj<-densityMap(tr.fit,states=levels(dat)[1:2],plot=FALSE)
	save(tr.fit,dat,tr,obj,file=paste("...SIMMAP_1K_RootFitzJohn_empirical_Glomeromycotina2state_ARDvsER_wDensityMap_Tree_",z,".RData",sep=""))
}



###################
####GENERATE CTT's
###################
require(phytools)

#load ctt.mod2 function below

#make an empty list for results of analyses done on sample of trees
ctt.to.plot.mod.muc<-vector(mode="list",length=101)

#Get analyses done on sample of trees and ML tree...use tr.fit, rather than obj (which is the density map)
for(x in 1:101){
	tr.fit<-NULL
	load(paste("...SIMMAP_1K_RootFitzJohn_empirical_Glomeromycotina2state_ARDvsER_wDensityMap_Tree_",x,".RData",sep=""))
	ctt.to.plot.mod.muc[[x]]<-ctt.mod2(tr.fit,timebin=10,direction="AM->Cyanobacterial_Endosymbiont")
	tr.fit<-NULL
}
#ctt.to.plot.mod.muc[[2]][2]
#ctt.to.plot.mod.muc[[101]][2]

#make two empty lists for results of analyses done on sample of trees; one for changes, and one for change rates
ch.muc<-vector(mode="list",length=101)
ch.r.muc<-vector(mode="list",length=101)

for(x in 1:101){
	#changes
	ch.muc[[x]]<-as.data.frame(cbind(max(ctt.to.plot.mod.muc[[x]]$segments)-as.vector(t(ctt.to.plot.mod.muc[[x]]$segments)),as.vector(rbind(ctt.to.plot.mod.muc[[x]]$nchanges,ctt.to.plot.mod.muc[[x]]$nchanges))))
	colnames(ch.muc[[x]])<-c("Time","Change")
	ch.muc[[x]]$Kind<-x
	
	#change rate
	ch.r.muc[[x]]<-as.data.frame(cbind(max(ctt.to.plot.mod.muc[[x]]$segments)-as.vector(t(ctt.to.plot.mod.muc[[x]]$segments)),as.vector(rbind(ctt.to.plot.mod.muc[[x]]$nchanges/ctt.to.plot.mod.muc[[x]]$edge.length,ctt.to.plot.mod.muc[[x]]$nchanges/ctt.to.plot.mod.muc[[x]]$edge.length))))
	colnames(ch.r.muc[[x]])<-c("Time","Change")
	ch.r.muc[[x]]$Kind<-x
}

require(dplyr)
ch.muc<-bind_rows(ch.muc)
ch.r.muc<-bind_rows(ch.r.muc)

ch.comb.muc<-ch.muc
colnames(ch.comb.muc)<-c("Time","ChangeNo","Kind")
ch.comb.muc<-cbind(ch.comb.muc,ch.r.muc$Change)

colnames(ch.comb.muc)<-c("Time","ChangeNo","Kind","ChangeRate")

ch.comb.muc$Kind[ch.comb.muc$Kind %in% "101"]<-"MCC"
mycols.muc<-c(rep("grey",100),"black")

write.csv(ch.comb.muc,file="...ch.comb.muc.csv")

###################
####Jointly PLOT CTT's TOGETHER ON ONE PAGE NICELY
###################

require(ggplot2)
require(gridExtra)
require(phytools)

ch.comb.muc<-read.csv(file="...ch.comb.muc.csv",stringsAsFactors=FALSE,row.names=1)
mycols.muc<-c(rep("grey",100),"black")

#read timescale in
timescale_ics2020<-read.csv(file="...timescale_ics2020.csv")

#reduce to just periods of interest
timescale_ics2020<-timescale_ics2020[timescale_ics2020$Type %in% "Period",]
timey<-timescale_ics2020[timescale_ics2020$Start<550,]

#add Ediacaran, because some extend into it
timey[13,"End"]<-541.000
timey[13,"Start"]<-635.000
timey[13,"Midpoint"]<-588.000
timey[13,"Col_R"]<-254
timey[13,"Col_G"]<-217
timey[13,"Col_B"]<-106

timey$Name<-c(NA,"Ng","Pg","K","J","T","P","C","D","S","O","Ca","Ed")

#make rgb color
for(x in 1:nrow(timey)){
	timey$RGB[x]<-rgb(timey$Col_R[x]/255,timey$Col_G[x]/255,timey$Col_B[x]/255)
}


timeplot.muc<-ggplot(ch.comb.muc,aes(x=Time,ymin=-0.025))+scale_x_reverse(expand=c(0,0),limits=c(560,0), breaks=c(560,0))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x.top=element_blank(), axis.title.x.bottom=element_blank(), axis.line.x=element_blank(), panel.background=element_blank(), axis.text.x=element_text(color="black"), axis.text.y=element_text(color="black"),axis.ticks.x = element_blank(), axis.ticks.y = element_line(color="black")) + annotate(geom='segment',y=-Inf,yend=Inf,x=Inf,xend=Inf) + annotate(geom='segment',y=-Inf,yend=Inf,x=-Inf,xend=-Inf)+annotate("rect", xmin=timey$End[1],xmax=timey$Start[1],ymin=-Inf,ymax=-0.01,fill=timey$RGB[1],color="black",alpha=1)+annotate("rect",xmin=timey$End[2],xmax=timey$Start[2],ymin=-Inf,ymax=-0.01,fill=timey$RGB[2],color="black",alpha=1)+annotate("rect",xmin=timey$End[3],xmax=timey$Start[3],ymin=-Inf,ymax=-0.01,fill=timey$RGB[3],color="black",alpha=1)+annotate("rect",xmin=timey$End[4],xmax=timey$Start[4],ymin=-Inf,ymax=-0.01,fill=timey$RGB[4],color="black",alpha=1)+annotate("rect",xmin=timey$End[5],xmax=timey$Start[5],ymin=-Inf,ymax=-0.01,fill=timey$RGB[5],color="black",alpha=1)+annotate("rect",xmin=timey$End[6],xmax=timey$Start[6],ymin=-Inf,ymax=-0.01,fill=timey$RGB[6],color="black",alpha=1)+annotate("rect",xmin=timey$End[7],xmax=timey$Start[7],ymin=-Inf,ymax=-0.01,fill=timey$RGB[7],color="black",alpha=1)+annotate("rect",xmin=timey$End[8],xmax=timey$Start[8],ymin=-Inf,ymax=-0.01,fill=timey$RGB[8],color="black",alpha=1)+annotate("rect",xmin=timey$End[9],xmax=timey$Start[9],ymin=-Inf,ymax=-0.01,fill=timey$RGB[9],color="black",alpha=1)+annotate("rect",xmin=timey$End[10],xmax=timey$Start[10],ymin=-Inf,ymax=-0.01,fill=timey$RGB[10],color="black",alpha=1)+annotate("rect",xmin=timey$End[11],xmax=timey$Start[11],ymin=-Inf,ymax=-0.01,fill=timey$RGB[11],color="black",alpha=1)+annotate("rect",xmin=timey$End[12],xmax=timey$Start[12],ymin=-Inf,ymax=-0.01,fill=timey$RGB[12],color="black",alpha=1)+annotate("rect",xmin=timey$End[13],xmax=560,ymin=-Inf,ymax=-0.01,fill=timey$RGB[13],color="black",alpha=1)+annotate("text",x=timey$Midpoint[12],y=-0.025,label=timey$Name[12],size=3,colour="black")+annotate("text",x=timey$Midpoint[11],y=-0.025,label=timey$Name[11],size=3,colour="black")+annotate("text",x=timey$Midpoint[10],y=-0.025,label=timey$Name[10],size=3,colour="black")+annotate("text",x=timey$Midpoint[9],y=-0.025,label=timey$Name[9],size=3,colour="black")+annotate("text",x=timey$Midpoint[8],y=-0.025,label=timey$Name[8],size=3,colour="black")+annotate("text",x=timey$Midpoint[7],y=-0.025,label=timey$Name[7],size=3,colour="black")+annotate("text",x=timey$Midpoint[6],y=-0.025,label=timey$Name[6],size=3,colour="black")+annotate("text",x=timey$Midpoint[5],y=-0.025,label=timey$Name[5],size=3,colour="black")+annotate("text",x=timey$Midpoint[4],y=-0.025,label=timey$Name[4],size=3,colour="black")+annotate("text",x=timey$Midpoint[3],y=-0.025,label=timey$Name[3],size=3,colour="black")+annotate("text",x=timey$Midpoint[2],y=-0.025,label=timey$Name[2],size=3,colour="black")+annotate("text",x=1.3,y=-0.025,label=timey$Name[1],size=3,colour="black") 

timeplot.muc<-timeplot.muc+ annotate(geom='segment',y=0.15,yend=0.15,x=475,xend=450,color=mycols.muc[1],linetype=1) + annotate(geom='segment',y=0.2,yend=0.2,x=475,xend=450,color=mycols.muc[101],linetype=1) + annotate(geom="text", x=405, y=0.15, label="Posterior",color="grey",hjust=0,size=2.5)+ annotate(geom="text", x=405, y=0.2, label="MCC",color="black",hjust=0,size=2.5)+ annotate(geom="text", x=410, y=0.25, label="Evolution of cyanobacterial associations in Glomeromycotina",color="black",hjust=0,size=3.5,fontface=2)

timeplot.muc.comb<-timeplot.muc+geom_step(data=subset(ch.comb.muc[ch.comb.muc$Kind=="1",]),linetype=1,alpha=1,color=mycols.muc[1],aes(y=ChangeNo))
timeplot.muc.combb<-timeplot.muc.comb+geom_step(data=subset(ch.comb.muc[ch.comb.muc$Kind=="2",]),linetype=1,alpha=1,color=mycols.muc[2],aes(y=ChangeNo))
for(x in 3:100){
	timeplot.muc.combb<-timeplot.muc.combb+geom_step(data=subset(ch.comb.muc[ch.comb.muc$Kind==x,]),linetype=1,alpha=1,color=mycols.muc[x],aes(y=ChangeNo))
}
timeplot.muc.combb<-timeplot.muc.combb+geom_step(data=subset(ch.comb.muc[ch.comb.muc$Kind=="MCC",]),linetype=1,alpha=1,color=mycols.muc[101],aes(y=ChangeNo))


timeplot.mucd<-timeplot.muc.combb+scale_y_continuous(name="Mean number of gains /\n10 my time bin")+theme(axis.title.y = element_text(size=8))
timeplot.muce<-timeplot.mucd
pdf(file="...mucoromycota_test_ARDvsER.pdf",width=7,height=2)
timeplot.muce
dev.off()






timeplot.muc.combbb<-timeplot.muc.combb+geom_step(data=subset(ch.comb.muc[ch.comb.muc$Kind=="1",]),linetype=3,alpha=1,color=mycols.muc[1],aes(y=ChangeRate*666))
timeplot.mucc<-timeplot.muc.combbb+geom_step(data=subset(ch.comb.muc[ch.comb.muc$Kind=="2",]),linetype=3,alpha=1,color=mycols.muc[2],aes(y=ChangeRate*666))
for(x in 3:100){
	timeplot.mucc<-timeplot.mucc+geom_step(data=subset(ch.comb.muc[ch.comb.muc$Kind==x,]),linetype=3,alpha=1,color=mycols.muc[x],aes(y=ChangeRate*666))
}
timeplot.mucc<-timeplot.mucc+geom_step(data=subset(ch.comb.muc[ch.comb.muc$Kind=="MCC",]),linetype=3,alpha=1,color=mycols.muc[101],aes(y=ChangeRate*666))

timeplot.mucd<-timeplot.mucc+scale_y_continuous(name="Mean number of gains /\n10 my time bin", sec.axis = sec_axis(~./666, name = "Mean number of gains /\n total branch length / 10 my time bin"))+theme(axis.title.y = element_text(size=8))
timeplot.muce<-timeplot.mucd
pdf(file="...mucoromycota_ARDvsER.pdf",width=7,height=2)
timeplot.muce
dev.off()


#############################################Did this on PC
require(colorspace)
require(plotrix)
require(phytools)

draw.half.circle<-function (x, y, radius, nv = 100, border = NULL, col = NA, lty = 1, 
    density = NULL, angle = 0, lwd = 1) 
{
    xylim <- par("usr")
    plotdim <- par("pin")
    ymult <- getYmult()
    angle.inc <- 2 * pi/nv
    angles <- seq(0, 2 * pi - angle.inc, by = angle.inc)
    if (length(col) < length(radius)) 
        col <- rep(col, length.out = length(radius))
    for (circle in 1:length(radius)) {
       # xv <- cos(angles) * radius[circle] + x
       # yv <- sin(angles) * radius[circle] * ymult + y
        xv <- cos(angles) * radius[circle] + x
        yv <- sin(angles) * radius[circle] * ymult + y
        polygon(xv[1:51], yv[1:51], border = border, col = col[circle], lty = lty, 
            density = density, angle = angle, lwd = lwd)
    }
    invisible(list(x = xv[1:51], y = yv[1:51]))
}

load('...SIMMAP_1K_RootFitzJohn_empirical_Glomeromycotina2state_ARDvsER_wDensityMap_Tree_101.RData')
# tr.fit, and obj (which is the density map)

#change color gradient
obj<-setMap(obj,c("darkgrey","deepskyblue"))
plot(obj$tree,obj$cols,0.5,lwd=1,mar=par()$mar,type="fan",part=0.5,ftype="off")

#modified from http://blog.phytools.org/2018/01/adding-geological-legend-to-fan-style.html
#and http://blog.phytools.org/2015/08/more-updates-associated-with-mapped.html
obj2<-geo.legend()

#read in ICS 2020
timescale_ics2020<-read.csv(file="...timescale_ics2020.csv")

for(x in 1:nrow(obj2$leg)){
	obj2$leg[x,1]<-timescale_ics2020$Start[x]
	obj2$leg[x,2]<-timescale_ics2020$End[x]
	obj2$colors[x]<-rgb(timescale_ics2020$Col_R[x], timescale_ics2020$Col_G[x], timescale_ics2020$Col_B[x],max=255)
}

r2<-max(obj2$leg[,1])-obj2$leg[,2]

plot(obj$tree,obj$cols,0.5,lwd=1,mar=par()$mar,type="fan",part=0.5,ftype="off")
for(i in 1:nrow(obj2$leg)){
    color<-paste(strsplit(obj2$colors[i],"")[[1]][1:7],collapse="")
    print(color)
    draw.half.circle(0,0,radius=r2[i],col=lighten(color,0.5,space="combined"),border="transparent")
}
par(fg="transparent")
plot(obj$tree,obj$cols,0.5,lwd=1.25,mar=par()$mar,type="fan",part=0.5,ftype="off",add=TRUE)
par(fg="black")

#mid point - skip quaternary
mids<-r2[[1]]-timescale_ics2020$Midpoint[2:11]
for(x in 2:11){
	text(mids[x-1],-10, timescale_ics2020$Abbrev[x],col="black",cex=0.5)
}

legend(350,400,c("AMF-only",expression(paste(italic("Geosiphon"),"-like",sep=""))),pch=15,col=c("darkgrey","deepskyblue"),pt.cex=1.5,cex=0.5,bty="n",text.col="black")

#read taxonomy
tax<-read.csv(file="...Mucoromycota_States.csv",stringsAsFactors=FALSE)

#drop outgroups
tax<-tax[tax$Ingroup=="Yes",]

table(tax$Order)
ords<-unique(tax$Order)


tax$Class<-NA
tax$Class[tax$Order=="Paraglomerales"]<-"Paraglomeromycetes"
tax$Class[tax$Order=="Archaeosporales"]<-"Archaeospomeromycetes"
tax$Class[tax$Order=="Glomerales"]<-"Glomeromycetes"
tax$Class[tax$Order=="Diversisporales"]<-"Glomeromycetes"


arch<-tax$Taxon[tax$Order %in% "Archaeosporales"]
arc.cladelabels(text="Archaeosporales",node=findMRCA(obj$tree,arch),mark.node=FALSE,col="black",cex=0.5)

#########################
#########################BiSSE in castor
#########################
#########################


conda activate r4-base
require(castor)
require(ape)
require(MuMIn)

#read in trait data
tax<-read.csv(file="...Mucoromycota_States.csv",stringsAsFactors=TRUE,row.names=1)
dat<-tax$State
names(dat)<-tax$Taxon
head(dat)

#remove outgroups
outs<-c("Endogone_pisformis_DQ322628_Endogone_Endogonaceae_Endogonales","Phycomyces_blakesleeanus_AY635837_Phycomyces_Phycomycetaceae_Mucorales","Rhizopus_oryzae_NG_062621_Rhizopus_Rhizopodaceae_Mucorales","Umbelopsis_ramanniana_DQ322627_Umbelopsis_Umbelopsidaceae_Umbelopsidales")
dat<-dat[!names(dat) %in% outs]
dat<-droplevels(dat)
states<-rep(1,length(dat))
names(states)<-names(dat)
states[names(states) %in% "VTX00241_Y15904_Geosiphon_Geosiphonaceae_Archaeosporales"]<-2
table(states)

#read in trees
trs.cleaner<-read.tree(file="...nuSSU_Parsed_rc_renamed_reduced_ginsi_trimmed_to_200m_10pcntburnin_100randomsamples_Tree101isMCC_noOuts.trees")

#make data frame and fit under Yule model
gens<-c("Lik","AIC","AM_to_Geo","Geo_to_AM","Birth_AM","Birth_Geo","Death_AM","Death_Geo","Root_is_AM","Root_is_Geo")
cols<-c("Tree",paste("MFRoot",gens,sep="_"))
df<-as.data.frame(matrix(nrow=101,ncol=length(cols)))
colnames(df)<-cols

for(x in 1:101){
	#blank vectors
	states.ordered<-NULL
	root.mf<-NULL
	root.geo<-NULL
	df$Tree[x]<-x
	
	#re-order to match order in phylogeny under investigation
	states.ordered<-states[order(match(names(states),trs.cleaner[[x]]$tip.label))]
	
	#fit w root conditioned using madfitz
	root.mf<-fit_musse(tree=trs.cleaner[[x]],Nstates=2,NPstates=0,tip_pstates=states.ordered,transition_rate_model="ARD",birth_rate_model="ARD",death_rate_model="ER",death_rates=c(0,0),root_prior="likelihoods",root_conditioning="madfitz",include_ancestral_likelihoods=TRUE)
	df[x,"MFRoot_Lik"]<-root.mf$loglikelihood
	df[x,"MFRoot_AIC"]<-root.mf$AIC
	df[x,"MFRoot_AM_to_Geo"]<-root.mf$parameters$transition_matrix[1,2]
	df[x,"MFRoot_Geo_to_AM"]<-root.mf$parameters$transition_matrix[2,1]
	df[x,"MFRoot_Birth_AM"]<-root.mf$parameters$birth_rates[1]
	df[x,"MFRoot_Birth_Geo"]<-root.mf$parameters$birth_rates[2]
	df[x,"MFRoot_Death_AM"]<-root.mf$parameters$death_rates[1]
	df[x,"MFRoot_Death_Geo"]<-root.mf$parameters$death_rates[2]
	df[x,"MFRoot_Root_is_AM"]<-root.mf$ancestral_likelihoods[1,1]
	df[x,"MFRoot_Root_is_Geo"]<-root.mf$ancestral_likelihoods[1,2]
}

summary(df$MFRoot_Root_is_AM)
head(df)
summary(df[1:100,])
df[101,]
write.csv(df,"...geosiphon_bisse_results_birthonly.csv")

#make new data frame and fit under birth-death model
df2<-as.data.frame(matrix(nrow=101,ncol=length(cols)))
colnames(df2)<-cols
for(x in 1:101){
	#blank vectors
	states.ordered<-NULL
	root.mf<-NULL
	root.geo<-NULL
	df2$Tree[x]<-x
	
	#re-order to match order in phylogeny under investigation
	states.ordered<-states[order(match(names(states),trs.cleaner[[x]]$tip.label))]
	
	#fit w root conditioned using madfitz
	root.mf<-fit_musse(tree=trs.cleaner[[x]],Nstates=2,NPstates=0,tip_pstates=states.ordered,transition_rate_model="ARD",birth_rate_model="ARD",death_rate_model="ARD",root_prior="likelihoods",root_conditioning="madfitz",include_ancestral_likelihoods=TRUE)
	df2[x,"MFRoot_Lik"]<-root.mf$loglikelihood
	df2[x,"MFRoot_AIC"]<-root.mf$AIC
	df2[x,"MFRoot_AM_to_Geo"]<-root.mf$parameters$transition_matrix[1,2]
	df2[x,"MFRoot_Geo_to_AM"]<-root.mf$parameters$transition_matrix[2,1]
	df2[x,"MFRoot_Birth_AM"]<-root.mf$parameters$birth_rates[1]
	df2[x,"MFRoot_Birth_Geo"]<-root.mf$parameters$birth_rates[2]
	df2[x,"MFRoot_Death_AM"]<-root.mf$parameters$death_rates[1]
	df2[x,"MFRoot_Death_Geo"]<-root.mf$parameters$death_rates[2]
	df2[x,"MFRoot_Root_is_AM"]<-root.mf$ancestral_likelihoods[1,1]
	df2[x,"MFRoot_Root_is_Geo"]<-root.mf$ancestral_likelihoods[1,2]
}

summary(df2$MFRoot_Root_is_AM)
head(df2)
summary(df2[1:100,])
df2[101,]
write.csv(df2,"...geosiphon_bisse_results_birthdeath.csv")


#########################
####EXTRA FUNCTIONS
#########################
#modified in Boyce et al. 2023 Geobiology from the ctt function in phytools.  

#modified so user can specify timebin duration - truncates oldest one near root.
#standardizing timebins allows users to plot trees of different heights together with standardized time windows.
#also modified so use can specify which changes specifically to look at.
ctt.mod2<-function (trees, timebin = 5,direction=NULL){
	#modified from phytools 0.7.70
    if (!(inherits(trees, "multiSimmap"))) 
        stop("trees should be an object of class \"multiSimmap\".")
    tree <- as.phylo(trees[[1]])
    changes <- sapply(trees, phytools:::getChanges)
    #need this to be a list...if only one change/tree, it will be numeric"
    if(class(changes)!="list"){
    	ch2<-vector(mode="list",length=length(changes))
    	for(q in 1:length(ch2)){
    		ch2[[q]]<-changes[q]
    	}
    	changes<-ch2
    } 
    #if wanting to focus on specific changes, keep changes from the changes list that are the desired changes
    if(!is.null(direction)){
		#make an empty changes list
		changes.mod <- vector(mode = "list", length = length(changes))
		for(chng in 1:length(changes)){
			changes.mod[[chng]]<-changes[[chng]][names(changes[[chng]]) %in% direction]
		}
		changes<-changes.mod	    
    }
    h <- max(nodeHeights(tree))
    b <- ceiling(h/timebin)
    segs <- cbind(rev(c(seq(h,0,-timebin),0))[1:b], rev(c(seq(h,0,-timebin),0))[1:b+1])
    #first cell needs to be partial because 0 is actually root in this situation
    #segs <- cbind(seq(0, h, timebin), seq(timebin, timebin*b, timebin))
    #change end of max segment to max height of tree
    #segs[nrow(segs),2]<-h    
    nchanges <- rep(0, b)
    for (i in 1:length(changes)) {
        for (j in 1:length(changes[[i]])) {
            ind <- which((changes[[i]][j] > segs[, 1]) + (changes[[i]][j] <= 
                segs[, 2]) == 2)
            nchanges[ind] <- nchanges[ind] + 1/length(changes)
        }
    }
    LTT <- ltt(tree, plot = FALSE)
    LTT <- cbind(LTT$ltt[2:(length(LTT$ltt) - 1)], LTT$times[2:(length(LTT$ltt) - 
        1)], LTT$times[3:length(LTT$ltt)])
    ii <- 1
    edge.length <- rep(0, b)
    for (i in 1:nrow(segs)) {
        done.seg <- FALSE
        while (LTT[ii, 2] <= segs[i, 2] && done.seg == FALSE) {
            edge.length[i] <- edge.length[i] + LTT[ii, 1] * (min(segs[i, 
                2], LTT[ii, 3]) - max(segs[i, 1], LTT[ii, 2]))
            if (LTT[ii, 3] >= segs[i, 2]) 
                done.seg <- TRUE
            if (LTT[ii, 3] <= segs[i, 2]) 
                ii <- if (ii < nrow(LTT)) 
                  ii + 1
                else ii
        }
    }
    object <- list(segments = segs, nchanges = nchanges, edge.length = edge.length, 
        tree = tree)
    class(object) <- "ctt"
    object
}
