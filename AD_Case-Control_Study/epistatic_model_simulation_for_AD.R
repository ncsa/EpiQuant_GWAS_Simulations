rm(list = ls())
#########################################

#This function will create simulated data for a trait with genes that following a geometric series.
# The user will input the number of QTL, the heritabilities to be explored, the additive effect,
# and the output directory

create.simluated.data <- function(){
 
  #Filter out SNPs of the 5K Subset that are monomorphic
  #genotypes.5K.subset.n.equals.300 <-genotypes.5K.subset[,c(1:5,which(colnames(genotypes.5K.subset) %in% as.character(simulated.data.300.subset.of.inds$`<Trait>`)))]
  
  #genotypes.5K.subset.only <- genotypes.5K.subset.n.equals.300[,-c(1:5)]
  
  #Total number of lines
  #ns <- ncol(genotypes.5K.subset.only)
  
  #Sum of the allele scores for each SNP
  #ss <- apply(genotypes.5K.subset.only, 1, sum)
  
  #Combine two situations: one where the allele coded as "2" is major; one where "0" is coded as major.
  #maf.matrix <- rbind((.5*ss/ns), (1-(0.5*ss/ns)))
  
  #Copy the minor allele frequencies for all SNPs
  #maf <- apply(maf.matrix, 2, min)
  
  #Find out which SNPs have MAF < 0.05
  #non.monomorphic.SNPs <- which(maf != 0)
  
  
  #genotypes.non.monomorphic.5K.subset <- genotypes.5K.subset[non.monomorphic.SNPs,]
  
  #Create a working directory for the output results:
  dir.create(paste(home.dir,"/", output.dir, sep = ""))
  
  #Set the working directory
  setwd(paste(home.dir,"/", output.dir, sep = ""))
  
  #Randomly select (without replacement) k additive QTN, and assign an effect size
  seed.number <- sample(-1000000:1000000, 1)
  #seed.number <- 130164
  
  # Use the seed number to generate uniform random variables
  set.seed(seed.number)
  
  vector.of.add.QTN <- sample(2:ncol(genotypes.for.n.2099.5000.marker.subset), Additive.QTN.number, replace = FALSE)
  Add.QTN.genotypic.information <- genotypes.for.n.2099.5000.marker.subset[,c(1,vector.of.add.QTN)]
  
  
  #Create an output file that gives the chromosome, bp, and additive effect of the effect sizes, as well as the seed
  # number used
  write.table(paste("Here_is_the_seed_number:_",seed.number, sep = ""), paste("Seed.number.for.", Additive.QTN.number,"Add.QTN",
                                                                              ".txt", sep = ""), row.names = FALSE, col.names = FALSE, sep = "\t",  quote = FALSE)
  write.table(Add.QTN.genotypic.information, paste("Genotypic.information.for.", Additive.QTN.number,".Additive.QTN",
                                                   ".txt", sep = ""), row.names = FALSE, sep = "\t",  quote = FALSE)
  
  
  #Randomly select (without replacement) 2*k epistatic QTN, and assign an effect size
  seed.number <- sample(-1000000:1000000, 1)
  #seed.number <- 130164
  
  # Use the seed number to generate uniform random variables
  set.seed(seed.number)
  
  vector.of.epi.QTN <- sample(2:ncol(genotypes.for.n.2099.5000.marker.subset), (2*Epistatic.QTN.number), replace = FALSE)
  Epi.QTN.genotypic.information <- genotypes.for.n.2099.5000.marker.subset[,c(1,vector.of.epi.QTN)]
  
  
  #Create an output file that gives the chromosome, bp, and additive effect of the effect sizes, as well as the seed
  # number used
  write.table(paste("Here_is_the_seed_number:_",seed.number, sep = ""), paste("Seed.number.for.", Epistatic.QTN.number,"Epi.QTN",
                                                                              ".txt", sep = ""), row.names = FALSE, col.names = FALSE, sep = "\t",  quote = FALSE)
  write.table(Epi.QTN.genotypic.information, paste("Genotypic.information.for.", Epistatic.QTN.number,".Epistatic.QTN",
                                                   ".txt", sep = ""), row.names = FALSE, sep = "\t",  quote = FALSE)
  
  
  #Create a "base line" trait, which is basically just the additive effects; this is what we would see if 
  # the heritability of the simulated trait were 1
  additive.effect.trait.object <- Add.QTN.genotypic.information[,-1] #this was originally the base.line.trait.object
  
  epistatic.effect.trait.object <-Epi.QTN.genotypic.information[,-1]
  #AEL Changed: - epistatic.effect.trait.object<- epistatic.effect.trait.object[,number.of.epistasis]
  
  
  
  #make base.line.trait additive.component and epistatic.component
  additive.component<- try(as.data.frame(matrix(0, nrow = nrow(additive.effect.trait.object), ncol = 1)))
  if(inherits(additive.component, "try-error")){
    additive.component <- as.data.frame(matrix(0, nrow = length(additive.effect.trait.object), ncol = 1))
  }
  epistatic.component <- try(as.data.frame(matrix(0, nrow = nrow(epistatic.effect.trait.object), ncol = 1)))
  if(inherits(epistatic.component, "try-error")){
    epistatic.component <- as.data.frame(matrix(0, nrow = length(epistatic.effect.trait.object), ncol = 1))
  }
  #base.line.trait <- as.data.frame(matrix(0, nrow = nrow(base.line.trait.object), ncol = 1)) 
  if(Additive.QTN.number > 1){
    for(i in 1:Additive.QTN.number) additive.component <- additive.component + (additive.effect.trait.object[,i]*(additive.effect^i))
  }else{
    additive.component <- additive.component + (additive.effect.trait.object*(additive.effect^i))
  }
  
  rownames(additive.component) <- rownames(additive.effect.trait.object)
  colnames(additive.component) <- "Additive.effect"
  additive.genetic.variance <- var(additive.component)
  
  last.number.of.this.loop <- Epistatic.QTN.number - 1
  for(i in 0:last.number.of.this.loop) epistatic.component <- epistatic.component + ((epistatic.effect.trait.object[,((2*i)+1)]*epistatic.effect.trait.object[,((2*i)+2)])*(epistatic.effect^(i+1)))
  rownames(epistatic.component) <- rownames(epistatic.effect.trait.object)
  colnames(epistatic.component) <- "Epistatic.effect"
  epistatic.genetic.variance<- var(epistatic.component)
  
  #Set the row names of the base.line.trait object to the new names
  base.line.trait <- additive.component+epistatic.component
  
  base.line.trait.temp <- data.frame(genotypes.for.n.2099.5000.marker.subset[,1],base.line.trait)
  colnames(base.line.trait.temp)[1] <- "<Trait>"
  colnames(base.line.trait.temp)[2] <- "the.normal.random.variables.1"
  
  #write.table(base.line.trait.temp, "Baseline_Trait_n=2099.txt", row.names = FALSE, sep = "\t",  quote = FALSE) 
  #base.line.trait.with.new.taxa <- merge(base.line.trait, taxa.name.converter, by.x = "row.names", 
  #                                       by.y = "Old_Taxa_ID")
  
  #the.new.taxa.ids <- as.character(base.line.trait.with.new.taxa[,2])
  #base.line.trait <- as.matrix(base.line.trait.with.new.taxa[,2], nrow = nrow(base.line.trait.with.new.taxa))
 # rownames(base.line.trait) <- as.character(base.line.trait[,3])
  
  #Alex: begin here.
  #For loop through the vector of heritabilities
  for(i in heritabilities.vector){
    #If heritability is zero
    if(i == 0){ 
      #Simulate m replicates of n N(0,b) random variables, where b = additive.genetic.variance
      the.seed.number.vector <- NULL
      for(j in 1:replicates){
        seed.number <- sample(-1000000:1000000, 1)
        set.seed(seed.number)
        the.normal.random.variables <- rnorm(nrow(base.line.trait), mean = 0, sd = 1)
        if(j == 1){
          simulated.data <- the.normal.random.variables
        }else{
          simulated.data <- cbind(simulated.data, the.normal.random.variables)
          colnames(simulated.data)[j] <- paste(colnames(simulated.data)[j],".",j,sep = "")
        }
        the.seed.number.vector <- c(the.seed.number.vector, seed.number)
      }
      
      #Format the output file for the simulated phenotypes
      simulated.data <- data.frame(as.factor(Add.QTN.genotypic.information[,1]),simulated.data)
      colnames(simulated.data)[1] <- "<Trait>"
      colnames(simulated.data)[2] <- "the.normal.random.variables.1"
      
      
  
      #Output the m replicates and the seed numbers, formatted for TASSEL
      write.table(the.seed.number.vector, paste("Seed.number.for.", replicates,".Reps",
                                                ".Herit.",i,"n-2099.txt", sep = ""), row.names = FALSE, col.names = FALSE, sep = "\t",  quote = FALSE)
      write.table(simulated.data, paste("Simulated.Data.", replicates,".Reps",
                                        ".Herit.",i,"n=2099.txt", sep = ""), row.names = FALSE, sep = "\t",  quote = FALSE)
      
     
      #Take a subset of 300 of these individuals
      #seed.number.for.300.individuals <- sample(-1000000:1000000, 1) #On August 24, 2017 at 4:53 pm I ran this and obtained a seed number of -434778
      seed.number.for.300.individuals <- 680144
      set.seed(seed.number.for.300.individuals)
      
      vector.of.row.nums.for.300.subset <- sample(1:nrow(simulated.data), 300, replace = FALSE)
      simulated.data.300.subset.of.inds <- simulated.data[vector.of.row.nums.for.300.subset,]
      
      write.table(simulated.data.300.subset.of.inds , paste("Simulated.Data.", replicates,".Reps",
                                                            ".Herit.",i,"n=300.txt", sep = ""), row.names = FALSE, sep = "\t",  quote = FALSE)  
      
    }else{
      #Calcualte V_e, the residual variance
      #residual.variance <- (additive.genetic.variance*(1-i))/i
      residual.variance <-( ((additive.genetic.variance+epistatic.genetic.variance)/i) - additive.genetic.variance - epistatic.genetic.variance)
      total.variance <- additive.genetic.variance+epistatic.genetic.variance+residual.variance
      write.table(total.variance, "Total.variance.2099.inds.txt", row.names = FALSE, sep = "\t",  quote = FALSE)       
      #Use V_e to generate n replicates of N(0, Ve) random variables
      write.table(base.line.trait.temp, "Baseline_Trait_n=2099.txt", row.names = FALSE, sep = "\t",  quote = FALSE) 
      
      the.seed.number.vector <- NULL
      col.name.vector <- NULL
      for(j in 1:replicates){
        seed.number <- sample(-1000000:1000000, 1)
        set.seed(seed.number)
        the.normal.random.variables <- rnorm(nrow(base.line.trait), mean = 0, sd = sqrt(residual.variance))
        the.base.line.trait.plus.normal.random.variables <- base.line.trait+the.normal.random.variables
        if(j == 1){
          simulated.data <- the.base.line.trait.plus.normal.random.variables
        }else{
          simulated.data <- cbind(simulated.data, the.base.line.trait.plus.normal.random.variables)
          colnames(simulated.data)[j] <- paste(colnames(simulated.data)[j],".",j,sep = "")
        }
        the.seed.number.vector <- c(the.seed.number.vector, seed.number)
        col.name.vector <- c(col.name.vector, paste("Heritability_",i, "_Rep_", j, sep = ""))
      }
      
      colnames(simulated.data)  <- col.name.vector
      
      #Format the output file for the simulated phenotypes
      simulated.data <- data.frame(as.factor(Add.QTN.genotypic.information[,1]),simulated.data)
      colnames(simulated.data)[1] <- "<Trait>"
      
      #Output the m replicates and the seed numbers, formatted for TASSEL
      write.table(the.seed.number.vector, paste("Seed.number.for.", replicates,".Reps",
                                                ".Herit.",i,".txt", sep = ""), row.names = FALSE, col.names = FALSE, sep = "\t",  quote = FALSE)
      write.table(simulated.data, paste("Simulated.Data.", replicates,".Reps",
                                        ".Herit.",i,"n=2648.txt", sep = ""), row.names = FALSE, sep = "\t",  quote = FALSE)
      
    
      
      #Take a subset of 300 of these individuals
      
      
      
      #seed.number.for.300.individuals <- sample(-1000000:1000000, 1) #On August 24, 2017 at 4:53 pm I ran this and obtained a seed number of -434778
     # seed.number.for.300.individuals <- -434778
      seed.number.for.300.individuals <- 680144
      set.seed(seed.number.for.300.individuals)
      
      
      vector.of.row.nums.for.300.subset <- sample(1:nrow(simulated.data), 300, replace = FALSE)
      #Take the observed "trait" values for the additive, epistatic, and "base line trait" vlaues
      
      additive.component.300.inds <- data.frame(additive.component[vector.of.row.nums.for.300.subset,]) 
      rownames(additive.component.300.inds) <- rownames(additive.component)[vector.of.row.nums.for.300.subset]                                            
      additive.genetic.variance.300.inds <- var(additive.component.300.inds)
      
      epistatic.component.300.inds <- data.frame(epistatic.component[vector.of.row.nums.for.300.subset,]) 
      rownames(epistatic.component.300.inds) <- rownames(epistatic.component)[vector.of.row.nums.for.300.subset]                                            
      epistatic.genetic.variance.300.inds <- var(epistatic.component.300.inds)      
      
      base.line.trait.300.inds <- additive.component.300.inds+epistatic.component.300.inds
      
      base.line.trait.300.inds.temp <- data.frame(as.factor(Add.QTN.genotypic.information[vector.of.row.nums.for.300.subset,1]),base.line.trait.300.inds)
      colnames(base.line.trait.300.inds.temp)[1] <- "<Trait>"
      colnames(base.line.trait.300.inds.temp)[2] <- "the.normal.random.variables.1"
      
      write.table(base.line.trait.300.inds.temp, "Baseline_Trait_n=300.txt", row.names = FALSE, sep = "\t",  quote = FALSE) 
      
      residual.variance.300.inds <-( ((additive.genetic.variance.300.inds+epistatic.genetic.variance.300.inds)/i) - additive.genetic.variance.300.inds - epistatic.genetic.variance.300.inds)
      
      total.variance.300.inds <- additive.genetic.variance.300.inds+epistatic.genetic.variance.300.inds
      write.table(total.variance.300.inds, "Total.variance.300.inds.txt", row.names = FALSE, sep = "\t",  quote = FALSE) 
      
      
      
      the.seed.number.vector.300.subset <- NULL
      col.name.vector <- NULL
      for(j in 1:replicates){
        seed.number <- sample(-1000000:1000000, 1)
        set.seed(seed.number)
        the.normal.random.variables.300.inds <- rnorm(nrow( base.line.trait.300.inds), mean = 0, sd = sqrt(residual.variance.300.inds))
        the.base.line.trait.plus.normal.random.variables.300.inds <-  base.line.trait.300.inds+the.normal.random.variables.300.inds 
        if(j == 1){
          simulated.data.300.subset.of.inds <- the.base.line.trait.plus.normal.random.variables.300.inds
        }else{
          simulated.data.300.subset.of.inds <- cbind(simulated.data.300.subset.of.inds, the.base.line.trait.plus.normal.random.variables.300.inds)
          colnames(simulated.data.300.subset.of.inds)[j] <- paste(colnames(simulated.data.300.subset.of.inds)[j],".",j,sep = "")
        }
        the.seed.number.vector.300.subset <- c(the.seed.number.vector.300.subset, seed.number)
        col.name.vector <- c(col.name.vector, paste("Heritability_",i, "_Rep_", j, sep = ""))
      }#end for(j in 1:replicates)
      
      colnames(simulated.data.300.subset.of.inds)  <- col.name.vector
      
      #Format the output file for the simulated phenotypes
      simulated.data.300.subset.of.inds <- data.frame(base.line.trait.300.inds.temp[,1], simulated.data.300.subset.of.inds)
      colnames(simulated.data.300.subset.of.inds)[1] <- "<Trait>"
      
      
      
      write.table(the.seed.number.vector.300.subset, paste("Seed.number.for.", replicates,".Reps",
                                                           ".Herit.",i,"n=300.txt", sep = ""), row.names = FALSE, col.names = FALSE, sep = "\t",  quote = FALSE)
      write.table(simulated.data.300.subset.of.inds, paste("Simulated.Data.", replicates,".Reps",
                                                           ".Herit.",i,"n=300.txt", sep = ""), row.names = FALSE, sep = "\t",  quote = FALSE) 
      
      #vector.of.row.nums.for.300.subset <- sample(1:nrow(simulated.data.1K.subset.of.inds), 300, replace = FALSE)
      #simulated.data.300.subset.of.inds <- simulated.data.1K.subset.of.inds[vector.of.row.nums.for.300.subset,]
      
      #write.table(simulated.data.300.subset.of.inds , paste("Simulated.Data.", replicates,".Reps",
      #                                                    ".Herit.",i,"n=300.txt", sep = ""), row.names = FALSE, sep = "\t",  quote = FALSE) 
      
    }   
  }#End for(i in heritabilities.vector)
  
}#end "create.simluated.data()"





###########################################################################################
###########################################################################################
###########################################################################################
#setwd("/Users/adminuser/Box Sync/Lipka_Mainzer_Chen_Epistasis_Shared_Folder/Simulation_Study")
setwd("/Users/adminuser/Box Sync/Lipka_Mainzer_Chen_Epistasis_Shared_Folder/Simulation_Study/AD/")
home.dir <- getwd()
#dir.of.GBS.SNPs <- "/Users/adminuser/Desktop/Work/Tocos_NAM_2009_2010/Joint_Linkage_Analysis/GBS_SNPs/"
#setwd("NCRPIS_Ames_Genotypes")
load("All.AD.Genotypes.and.Markers.Rdata")


#Below is some work I did to i.) output the 50K SNPs and ii.) obtain a subset of 5,000 SNPs
#setwd("NCRPIS_Ames_Genotypes")

#i.) output the 50K SNPs and the seed number used to generate this particular subset
#write.table(genotypes, "Subset.of.50K.SNPs.txt", row.names = FALSE, sep = "\t",  quote = FALSE)
#write.table(vector.of.seed.numbers, "Seed.Numbers.For.Each.Chr.For.Obtaining.50K.SNPs.txt", row.names = FALSE, sep = "\t",  quote = FALSE)

#ii.) obtain a subset of 5,000 SNP
#seed.number <- sample(-1000000:1000000, 1) #I ran this line of code on 8/24 at 1:02 pm Central time, and I obtained a value of -984153
#seed.number.for.generating.5K.subset <- -984153
#set.seed(seed.number)

#vector.of.row.IDs.for.5K.subset <- sample(1:nrow(genotypes), 5000, replace = FALSE)
#genotypes.5K.subset <- genotypes[vector.of.row.IDs.for.5K.subset,]

#setwd(home.dir)
#save.image("50K.and.5K.GBS.SNPs.for.Ames.2648.Rdata")

#write.table(genotypes.5K.subset, "Subset.of.5K.SNPs.txt", row.names = FALSE, sep = "\t",  quote = FALSE)

#setwd(home.dir)
#genotypes <- read.table("4K_SNPsmdp_genotype_test1.hmp.txt", head = TRUE)
#setwd(home.dir)
library(rrBLUP)
library('MASS')
#library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(EMMREML)
library(compiler) #this library is already installed in R
library("scatterplot3d")
source("http://zzlab.net/GAPIT/previous/gapit_functions20160408.txt")
#source("GAPIT_Code_from_Internet_20120411_Allelic_Effect.R")
source("http://zzlab.net/GAPIT/emma.txt")
#library(rrBLUP)


###############
#User input

#Number of additive QTN (k)
Additive.QTN.number <- 4


#Number of epistatic QTN (m)
Epistatic.QTN.number <- 4

#Vector of heritabilities to investigate
heritabilities.vector <- 0.99

#Size of the additive effect of the largest QTL (must be (-1,1) but preferably (0,1))
additive.effect <- 0.95


#Size of the epistatic effect of the largest QTL (must be (-1,1) but preferably (0,1)|
epistatic.effect <- 0.95

#Number of replicates of simulated phenotypes for each heritability (m)
replicates <- 100


#file.G="NCRPIS_n=2648_Chr_" 
#file.Ext.G = "hmp.txt"


#file.from = 9
#file.to = 10
#file.fragment = 1000000


#Output directory
output.dir <- paste(Additive.QTN.number,"_Add_QTN",Epistatic.QTN.number,"_Epi_QTN_h.2_",
                    heritabilities.vector,"_add.eff_", additive.effect,"_epis.eff_", epistatic.effect,"_reps_", replicates, sep = "")

################
#Create the simulated data
create.simluated.data()


#################################################################################################################################
#Old code, which can be discarded
#Below is some code for reformatting the 5K and 50K marker subsets so that they can be efficiently read into TASSEL
setwd("NCRPIS_Ames_Genotypes")
View(genotypes.5K.subset)
#Reforamt the 5K subset of markers so that it can be read into TASSEL
genotypes.5K.subset.for.TASSEL.temp.1 <- genotypes.5K.subset[,-c(1:5)]/2
genotypes.5K.subset.for.TASSEL.temp.2 <- data.frame(genotypes.5K.subset[,1], genotypes.5K.subset.for.TASSEL.temp.1)

genotypes.5K.subset.for.TASSEL <- t(genotypes.5K.subset.for.TASSEL.temp.2)

colnames(genotypes.5K.subset.for.TASSEL) <- genotypes.5K.subset.for.TASSEL[1,] 

genotypes.5K.subset.for.TASSEL <- genotypes.5K.subset.for.TASSEL[-1,]

genotypes.5K.subset.for.TASSEL <- data.frame(colnames(genotypes.5K.subset[-c(1:5)]), genotypes.5K.subset.for.TASSEL)


colnames(genotypes.5K.subset.for.TASSEL)[1] <- "<Marker>"

write.table("<Numeric>" , "5K.Marker.Subset.for.TASSEL.txt", col.names = FALSE,row.names = FALSE, sep = "\t",  quote = FALSE) 
write.table(genotypes.5K.subset.for.TASSEL, "5K.Marker.Subset.for.TASSEL.txt", row.names = FALSE, sep = "\t",  quote = FALSE, append = TRUE)


#Reforamt the 50K subset of markers so that it can be read into TASSEL
genotypes.for.TASSEL.temp.1 <- genotypes[,-c(1:5)]/2
genotypes.for.TASSEL.temp.2 <- data.frame(genotypes[,1], genotypes.for.TASSEL.temp.1)

genotypes.for.TASSEL <- t(genotypes.for.TASSEL.temp.2)

colnames(genotypes.for.TASSEL) <- genotypes.for.TASSEL[1,] 

genotypes.for.TASSEL <- data.frame(rownames(genotypes.for.TASSEL), genotypes.for.TASSEL)

genotypes.for.TASSEL <- genotypes.for.TASSEL[-1,]

colnames(genotypes.for.TASSEL)[1] <- "<Marker>"

write.table("<Numeric>" , "50K.Marker.Subset.for.TASSEL.txt", col.names = FALSE,row.names = FALSE, sep = "\t",  quote = FALSE) 
write.table(genotypes.for.TASSEL, "50K.Marker.Subset.for.TASSEL.txt", row.names = FALSE, sep = "\t",  quote = FALSE, append = TRUE)


