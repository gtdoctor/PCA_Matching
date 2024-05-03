# Ancestry matching algorithm first developedby AP Levine.
# updated 2024 by GT Doctor


library(ggplot2)
library(parallel)


###### CHOICES #######

# PATHS # 
DATADIR="my/pca/data"
OUTDIR="my/output/dir" # output directory. 

# INPUTS #
PCA.eigenval="casectrl_forPCA.eigenval" # headerless list of eigenvals 1 per line. 
PCA.eigenvec="casectrl_forPCA.eigenvec" # plink2 --PCA output, headed, with FID, IID, PCAs
PSAM="plink.psam/orfam/" # plinkfile: headerless table with FID, IID columns (col1 is ignored entirely); or headed plink2 psam

# OUTPUTS #

FILEOUTSTEM= "distancematrix" # save the R workspace image. euc/eig/wtd will be appended
screename="Scree plot of principal components"    

# Parameters #
plinkversion="v2" ## choose "v1.9" or "v2" if v2 expects a column called Cohort - might be PHENO1
smallsample <- "n" # subsets data to a small proportion of cases and controls for testing
euclidean <- "y" # default manhattan. apply "y" to calculate unweighted euclidean distances  
eigvalweight <- "y" # should distances be weighted by the eigvalue of the PC.
nchunks = 40 # number of case chunks for parallel processing.

PCmax=10 # if graphing is to be done, can only take value >=8 at the moment. 

  
######### CODE ############
an = function(x){as.numeric(as.character(x))}
setwd(DATADIR)
FILEOUTSTEM <- paste0(FILEOUTSTEM,"_",PCmax,"pcs")

#Load PC data
d=read.table(PCA.eigenvec,header=T, comment.char = "")
if (plinkversion=="v1.9"){
names(d)=c("FID","IID",paste("PC",1:10,sep=""))
}

#Load weightings
ev = an(read.table(PCA.eigenval, header=F)[,1])
ev.weight = ev/sum(ev) # weighting as proportion of top 10 eigenvals


# read case control. check here that the correct columns (e.g. XIID, PC1-10 are being selected)
if (plinkversion == "v1.9") {
cohort=read.table(PSAM, header=F)
cases <- d[d$IID %in% cohort$V2[cohort$V6 == 2],2:12]
controls  <- d[d$IID %in% cohort$V2[cohort$V6 == 1],2:12]
}
if (plinkversion =="v2"){
cohort=read.table(PSAM, header = T, comment.char = "", sep = "\t")
cases <- d[d$IID %in% cohort$IID[cohort$PHENO1 == 2],2:12] # check column names
controls <- d[d$IID %in% cohort$IID[cohort$PHENO1 == 1],2:12]
}

# small sample
if (smallsample == "y") {
  cases= cases[1:50,]
  controls=controls[1:2000,]
  controlnames=c(paste0("C",1:nrow(controls)))
  controls$IID=controlnames
}


## Scree Plot - to decide on PCs needed
scree <- data.frame(PC = 1:length(ev.weight), VarianceExplained = ev.weight * 100)

# Plotting the scree plot
ggplot(scree, aes(x = factor(PC), y = VarianceExplained)) +
  geom_bar(stat = "identity") +
  xlab("Principal Component") +
  ylab("Variance Explained (%)") +
  ggtitle("Scree Plot") +
  ylim(0, 100)

# Create case and control matrices
cases_matrix <- as.matrix(cases[, 2:(PCmax+1)])  # Exclude ID columns
controls_matrix <- as.matrix(controls[, 2:(PCmax+1)])  # Exclude ID columns
ncases=nrow(cases_matrix)
nctrls=nrow(controls_matrix)

# Initialize a matrix to store the distances
ccdistances <- matrix(nrow = ncases, ncol = nctrls)


# Split cases_matrix into chunks for parallel
chunk_size = ceiling(ncases / nchunks)
casechunks = lapply(1:nchunks, function(i) {
  start_row = (i - 1) * chunk_size + 1
  end_row = min(i * chunk_size, ncases)
  cases_matrix[start_row:end_row, ]
})

## Distance calculations

# Compute weighted Manhattan distances
if (euclidean == "n") {
# Create  (weighted) distances matrices
if (eigvalweight == "y") {
  filesave=paste0(FILEOUTSTEM,"_manwtd.Rdata")
  ev.weight = ev.weight[1:PCmax]

  print("Computing weighted Manhattan distances")

  distancefunction = function(casechunk, nctrls, ctrls_index, controls_matrix, ev.weight) {
    # Initialize matrix to store results for this chunk with appropriate dimensions
    chunk_results = matrix(0, nrow(casechunk), nctrls)
    
    for (i in 1:nrow(casechunk)) {
      for (j in 1:nctrls) {
        # Calculate abs weighted difference for all PCs between case i and control j
        abs_diffs <- abs(casechunk[i, ] - controls_matrix[j, ]) * ev.weight
        chunk_results[i, j] <- sum(abs_diffs)
      }
    }
    chunk_results
  }
  results = mclapply(casechunks, distancefunction, nctrls = nctrls, controls_matrix = controls_matrix, ev.weight = ev.weight, mc.cores = detectCores())
  ccdistances = do.call(rbind, results)
  
  } # end Man weighted 

if (eigvalweight=="n"){
  filesave=paste0(FILEOUTSTEM,"_manunwtd.Rdata")  
# Compute uneighted Manhattan distances
print("Computing absolute differences")


distancefunction = function(casechunk, nctrls, ctrls_index, controls_matrix) {
  # Initialize matrix to store results for this chunk with appropriate dimensions
  chunk_results = matrix(0, nrow(casechunk), nctrls)
  
  for (i in 1:nrow(casechunk)) {
    for (j in 1:nctrls) {
      # Calculate abs  difference for all PCs between case i and control j
      abs_diffs <- abs(casechunk[i, ] - controls_matrix[j, ]) 
      chunk_results[i, j] <- sum(abs_diffs)
    }
  }
  chunk_results
}
results = mclapply(casechunks, distancefunction, nctrls = nctrls, controls_matrix = controls_matrix,  mc.cores = detectCores())
ccdistances = do.call(rbind, results)


} # End man weighted
} # end euclidean == "n"

if (euclidean == "y") {  

if (eigvalweight == "n") {
  print("Computing unweighted Euclidean distances")
    filesave=paste0(FILEOUTSTEM,"_eucunwtd.Rdata")
    

for (i in 1:nrow(cases_matrix)) {
  for (j in 1:nrow(controls_matrix)) {
    # Calculate squared difference for all PCs between case i and control j
    squared_distances <- (cases_matrix[i, ] - controls_matrix[j, ])^2
    # Sum and take the square root to get weighted Euclidean distance
    ccdistances[i, j] <- sqrt(sum(squared_distances))
  }
  }


    distancefunction = function(casechunk, nctrls, ctrls_index, controls_matrix, ) {
      # Initialize matrix to store results for this chunk with appropriate dimensions
      chunk_results = matrix(0, nrow(casechunk), nctrls)
      for (i in 1:nrow(casechunk)) {
        for (j in 1:nctrls) {
          # Calculate squared difference for all PCs between case i and control j
          squared_distances <- (casechunk[i, ] - controls_matrix[j, ])^2 # squared_distances is a vector as is ev weight; multiples matching elements
          # Sum and take the square root to get weighted Euclidean distance
          chunk_results[i, j] <- sqrt(sum(squared_distances))
        }
      }
      chunk_results
    }
    results = mclapply(casechunks, distancefunction, nctrls = nctrls, controls_matrix = controls_matrix, mc.cores = detectCores())
    ccdistances = do.call(rbind, results)
    
    } # end euc unweighted

if (eigvalweight =="y") {
  print("Computing weighted Euclidean distances")
  filesave=paste0(FILEOUTSTEM,"_eucwtd.Rdata")
  ev.weight = ev.weight[1:PCmax]
  
    distancefunction = function(casechunk, nctrls, ctrls_index, controls_matrix,ev.weight) {
    # Initialize matrix to store results for this chunk with appropriate dimensions
    chunk_results = matrix(0, nrow(casechunk), nctrls)
    for (i in 1:nrow(casechunk)) {
      for (j in 1:nctrls) {
        # Calculate squared weighted difference for all PCs between case i and control j
        squared_distances <- (casechunk[i, ] - controls_matrix[j, ])^2 * ev.weight # squared_distances is a vector as is ev weight; multiples matching elements
        # Sum and take the square root to get weighted Euclidean distance
        chunk_results[i, j] <- sqrt(sum(squared_distances))
      }
    }
    chunk_results
  }
  results = mclapply(casechunks, distancefunction, nctrls = nctrls, controls_matrix = controls_matrix, ev.weight = ev.weight, mc.cores = detectCores())
  ccdistances = do.call(rbind, results)
}  
}  # end euclidean == "y"


setwd(OUTDIR)
ggsave(plot=last_plot(), filename = paste0(screename, ".jpg"), device = "jpeg")
print(paste0("Saving ",filesave))
save(ccdistances, file=filesave)
