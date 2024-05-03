# Ancestry matching algorithm first developedby AP Levine.
# updated 2024 by GT Doctor g.doctor@ucl.ac.uk UCL Centre for Genetics and Genomics

library(ggplot2)
library(grid)
library(gridExtra)


# loading the distmatrix #
setwd("/home/gabriel/nasgab/INS/Dec2023/Merged/Ancestry_match16042024")
shouldload <- "y" # if distance matrix should be loaded (if matrix part 1 run in seperate session)
### Loading matrices #####
RWORKSPACE="pathto/distancematrix.Rdata" ### distance matrices, case cohort PC matrices

if (shouldload == "y") {load(file=RWORKSPACE)}

###### CHOICES #######

# INPUTS #
COHORT <-  "/mnt/nas1/projects/gabriel/INS/Dec2023/Merged/mergedmegaUkbforpca.23feb.psam"
PCA.eigenvec <- "casectrl_forPCA.eigenvec" # plink2 --PCA output, headed, with FID, IID, PCAs

## CHECK if the dataframes have ID as X.IID or IID - may need to change throughout. 
plinkversion <- "v2"

# OUTPUTS #
FILEOUTSTEM="megains_ukb1602024" 

# PATHS # 
DATADIR=file.path("/mnt/nas1/projects/gabriel/INS/Dec2023/Merged/")
OUTDIRP=file.path("/mnt/nas1/projects/gabriel/INS/Dec2023/Merged/Ancestry_match16042024/") # output directory parent. No final slash 



# PARAMETERS #
threshcutoff = 0.3 # number 0-1 # smaller == tighter
nc = 7 #  max number of controls that permitted per case within distance t for case to be included
ncmin= 7  # <=nc case excluded if fewer controls found. 

######### CODE ############
setwd(DATADIR)
COHORT  <-  read.table(COHORT, header = T, comment.char = "")
PCA.eigenvec <- read.table(PCA.eigenvec, header=T, comment.char = "")

setwd(OUTDIRP)
    
text_choice=paste0("Max controls = ", nc, ", Min controls = ", ncmin, ", Thresholding = ", threshcutoff, " of median")

fileoutkept=paste0("kept.", FILEOUTSTEM)
fileoutremoved=paste0("removed.",FILEOUTSTEM)
graphstem <- paste0("plot.", FILEOUTSTEM)


ac = function(x){as.character(x)}
an = function(x){as.numeric(ac(x))}



######------------- Operations on distance matrix ------- #########


# calculate cutoff based on median of all case-control differences
median_distance <- median(ccdistances, na.rm = TRUE)
threshold=threshcutoff * median_distance

# create matrix of ids
if (plinkversion=="1.9") {
case_ids <- cases$IID  # Assuming the case IDs are stored in 'IID' column of 'cases_use'
control_ids <- controls$IID
}
if (plinkversion=="v2") {
  case_ids <- COHORT$IID[COHORT$PHENO1==2]  # Assuming the case IDs are stored in 'IID' column of 'cases_use'
  control_ids <- COHORT$IID[COHORT$PHENO1==1]  # Assuming the case IDs are stored in 'IID' column of 'cases_use'
}

cases <- PCA.eigenvec[PCA.eigenvec$IID %in% case_ids,]
controls <- PCA.eigenvec[PCA.eigenvec$IID %in% control_ids,]
  
# Create a list with empty character vectors named by index i

matched_casescontrols <- vector("list", length(case_ids))
names(matched_casescontrols) <- case_ids

# Initialize lists for keeping and removing
cases_kept <- character()
cases_removed <- character()
controls_kept <- character()  # If using lists, initialize as vector() and then convert to unique values at the end
controls_removed <- character()  # Similarly, initialize appropriately

# Create "shortest_distances" matrix with i rows and three cols (d, j, count)
shortest_distances <- matrix(NA, nrow = nrow(ccdistances), ncol = 3)
colnames(shortest_distances) <- c("d", "j", "count")

# Create a thresholded version of the ccdistances matrix: 
# Replace NA for any D > threshold in ccdistances:

wdt = ccdistances
wdt[wdt > threshold] <- NA

# Assign NA to cases where there aren't enough controls left.
for (i in 1:nrow(wdt)) {
  if (sum(!is.na(wdt[i, ])) < ncmin) {
    wdt[i, ] <- NA
  }
}

initial_cases_removed = nrow(wdt) - sum(apply(!is.na(wdt), 1, any))

# For each case i, find smallest D and its j index
min_d = numeric()
j_index = numeric()
index_largest_d = numeric()
j_value = numeric()
replace_i = numeric()
j.rmNA=numeric() # a caddy to update all newly found controls with NA in the main matrix

for (i in 1:nrow(wdt)) {
  # Find the minimum distance and its index if there are non-NA values
     j_index <- which.min(wdt[i, ])
     if (length(j_index)>0) {
    min_d <- wdt[i, j_index]
    j.rmNA <- c(j.rmNA, j_index)
    shortest_distances[i, ] <- c(min_d, j_index, 0) # update the shortest_distances matrix
  } else {
    # Assign NA to d and j but keep the row consistent - ie no matching controls at the start. 
    shortest_distances[i, ] <- c(NA, NA, 0) # Assign 0 or NA to 'count' 
  }
}

# block out used controls from 1st round. 
wdt[,unique(j.rmNA)] <- NA 



# Loop until conditions are met
repeat {
  
  # find the case with the largest smallest distance
  index_largest_d <- which.max(shortest_distances[, 1])  
  
  
  # Check if the process should continue
  # if there are no matches left, index_largest_d will be NA. 
  # If in a previous round a case has matched nc times, it will have been left as NA
  if (length(index_largest_d) == 0) {
    break # Exit the loop if all cases have nc matches or no controls are available
  }
  
  #add the largest_d to the matched_casecontrol list and increment the count for that case
  matched_casescontrols[[index_largest_d]] <- c(matched_casescontrols[[index_largest_d]], 
                                                an(shortest_distances[index_largest_d, 2]))
  shortest_distances[index_largest_d, 3] <- shortest_distances[index_largest_d, 3] + 1
  
  # identify the control used in prev step, and all the cases matched to it. 
  j_value <- shortest_distances[index_largest_d, 2]
  replace_i <- which(shortest_distances[, 2 ] == j_value )
  
  # add NA to the cases and matched controls that are to be replaced. These will remain NA if no further match.
  shortest_distances[replace_i, 1:2] <- NA
  
  
  # if case has been fully matched , do not include it in search for new controls
  # nb if there are no other cases sharing controls and nc has been reached, replace_i will be empty 
  # and the loop should restart
  if (shortest_distances[index_largest_d,3] >= nc) {
    # Remove index_largest_d from the replace_i if nc has been met
    replace_i <- setdiff(replace_i, index_largest_d)
  }  
  
  # find replacement controls j and d, if needed:
  if (length(replace_i) > 0) {
    has_non_na <- numeric()
    j.rmNA=numeric() 
    
    # Identify rows in replace_i with any non-NA value
    has_non_na <- apply(!is.na(wdt[replace_i, , drop = FALSE]), 1, any)
    
    # Check if there are any rows with non-NA values
    if (any(has_non_na)) {
      # Step 2: For rows with any non-NA values, find the minimum and its index
      for (i in replace_i[has_non_na]) {
        r.j_index <- which.min(wdt[i, ])
        r.min_d <- wdt[i,r.j_index]
        j.rmNA = c(j.rmNA, r.j_index)
        shortest_distances[i, 1] <- r.min_d
        shortest_distances[i, 2] <- r.j_index
             }
      # Block out these controls (cols) for further use
      wdt[,unique(j.rmNA)] <- NA
    }
  }
}    # end repeat block 



## identify cases_kept and cases_removed:
controls_kept_index = numeric()
for (i in seq_along(matched_casescontrols)) {
  # Check if the length of the vector is greater than nc
  if (length(matched_casescontrols[[i]]) >= ncmin) {
    # Add the name of the vector to cases_kept
    cases_kept <- c(cases_kept, names(matched_casescontrols)[i])
    controls_kept_index <- an(c(controls_kept_index, matched_casescontrols[[i]]))
  }
}

cases_removed <- setdiff(case_ids, cases_kept) 

controls_kept <- control_ids[controls_kept_index]
controls_removed <- setdiff(control_ids, controls_kept)


## Update dataframes
cases$Group <- ifelse(cases$IID %in% cases_kept, "Case kept", "Case removed")
controls$Group <- ifelse(controls$IID %in% controls_kept, "Control kept", "Control removed")


# create merged list of IIDs to keep. 
outputtable1=cbind(ac(c(cases_kept, controls_kept)))
outputtable2=cbind(ac(c(cases_removed, controls_removed)))



## plotting
mycols_before <- c('Controls' = "blue",
                    'Cases' = "red")
mycolors <- c('Control kept' = "blue",
              'Control removed' ="black",
              'Case kept' = "red",
              'Case removed' = "pink")


theme_set(theme_minimal(base_size = 6) +
            theme(
              plot.title = element_text(size = 8),
              axis.title = element_text(size = 8),
              legend.title = element_text(size = 8),
              legend.text = element_text(size = 8),
              legend.key.size = unit(0.2, "lines"),  # Adjusts the size of the legend keys
              legend.spacing.x = unit(0.2, "cm"),  # Adjusts spacing between legend entries (horizontally)
              legend.spacing.y = unit(0.5, "cm"),  # Adjusts spacing between legend entries (vertically)
              legend.box.margin = margin(1, 1, 1, 1),  # Adjusts the margin around the entire legend box
                                axis.text = element_blank(), # Removes axis text
              axis.ticks = element_blank() # Removes axis ticks
            ))



g1a <- ggplot() +
  geom_point(data = controls, aes(x = PC1, y = PC2, color = "Controls"), size=0.4) +
  geom_point(data = cases, aes(x = PC1, y = PC2, color = "Cases"), size=0.4) +
  scale_color_manual(values = mycols_before) +
  labs(x = "PC1", y = "PC2", color = "Group") +
  coord_fixed()


g1b <-ggplot() +
  geom_point(data=controls, aes(x = PC1, y = PC2, color = Group), size=0.4) +
  geom_point(data=cases, aes(x = PC1, y = PC2, color = Group), size=0.4) +
  scale_color_manual(values = mycolors) +
  labs( x = "PC1", y = "PC2") +
  coord_fixed()

g1c <-ggplot() +
  geom_point(data = subset(controls, Group == "Control kept"), aes(x = PC1, y = PC2, color = "Controls"), size=0.4) +
  geom_point(data = subset(cases, Group == "Case kept"), aes(x = PC1, y = PC2, color = "Cases"), size=0.4) +
  scale_color_manual(values = mycols_before) +
  labs(x = "PC1", y = "PC2", color = "Group") +
  coord_fixed()


g2a <- ggplot() +
  geom_point(data = controls, aes(x = PC1, y = PC3, color = "Controls"), size=0.4) +
  geom_point(data = cases, aes(x = PC1, y = PC3, color = "Cases"), size=0.4) +
  scale_color_manual(values = mycols_before) +
  labs(x = "PC1", y = "PC3", color = "Group") +
  coord_fixed()

g2b <- ggplot() +
  geom_point(data = controls, aes(x = PC1, y = PC3, color = Group), size=0.4) +
  geom_point(data = cases, aes(x = PC1, y = PC3, color = Group), size=0.4) +
  scale_color_manual(values = mycolors) +
  labs( x = "PC1", y = "PC3") +
  coord_fixed()

g2c <-ggplot() +
  geom_point(data = subset(controls, Group == "Control kept"), aes(x = PC1, y = PC3, color = "Controls"), size=0.4) +
  geom_point(data = subset(cases, Group == "Case kept"), aes(x = PC1, y = PC3, color = "Cases"), size=0.4) +
  scale_color_manual(values = mycols_before) +
  labs(x = "PC1", y = "PC3", color = "Group")  +
  coord_fixed()

g3a <- ggplot() +
  geom_point(data = controls, aes(x = PC1, y = PC4, color = "Controls"), size=0.4) +
  geom_point(data = cases, aes(x = PC1, y = PC4, color = "Cases"), size=0.4) +
  scale_color_manual(values = mycols_before) +
  labs(x = "PC1", y = "PC4", color = "Group") +
  coord_fixed()

g3b <- ggplot() +
  geom_point(data = controls, aes(x = PC1, y = PC4, color = Group), size=0.4) +
  geom_point(data = cases, aes(x = PC1, y = PC4, color = Group), size=0.4) +
  scale_color_manual(values = mycolors) +
  labs( x = "PC1", y = "PC4") +
  coord_fixed()

g3c <-ggplot() +
  geom_point(data = subset(controls, Group == "Control kept"), aes(x = PC1, y = PC4, color = "Controls"), size=0.4) +
  geom_point(data = subset(cases, Group == "Case kept"), aes(x = PC1, y = PC4, color = "Cases"), size=0.4) +
  scale_color_manual(values = mycols_before) +
  labs(x = "PC1", y = "PC4", color = "Group") +
  coord_fixed()


g4a <- ggplot() +
  geom_point(data = controls, aes(x = PC1, y = PC5, color = "Controls"), size=0.4) +
  geom_point(data = cases, aes(x = PC1, y = PC5, color = "Cases"), size=0.4) +
  scale_color_manual(values = mycols_before) +
  labs(x = "PC1", y = "PC5", color = "Group") +
  coord_fixed()

g4b <- ggplot() +
  geom_point(data = controls, aes(x = PC1, y = PC5, color = Group), size=0.4) +
  geom_point(data = cases, aes(x = PC1, y = PC5, color = Group), size=0.4) +
  scale_color_manual(values = mycolors) +
  labs( x = "PC1", y = "PC5") +
  coord_fixed()

g4c <-ggplot() +
  geom_point(data = subset(controls, Group == "Control kept"), aes(x = PC1, y = PC5, color = "Controls"), size=0.4) +
  geom_point(data = subset(cases, Group == "Case kept"), aes(x = PC1, y = PC5, color = "Cases"), size=0.4) +
  scale_color_manual(values = mycols_before) +
  labs(x = "PC1", y = "PC5", color = "Group") +
  coord_fixed()

#### RESULTS #####
current_date <- Sys.Date()
filename_append = paste0("_threshold", threshcutoff, "_maxmatches", nc, "_", current_date)
fileoutkept=paste0(fileoutkept, filename_append, ".txt")
exp_out = paste0("samplesbypop", filename_append,".txt")
fileoutremoved=paste0(fileoutremoved, filename_append, ".txt")
GRAPHPRE=paste0(graphstem, ".allsamples.", filename_append, ".jpg")
GRAPHSPLIT=paste0(graphstem, ".splitsamples.", filename_append, ".jpg")
GRAPHKEPT=paste0(graphstem, ".onlykept.", filename_append, ".jpg")

nctrlskept = length(controls_kept)
ncaseskept = length(cases_kept)
print(text_choice)
text_results=(paste0("# Results: ", nrow(controls), " controls and ",nrow(cases), " cases processed. ", nctrlskept , " controls matched to ", ncaseskept, " cases."))
print(text_results)
print(paste0("Cases removed initially for <", ncmin, " possible matches: ", initial_cases_removed))
#### --- saving out
#
before= arrangeGrob(g1a, g2a, g3a, g4a, nrow = 2)
after = arrangeGrob(g1b, g2b, g3b, g4b, nrow = 2 )
after2 = arrangeGrob(g1c, g2c, g3c, g4c, nrow = 2 )

# writing out files
outdir=paste0(OUTDIRP,"/", FILEOUTSTEM, "_",filename_append)
dir.create(outdir, recursive = F)
setwd(outdir)
ggsave(GRAPHPRE, plot = before,device = "jpeg", path = outdir,  width = 20, height = 16, dpi = 300)
ggsave(GRAPHSPLIT,
        grid.arrange(after,top= textGrob(paste0(text_choice, " ", text_results),
                                        gp = gpar(fontsize = 10, font =2))),
        device = "jpeg", path = outdir, width = 20, height = 16, dpi = 300)
ggsave(GRAPHKEPT, plot = after2 ,device = "jpeg", path = outdir,  width = 20, height = 16, dpi = 300)

write(paste0("# Ancestry Matching: ", text_choice, "\n", text_results), file = fileoutkept )
write.table(outputtable1,file=fileoutkept, quote=F,row.names=F,col.names=F,sep="\t", append = TRUE)
write.table(outputtable2,file=fileoutremoved,quote=F,row.names=F,col.names=F,sep="\t", append = TRUE)


