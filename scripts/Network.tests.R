install.packages("igraph")
install.packages("network") 
install.packages("sna") 
install.packages("ndtv")
install.packages("ergm")
install.packages("plyr")
library(igraph)
library(network)
library(sna)
library(ergm)
library(plyr)


 #### Load the data ####


Rel.Start <- read.csv('./input.files/RelationshipStart.csv', header=T, as.is=T) 
#Rel.End <- read.csv("RelationshipEnd.csv", header=T, as.is=T) 

# Rel.Start$Rel_type..0...TRANSITORY..1...INFORMAL..2...MARITAL..3...COMMERCIAL
# Relationship types:
# 0 = TRANSITORY
# 1 = INFORMAL
# 2 = MARITAL
# 3 = COMMERCIAL

Rel.Start$Rel_type <- Rel.Start$Rel_type..0...TRANSITORY..1...INFORMAL..2...MARITAL..3...COMMERCIAL.
Rel.Start$Rel_start_time_YEARS <- Rel.Start$Rel_start_time/365

Rel.Start.list <- subset(Rel.Start, select = c("A_ID", "B_ID", "Rel_type", "Rel_start_time_YEARS", "A_risk"))
Rel.Start.list.1year <- Rel.Start.list[Rel.Start.list$Rel_start_time_YEARS > 34.0 & 
                                         Rel.Start.list$Rel_start_time_YEARS < 35.0, ]
#write.csv(Rel.Start.list.1year, file = "Eswatini.one.year.csv")

aaa <- duplicated(Rel.Start.list.1year[, 1:2])
Rel.Start.list.1year.test <- Rel.Start.list.1year[!aaa, ]
#write.csv(Rel.Start.list.1year.test, file = "Eswatini.one.year.deduplicated.csv")


#### Network visualization ####

  
# Relationship types:
# 0 = TRANSITORY
# 1 = INFORMAL
# 2 = MARITAL
# 3 = COMMERCIAL

net0 <- as.network(Rel.Start.list.1year.test[Rel.Start.list.1year.test$Rel_type == 0,], directed = FALSE)
net1 <- as.network(Rel.Start.list.1year.test[Rel.Start.list.1year.test$Rel_type == 1,], directed = FALSE)
net2 <- as.network(Rel.Start.list.1year.test[Rel.Start.list.1year.test$Rel_type == 2,], directed = FALSE)
net3 <- as.network(Rel.Start.list.1year.test[Rel.Start.list.1year.test$Rel_type == 3,], directed = FALSE)

# Transitory
comp.dist.net0 <- sna:::component.dist(net0)
length(comp.dist.net0$cdist)

# Informal
comp.dist.net1 <- sna:::component.dist(net1)
length(comp.dist.net1$cdist)

# Marital
comp.dist.net2 <- sna:::component.dist(net2)
length(comp.dist.net2$cdist)

# Commercial
comp.dist.net3 <- sna:::component.dist(net3)
length(comp.dist.net3$cdist)

par(mfrow=c(2,2), mar=c(2,2,2,2))
plot.network(net0, 
             vertex.cex = 0.07,
             displaylabels = F,
             jitter = TRUE,
             edge.len = 0.1, 
             edge.curve = 0.001, 
             object.scale = 0.05, 
             edge.lwd=0.2,
             main = paste(length(comp.dist.net0$cdist)," transitory partnerships"))
plot.network(net1, 
             vertex.cex = 0.07,
             displaylabels = F,
             jitter = TRUE,
             edge.len = 0.1, 
             edge.curve = 0.001, 
             object.scale = 0.05, 
             edge.lwd=0.2,
             main = paste(length(comp.dist.net1$cdist)," informal partnerships"))
plot.network(net2, 
             vertex.cex = 0.07,
             displaylabels = F,
             jitter = TRUE,
             edge.len = 0.1, 
             edge.curve = 0.001, 
             object.scale = 0.05, 
             edge.lwd=0.2,
             main = paste(length(comp.dist.net2$cdist)," marital partnerships"))
plot.network(net3, 
             vertex.cex = 0.07,
             displaylabels = F,
             jitter = TRUE,
             edge.len = 0.1, 
             edge.curve = 0.001, 
             object.scale = 0.1, 
             edge.lwd=0.2,
             main = paste(length(comp.dist.net3$cdist)," commercial partnerships"))


#### Partners per year ####

  
# from "Rel.Start.list.1year.test" from above

bbb = plyr::count(Rel.Start.list.1year.test$A_ID)
# counts the number of unique individual IDs (A_ID), which is the number
# of partners each unique individual has per year (well, in the year chosen on line 25)

Eswatini.counts <- plyr::count(bbb$freq)
# counts the number of individuals with the same number of partners per year

plot(Eswatini.counts$x, log10(Eswatini.counts$freq/sum(Eswatini.counts$freq)), 
     col = "red",
     ylab = "log10 frequency",
     xlab = "Partners per year",
     pch = 16)

### Atlanta

Atlanta.counts <- read.csv("Atlanta/edge_counts.csv", header=T) 

# Log10 frequency

plot(Atlanta.counts$count, log10(Atlanta.counts$sim1_1/sum(Atlanta.counts$sim1_1)), 
     col = "red",
     ylab = "log10 frequency",
     xlab = "Partners per year",
     pch = 16,
     ylim = c(-5, 0),
     main = "Partners per year")
text(20, log10(0.5), paste(sum(Atlanta.counts$sim1_1), " total individuals with relationships, Atlanta"))

points(Eswatini.counts$x, log10(Eswatini.counts$freq/sum(Eswatini.counts$freq)), 
       col = "blue",
       pch = 16)
text(20, log10(0.0001), paste(sum(Eswatini.counts$freq), " total individuals with relationships, Eswatini"))

legend(30, -0.5, 
       c("Eswatini", "Atlanta"),
       pch = 16,
       col = c("blue", "red"))

# Un-transformed frequency
par(mfrow=c(1,1))
plot(Atlanta.counts$count, Atlanta.counts$sim1_1/sum(Atlanta.counts$sim1_1), 
     col = "red",
     ylab = "Frequency",
     xlab = "Partners per year",
     pch = 16,
     ylim = c(0, 1),
     main = "Partners per year")
  
points(Eswatini.counts$x, Eswatini.counts$freq/sum(Eswatini.counts$freq),
       col = "blue",
       pch = 16)

legend(30, 1, 
       c("Eswatini", "Atlanta"),
       pch = 16,
       col = c("blue", "red"))

### Kind of adjusted by rounding and using 50000 total N for both

Eswatini.adjusted.counts <- Eswatini.counts$freq/sum(Eswatini.counts$freq) * 50000
Eswatini.adjusted.counts.round <- round(Eswatini.adjusted.counts)

Atlanta.adjusted.counts <- Atlanta.counts$sim1_1/sum(Atlanta.counts$sim1_1) * 50000
Atlanta.adjusted.counts.round <- round(Atlanta.adjusted.counts)

plot(Atlanta.counts$count, log10(Atlanta.adjusted.counts.round/sum(Atlanta.adjusted.counts.round)), 
     col = "red",
     ylab = "log10 frequency",
     xlab = "Partners per year",
     pch = 16,
     ylim = c(-5, 0),
     main = "Partners per year")
text(20, log10(0.5), paste(sum(Atlanta.adjusted.counts.round), " total individuals with relationships, Atlanta"))

points(Eswatini.counts$x, log10(Eswatini.adjusted.counts.round/sum(Eswatini.adjusted.counts.round)), 
       col = "blue",
       pch = 16)
text(20, log10(0.0001), paste(sum(Eswatini.adjusted.counts.round), " total individuals with relationships, Eswatini"))

legend(30, -0.5, 
       c("Eswatini", "Atlanta"),
       pch = 16,
       col = c("blue", "red"))


----------------------------------------------
#### Fit distribution to Eswatini counts ####
----------------------------------------------
  
library(fitdistrplus)

fit <- fitdist(bbb$freq, "nbinom")

#Fitting of the distribution ' nbinom ' by maximum likelihood 
#Parameters:
#     estimate  Std. Error
#size 0.8906601 0.04927525
#mu   2.8581882 0.12207539

hist(bbb$freq, breaks = 100, 
     col = "red",
     main = "Eswatini partners per year",
     xlab = "Partners per year",
     ylab = "Count")

fitD <- rnbinom(10000, size = 0.89, mu = 2.858)

fitD.counts <- plyr::count(fitD)

hist(fitD+1, breaks = 100)


#### Partner distributions by relational type ####


# https://docs.idmod.org/projects/emod-hiv/en/latest/software-report-relationship-start.html
# <A or B>_extra_relational_bitmask
# 8 = 1000
# 9 = 1001
# 11 = 1010
# 15 = 1111

Rel.Start <- read.csv("RelationshipStart.csv", header=T, as.is=T) 

Rel.Start$Rel_type <- Rel.Start$Rel_type..0...TRANSITORY..1...INFORMAL..2...MARITAL..3...COMMERCIAL.
Rel.Start$Rel_start_time_YEARS <- Rel.Start$Rel_start_time/365

Rel.Start.list <- subset(Rel.Start, 
                         select = c("A_ID", "B_ID", "Rel_type", "A_extra_relational_bitmask", "Rel_start_time_YEARS"))
Rel.Start.list.1year <- Rel.Start.list[Rel.Start.list$Rel_start_time_YEARS > 34.0 & 
                                         Rel.Start.list$Rel_start_time_YEARS < 35.0, ]
#write.csv(Rel.Start.list.1year, file = "Eswatini.one.year.csv")

aaa <- duplicated(Rel.Start.list.1year[, 1:2])
Rel.Start.list.1year.test <- Rel.Start.list.1year[!aaa, ]

Rel.Start.bitmask.8 <- Rel.Start.list.1year.test[Rel.Start.list.1year.test$A_extra_relational_bitmask == 8,]
Rel.Start.bitmask.9 <- Rel.Start.list.1year.test[Rel.Start.list.1year.test$A_extra_relational_bitmask == 9,]
Rel.Start.bitmask.11 <- Rel.Start.list.1year.test[Rel.Start.list.1year.test$A_extra_relational_bitmask == 11,]
Rel.Start.bitmask.15 <- Rel.Start.list.1year.test[Rel.Start.list.1year.test$A_extra_relational_bitmask == 15,]

bbb.8 = plyr::count(Rel.Start.bitmask.8$A_ID)
bbb.9 = plyr::count(Rel.Start.bitmask.9$A_ID)
bbb.11 = plyr::count(Rel.Start.bitmask.11$A_ID)
bbb.15 = plyr::count(Rel.Start.bitmask.15$A_ID)

# 8 = 1000
# 9 = 1001
# 11 = 1010
# 15 = 1111

par(mfrow=c(2,2), mar=c(4,4,1,1), oma = c(0, 0, 2, 0))
hist(bbb.8$freq, 
     #breaks = 100, 
     col = "red",
     main = "",
     xlab = "Partners per year",
     ylab = "Count")
text(1.4, 200, "8 / 1000")
hist(bbb.9$freq, 
     #breaks = 100, 
     col = "red",
     main = "",
     xlab = "Partners per year",
     ylab = "Count")
text(2, 100, "9 / 1001")
hist(bbb.11$freq, 
     breaks = 100, 
     col = "red",
     main = "",
     xlab = "Partners per year",
     ylab = "Count")
text(3, 40, "11 / 1010")
hist(bbb.15$freq, 
     breaks = 100, 
     col = "red",
     main = "",
     xlab = "Partners per year",
     ylab = "Count")
text(15, 20, "15 / 1111")
mtext("Partners per year by potential relationship types", outer = TRUE, cex = 1)

# counts the number of unique individual IDs (A_ID), which is the number
# of partners each unique individual has per year (well, in the year chosen on line 25)

Eswatini.counts.8 <- plyr::count(bbb.8$freq)
Eswatini.counts.9 <- plyr::count(bbb.9$freq)
Eswatini.counts.11 <- plyr::count(bbb.11$freq)
Eswatini.counts.15 <- plyr::count(bbb.15$freq)
# counts the number of individuals with the same number of partners per year

par(mfrow=c(2,2), mar=c(2,2,2,2))
plot(Eswatini.counts.8$x, log10(Eswatini.counts.8$freq/sum(Eswatini.counts.8$freq)), 
     col = "red",
     ylab = "log10 frequency",
     xlab = "Partners per year",
     pch = 16)


#### Partner distributions by risk ####


Rel.Start.list.1year <- Rel.Start.list[Rel.Start.list$Rel_start_time_YEARS > 34.0 & 
                                         Rel.Start.list$Rel_start_time_YEARS < 35.0, ]
#write.csv(Rel.Start.list.1year, file = "Eswatini.one.year.csv")

aaa <- duplicated(Rel.Start.list.1year[, 1:2])
Rel.Start.list.1year.test <- Rel.Start.list.1year[!aaa, ]

Rel.Start.low <- Rel.Start.list.1year.test[Rel.Start.list.1year.test$A_risk == "Risk-LOW",]
Rel.Start.medium <- Rel.Start.list.1year.test[Rel.Start.list.1year.test$A_risk == "Risk-MEDIUM",]
Rel.Start.high <- Rel.Start.list.1year.test[Rel.Start.list.1year.test$A_risk == "Risk-HIGH",]

bbb.low = plyr::count(Rel.Start.low$A_ID)
bbb.medium = plyr::count(Rel.Start.medium$A_ID)
bbb.high = plyr::count(Rel.Start.high$A_ID)

par(mfrow=c(2,2), mar=c(4,4,1,1), oma = c(0, 0, 2, 0))
hist(bbb.low$freq, 
     #breaks = 100, 
     col = "red",
     main = "",
     xlab = "Partners per year",
     ylab = "Count",
     xlim = c(1, 35))
text(10, 200, "Low risk")
hist(bbb.medium$freq, 
     #breaks = 100, 
     col = "red",
     main = "",
     xlab = "Partners per year",
     ylab = "Count",
     xlim = c(1, 35))
text(10, 100, "Medium risk")
hist(bbb.high$freq, 
     breaks = 100, 
     col = "red",
     main = "",
     xlab = "Partners per year",
     ylab = "Count",
     xlim = c(1, 35))
text(10, 6, "High risk")
mtext("Partners per year by individual risk", outer = TRUE, cex = 1.5)




#### Contact tracing ####

# Data manipulation to start

Rel.Start$Rel_type <- Rel.Start$Rel_type..0...TRANSITORY..1...INFORMAL..2...MARITAL..3...COMMERCIAL.
Rel.Start$Rel_start_time_YEARS <- Rel.Start$Rel_start_time/365

# A_IDs are males
# B_IDs are females

# Some possible behavioral, clinical, and demographic risk attributes of interest each indivdual
# "A_ or B_risk"
# "A_ or B_HIV_disease_stage"
# "A_ or B_is_infected"
# "A_ or B_is_superspreader"
# "A_ or B_HIV_Tested_Positive"

# test:  list all contacts for one individual (this is for the entire modeled lifespan)

#Rel.Start.list <- Rel.Start.list[order(Rel.Start.list$A_ID),]
contacts.ID.2 <- Rel.Start$B_ID[Rel.Start$A_ID == 2]

#> contacts.ID.2
#[1] 2980 3388 3181 4375 1485 3472 5819 3393 4249 3537



# Now:

# Pick a timespan of the simulated epidemic (e.g. 1 or 2 years?)
# (this N-year timespan is a proxy for a contact tracing and genetic sampling campaign)

Rel.Start.1years <- Rel.Start[Rel.Start$Rel_start_time_YEARS > 50.0 & 
                                Rel.Start$Rel_start_time_YEARS < 51.0, ]


A_ID_HIVinfected <- Rel.Start.1years[Rel.Start.1years$A_HIV_Tested_Positive == 1, ]
# str(A_ID_HIVinfected)
# 'data.frame':	315 obs. of  53 variables

B_ID_HIVinfected <- Rel.Start.1years[Rel.Start.1years$B_HIV_Tested_Positive == 1, ]
# str(B_ID_HIVinfected)
# 'data.frame':	2185 obs. of  53 variables



length(unique(Rel.Start.1years$A_ID))
# Number of total possible index individuals 

# Select only the HIV-infected and diagnosed individuals from those A_ID index possibilities
index.HIVinfected <- Rel.Start.1years[Rel.Start.1years$A_HIV_Tested_Positive == 1, ]

length(unique(index.HIVinfected$Rel_ID))
# 315 unique relationships during that time that involved at least one HIV+ individual
length(unique(index.HIVinfected$A_ID))
# 62 unique A_ID individuals
length(unique(index.HIVinfected$B_ID))
# 119 unique B_ID

# Select a set of index non-HIV-infected individuals
index.HIVuninfected <- Rel.Start.1years[Rel.Start.1years$A_HIV_Tested_Positive == 0, ]

length(unique(index.HIVuninfected$A_ID))
# 804 unique A_ID individuals
length(unique(index.HIVuninfected$B_ID))
# 798 unique B_ID 


# For each index in "index.HIVinfected", identify all contacts (all B_IDs)
# that are prior to a certain date (within the timespan)

partner.HIVinfected <- index.HIVinfected$B_ID
length(partner.HIVinfected)

# for all rows in which A_is_infected == 1, identify all B_IDs
# 

Rel.Start$B_ID[Rel.Start$A_ID == 2]

