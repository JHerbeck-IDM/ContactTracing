#setwd("~/Dropbox/EMOD.infections")

# Goal:  get a distribution of the number of infection-causing contacts for a population. 
# So the x-axis will be "infection causing contacts" and the y-axis is the count 
# (or frequency if we prefer that)

# The counts will be, for each individual that has a "NewInfectionEvent", 
# how many "RedundantInfectionEvent" they have + 1 (for the initial infection).

df <- read.csv("ReportEventRecorder (2).csv")

length(unique(df$Individual_ID))
#53,330

length(which(df$Event_Name == "NewInfectionEvent"))
# 25,102

# Number of individuals who were not infected (did not experience an 
# infection-causing exposures):

length(unique(df$Individual_ID)) - length(which(df$Event_Name == "NewInfectionEvent"))
# 28,228 with 0 infection-causing exposures


length(which(df$Event_Name == "RedundantTransmissionEvent"))
# 297,855 additional infection-causing exposures (spread out over 25,102 individuals)

df.Redundant <- df[df$Event_Name == "RedundantTransmissionEvent", ]
df.STIDebut <- df[df$Event_Name == "STIDebut", ]

aaa <- aggregate(df.Redundant$Event_Name == "NewInfectionEvent", by=list(df.Redundant$Individual_ID), FUN=length)
colnames(aaa)[colnames(aaa) == 'Group.1'] <- 'Individual_ID'
colnames(aaa)[colnames(aaa) == 'x'] <- 'Redundant_infections'
aaa$Infecting_exposures <- aaa$Redundant_infections + 1

# Number of individuals with just one infection-causing exposure = 
# Total infected individuals - Infected individuals that have a redundant infecting exposure
# 25,102 - length(aaa$Group.1)
# 25,102 - 20,995 = 4,107

Individual_ID <- seq(20996, (20996+4106), by=1)
Infecting_exposures <- replicate(4107, 1)
df.One.Infection <- data.frame(Individual_ID, Infecting_exposures)

Individual_ID <- seq(25103, (25103+28227), by=1)
Infecting_exposures <- replicate(28228, 0)
df.Zero.Infection <- data.frame(Individual_ID, Infecting_exposures)

aaa$Redundant_infections <- NULL
aaa <- rbind(aaa, df.One.Infection, df.Zero.Infection)

hist.out <- hist(aaa$Infecting_exposures, 
     breaks = 200, 
     main = "Total infecting exposures",
     xlab = "Number of infecting exposures",
     ylab = "Number of individuals",
     freq = FALSE,
     xlim = c(0, 25))

