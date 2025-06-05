library(ggplot2)

setwd("/home/ubuntu/")
female_rate = read.csv("~/Aging_thyroid/data/SEER_thyroidCancer_rate/SEER_thyroidCancer_rate_bySubtype_female.csv")[-1,]
male_rate = read.csv("~/Aging_thyroid/data/SEER_thyroidCancer_rate/SEER_thyroidCancer_rate_bySubtype_male.csv")[-1,]

rate = c(female_rate$Follicular.Subtype.Carcinomas,female_rate$Papillary.Subtype.Carcinomas,
  male_rate$Follicular.Subtype.Carcinomas,male_rate$Papillary.Subtype.Carcinomas)
fulldata = data.frame(rate)
fulldata$age = c(female_rate$X, female_rate$X, male_rate$X, male_rate$X)
fulldata$disease = c(rep("FTC",length(female_rate$Follicular.Subtype.Carcinomas)),
                     rep("PTC",length(female_rate$Papillary.Subtype.Carcinomas)),
                     rep("FTC",length(male_rate$Follicular.Subtype.Carcinomas)),
                     rep("PTC",length(male_rate$Papillary.Subtype.Carcinomas)))
fulldata$disease = factor(fulldata$disease)
fulldata$sex = c(rep("Female",length(female_rate$Follicular.Subtype.Carcinomas)),
                     rep("Female",length(female_rate$Papillary.Subtype.Carcinomas)),
                     rep("Male",length(male_rate$Follicular.Subtype.Carcinomas)),
                     rep("Male",length(male_rate$Papillary.Subtype.Carcinomas)))
fulldata$sex = factor(fulldata$sex, levels = c("Female", "Male"))

# Replace missing data designated by "^" by 0
fulldata$rate = as.numeric(fulldata$rate)
fulldata[fulldata == 0] = NA

fulldata$age[which(fulldata$age == "1-4")] = "01-04"
fulldata$age[which(fulldata$age == "5-9")] = "05-09"

# Plot rates per 100,000
maintitle = "Incidence rates of thyroid cancers by age at diagnosis and sex"
ggplot(fulldata, aes(x = age, y = rate, color = sex, group = interaction(disease, sex))) + geom_point() + geom_path(aes(linetype=disease)) +
  theme_bw() + ggtitle(maintitle) + xlab("Age") + ylab("Rate per 100,000")


# Plot each disease separately
fulldata_ptc = fulldata[which(fulldata$disease == "PTC"),]
maintitle = "Incidence rates of PTC by age at diagnosis and sex"
ggplot(fulldata_ptc, aes(x = age, y = rate, color = sex, group = sex)) + geom_point() + geom_path() +
  theme_bw() + ggtitle(maintitle) + xlab("Age") + ylab("Rate per 100,000")

fulldata_ftc = fulldata[which(fulldata$disease == "FTC"),]
maintitle = "Incidence rates of FTC by age at diagnosis and sex"
ggplot(fulldata_ftc, aes(x = age, y = rate, color = sex, group = sex)) + geom_point() + geom_path() +
  theme_bw() + ggtitle(maintitle) + xlab("Age") + ylab("Rate per 100,000")
