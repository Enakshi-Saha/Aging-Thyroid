library(RCurl)
url1 <- "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_TST.xpt"
destination1 <- "/home/esaha/NHANES_hormone/P_TST.xpt"

url2 <- "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_DEMO.xpt"
destination2 <- "/home/esaha/NHANES_hormone/P_DEMO.xpt"

# Download the file
download.file(url1, destfile = destination1, mode = "wb")
download.file(url2, destfile = destination2, mode = "wb")

library(haven)
library(ggplot2)
library(gghalves)

hormone = read_xpt("/home/esaha/NHANES_hormone/P_TST.xpt")
demography = read_xpt("/home/esaha/NHANES_hormone/P_DEMO.xpt")
demography = demography[match(hormone$SEQN,demography$SEQN),]

age = demography$RIDAGEYR
gender = demography$RIAGENDR
gender[which(gender == 1)] = "Male"
gender[which(gender == 2)] = "Female"
gender = factor(gender)

age_at_thyroidCancer = medical$MCQ240BB

# Hormone names
progesterone = hormone$LBXPG4 # unit = ng/dL
estradiol = hormone$LBXEST # unit = pg/mL
SHBG = hormone$LBXSHBG #sex_hormone_binding_globulin; unit = nmol/L

# Plot hormone level by age for males and females
fulldata = data.frame(cbind(age, progessterone, estradiol, SHBG, age_at_thyroidCancer))
fulldata$thyroid_cancer = rep("No", nrow(fulldata))
fulldata$thyroid_cancer[which(!is.na(age_at_thyroidCancer))] = "Yes"
fulldata_female = fulldata[which(gender == "Female"),]
fulldata_male = fulldata[which(gender == "Male"),]
head(fulldata_female)
head(fulldata_male)

maintitle = "NHANES female: log progesterone by age"
ggplot(fulldata_female, aes(x = age, y = log(1+testosterone))) + geom_point(colour = "grey") +
  theme_bw() + ggtitle(maintitle) + xlab("Age") + ylab("log ng/dL") + geom_smooth(color = "blue")

maintitle = "NHANES female: log estradiol by age"
ggplot(fulldata_female, aes(x = age, y = log(1+estradiol))) + geom_point(colour = "grey") +
  theme_bw() + ggtitle(maintitle) + xlab("Age") + ylab("log pg/mL") + geom_smooth(color = "red")

maintitle = "NHANES male: log progesterone by age"
ggplot(fulldata_male, aes(x = age, y = log(1+testosterone))) + geom_point(colour = "grey") +
  theme_bw() + ggtitle(maintitle) + xlab("Age") + ylab("log ng/dL") + geom_smooth(color = "blue")

maintitle = "NHANES male: log estradiol by age"
ggplot(fulldata_male, aes(x = age, y = log(1+estradiol))) + geom_point(colour = "grey") +
  theme_bw() + ggtitle(maintitle) + xlab("Age") + ylab("log pg/mL") + geom_smooth(color = "red")

# Compare thyroid samples with and without thyroid cancer (only 2 male samples, so we work with female)
thyroid_cancer_samples_F = data.frame(fulldata_female[!is.na(fulldata$age_at_thyroidCancer),])
normal_samples_F = data.frame(fulldata_female[is.na(fulldata$age_at_thyroidCancer),])
estro = c(thyroid_cancer_samples_F$estradiol, normal_samples_F$estradiol)
newdata = data.frame(estro)
newdata$thyroid_cancer = c(rep("Yes", nrow(thyroid_cancer_samples_F)), rep("No", nrow(normal_samples_F)))

maintitle = "NHANES female: log estrogen in normal vs thyroid cancer samples"
ggplot(newdata, aes(x = thyroid_cancer, y = log(1+estro))) + geom_half_violin(side = "r") + geom_half_boxplot(side = "l") +
  theme_bw() + ggtitle(maintitle) + xlab("ever had thyroid cancer") + ylab("log ng/dL")

# Test if estradiol level has any difference between disease and normal
wilcox.test(thyroid_cancer_samples_F$estradiol, normal_samples_F$estradiol)
