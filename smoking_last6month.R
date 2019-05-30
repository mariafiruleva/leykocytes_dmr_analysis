library(dplyr)
library(ggplot2)

#load data
read_data <- function(path, desired_pattern){
  xls_files <- dir(path, pattern = desired_pattern)
  files_path <- paste0(path, xls_files)
  all_data <- do.call(rbind, lapply(files_path, readxl::read_excel))
  return(all_data)
}
smoking_51 <- read_data("~/Bioinformatics/Smoking/Spring2019/", "smoking_51_subjects.xls")

#filter data
smoking_51$smoking_last6months <- smoking_51$SMOKEDAY
smoking_51$smoking_last6months[is.na(smoking_51$smoking_last6months)] <- 0
smoking_51$smoking_last6months <- factor(smoking_51$smoking_last6months, levels = c(0:5), labels = c("Don't smoke", "1 cigarette", "2-5 cigarette", 
                                                           "6-10 cigarette", "11-20 cigarette", ">20 cigarette"))
#build a table
smoking_last6months <- as.data.frame(table(smoking_51$smoking_last6months))
percentages_smoking_last6months <- as.data.frame(prop.table(table(smoking_51$smoking_last6months)))
colnames(smoking_last6months) <- c("Cigarette", "Frequency")
colnames(percentages_smoking_last6months) <- c("Cigarette", "Percentage")
percentages_smoking_last6months$Percentage <- unlist(lapply(percentages_smoking_last6months$Percentage, function(x) x * 100))

#visualize distribution
distribution_plot <- ggplot(data = smoking_last6months, aes(x = Cigarette, y = Frequency, fill = Cigarette)) + 
  geom_bar(stat="identity", position = "dodge", colour = "black")+
  theme_bw()+
  labs(x = "Cigarette number", y = "Participant number", title = "Smoking last 6 months per day")

percentages_distribution_plot <- ggplot(data = percentages_smoking_last6months, aes(x = Cigarette, y = Percentage, fill = Cigarette)) + 
  geom_bar(stat="identity", position = "dodge", colour = "black")+
  theme_bw()+
  labs(x = "Cigarette number", y = "Percentage of participant number", title = "Smoking last 6 months per day")

#save results
write.csv(smoking_51, "smoking_last6months.csv")
png("distribution_plot.png")
plot(distribution_plot)
dev.off()

png("percentage_distribution_plot.png")
plot(percentages_distribution_plot)
dev.off()
