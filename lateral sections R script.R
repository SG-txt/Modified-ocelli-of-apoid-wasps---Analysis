####set-up####

  library(ggplot2)
library ("plyr")
library("dplyr")
library("tidyr")
  
setwd("C:/Users/Stefan/Desktop/Dr/Analysen/Lateral sections/Rscript & analyse")
wd <- getwd()

if (dir.exists(paste0(wd,"/save/"))){
  print("folder already exists")
} else {
  dir.create(paste0(wd,"/save/"))}
####stats####
####!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NUMBER OF MISSING FILES!!!!!!!!!!!!!!!!!!!!!!!!!!!!####
m <- 2 #2 11/12_vt files of Stictiella pulchella, as there is no visible vt retina
m <- m/2
####script####
{

 
  wd <- getwd()
  
  dist_master <- data.frame(
    genus = character(),
    species = character(),
    sex = character(),
    ID = character(),
    type = character(),
    mean = numeric(),
    sd = numeric())  

length_master <- data.frame(
  genus = character(),
  species = character(),
  sex = character(),
  ID = character(),
  type = character(),
  length = numeric(),
  region = character())

Rblength_mean_sd <- data.frame(
  genus = character(),
  species = character(),
  sex = character(),
  ID = character(),
  type = character(),
  mean = numeric(),
  sd = numeric(),
  min = numeric(),
  max = numeric(),
  n = numeric(),
  region = character())


#Nr of .txt files in wd to use as limit for loop
files <- list.files(pattern = "\\.txt$")
txt_file_count <- length(files)


i<-1
while (i <= txt_file_count/2+m) {
  #data frame for the results 
  if (file.exists(paste0(i,"_do.txt"))){
    
  
   
  temp_table <- data.frame(
    genus = character(),
    species = character(),
    sex = character(),
    ID = character(),
    type = character(),
    length = character(),
    region = character())
  
  temp_table2 <- data.frame(
    genus = character(),
    species = character(),
    sex = character(),
    ID = character(),
    type = character(),
    dist = character())
  
  #process dorsal retina measurements
  RAW <- read.table(paste0(i,"_do.txt"), header = TRUE, skip = 1)
  RAW$X <- as.numeric(RAW$X)
  RAW$Y <- as.numeric(RAW$Y)
  RAW <- RAW[, !(names(RAW) %in% c("Area", "Mean","Min","Max","time"))]
  region <- "do"
  table <- read.table(paste0(i,"_do.txt"), header = FALSE, nrows = 1, sep = ",")
  genus <- table [,1]
  species <- table [,2]
  ID <- table [,3]
  type <- table [,4]
  sblength <- table [,5]
  sex <- table [,6]
  mod <- table[,7]
  
  #Factor to calculate length in micrometer from the measurement on the scalebar 
  # point 1 and point 2 are the respective points of the scalebar in every Dataset 
  sbfactor <- (sqrt((RAW[2,1]-RAW[1,1])^2+(RAW[2,2]-RAW[1,2])^2))
  
  #calculate the distance between lens and fovea 
  j <- 3
  while (j<=8) {
    dist  <- (sqrt((RAW[j+1,1]-RAW[j,1])^2+(RAW[j+1,2]-RAW[j,2])^2))*sblength/sbfactor
    
    temp_row <- c(genus,species,sex,ID,type,dist)
    
    temp_table2 <- rbind(temp_table2,temp_row)
    
    j <- j+2
  }
  colnames(temp_table2) <- c("genus","species","sex","ID","type","dist")
  temp_table2$dist <- as.numeric(temp_table2$dist)
  mean_dist <- mean(temp_table2$dist)
  sd_dist<- sd(temp_table2$dist)
  temp_row <- c(genus,species,sex,ID,type,mean_dist,sd_dist)
  dist_master <- rbind(dist_master,temp_row)
  
  #calculate the actual length of each Rhabdom 
  j <- 9
  while (j<=(length(RAW$X))) {
  length  <- (sqrt((RAW[j+1,1]-RAW[j,1])^2+(RAW[j+1,2]-RAW[j,2])^2))*sblength/sbfactor
  
  temp_row <- c(genus,species,sex,ID,type,length,region)
  
  temp_table <- rbind(temp_table,temp_row)
  length_master <- rbind(length_master,temp_row)
  j <- j+2
  }
  #calculate and save results
  colnames(temp_table) <- c("genus","species","sex","ID","type","length","region")
  temp_table$length <- as.numeric(temp_table$length)
  #save RAW length in .txt file 
    temp_name <- paste0("Rblength_lat_",genus,"_",ID,"_",sex,"_",type,"_",region) 
  assign(temp_name,temp_table)
  
  filename <- paste0(wd,"/save/",temp_name, ".txt")
  write.table(get(temp_name), file = filename, sep = "\t", quote = FALSE, row.names = FALSE)
  
  #take means from dorsal and ventral retina + sd and add to sepperate table 
  mean <- mean(temp_table$length)
  sd <- sd(temp_table$length)
  min <- min(temp_table$length)
  max <- max(temp_table$length)
  n <- length(temp_table$length)
  temp_row <- c(genus,species,sex,ID,type,mean,sd, min, max, n, region)
  Rblength_mean_sd <- rbind(Rblength_mean_sd,temp_row)
  #delete temp table
  rm(temp_table)
  
  } else {
    print(paste0(i,"_do.txt is missing"))}
  
  if (file.exists(paste0(i,"_vt.txt"))){
  #process ventral retina measurements
  RAW <- read.table(paste0(i,"_vt.txt"), header = TRUE, skip = 1)
  RAW$X <- as.numeric(RAW$X)
  RAW$Y <-  as.numeric(RAW$Y)
  RAW <- RAW[, !(names(RAW) %in% c("Area", "Mean","Min","Max","time"))]
  region <- "vt"
  table <- read.table(paste0(i,"_do.txt"), header = FALSE, nrows = 1, sep = ",")
  genus <- table [,1]
  species <- table [,2]
  ID <- table [,3]
  type <- table [,4]
  sblength <- table [,5]
  sex <- table [,6]
  mod <- table[,7]
  
  
  temp_table <- data.frame(
    genus = character(),
    species = character(),
    sex = character(),
    ID = character(),
    type = character(),
    length = character(),
    region = character())
  
  #Factor to calculate length in micrometer from the measurement on the scalebar 
  # point 1 and point 2 are the respective points of the scalebar in every Dataset 
  sbfactor <- (sqrt((RAW[2,1]-RAW[1,1])^2+(RAW[2,2]-RAW[1,2])^2))
  
  
  
  #calculate the actual length of each Rhabdom 
  j <- 3
  while (j<=(length(RAW$X))) {
    length  <- (sqrt((RAW[j+1,1]-RAW[j,1])^2+(RAW[j+1,2]-RAW[j,2])^2))*sblength/sbfactor
    
    temp_row <- c(genus,species,sex,ID,type,length,region)
    
    temp_table <- rbind(temp_table,temp_row)
    length_master <- rbind(length_master,temp_row)
    j <- j+2
  }
  
  #calculate and save results
  colnames(temp_table) <- c("genus","species","sex","ID","type","length","region")
  temp_table$length <- as.numeric(temp_table$length)
  #save RAW length in .txt file 
  temp_name <- paste0("Rblength_lat_",genus,"_",ID,"_",sex,"_",type,"_",region) 
  assign(temp_name,temp_table)
  
  filename <- paste0(wd,"/save/",temp_name, ".txt")
  write.table(get(temp_name), file = filename, sep = "\t", quote = FALSE, row.names = FALSE)
  
  #take means from dorsal and ventral retina + sd and add to sepperate table 
  mean <- mean(temp_table$length)
  sd <- sd(temp_table$length)
  min <- min(temp_table$length)
  max <- max(temp_table$length)
  n <- length(temp_table$length)
  temp_row <- c(genus,species,sex,ID,type,mean,sd, min, max,n,region)
  Rblength_mean_sd <- rbind(Rblength_mean_sd,temp_row)
  #delete temp table
  rm(temp_table)
  } else {
    print(paste0(i,"_vt.txt is missing"))}
  
  
  
  i <- i+1
  
}

colnames(Rblength_mean_sd) <- c("Genus","Species","Sex","ID","Type","mean","sd","min","max","n","Region")
colnames(length_master) <- c("Genus","Species","Sex","ID","Type","length","Region")
length_master$length <- as.numeric(length_master$length)
Rblength_mean_sd$mean <- as.numeric(Rblength_mean_sd$mean )
Rblength_mean_sd$sd <- as.numeric(Rblength_mean_sd$sd )
Rblength_mean_sd$min <- as.numeric(Rblength_mean_sd$min)
Rblength_mean_sd$max <- as.numeric(Rblength_mean_sd$max)
Rblength_mean_sd <- Rblength_mean_sd %>%
  mutate(
    mean = round(mean, 1),
    sd = round(sd, 1),
    min = round(min, 1),
    max = round(max, 1)
  )
write.table(length_master, file = paste0(wd,"/save/length_master.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(Rblength_mean_sd, file = paste0(wd,"/save/Rblength_mean_sd.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

colnames(dist_master) <- c("Genus","Species","Sex","ID","Type","mean_dist","sd_dist")

write.table(dist_master, file = paste0(wd,"/save/dist_master.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
}



####------------------------------------------------------------------####
####Table for publication####
pb <- Rblength_mean_sd
pb <- pb %>%
  # Zuerst nach Genus und Sex gruppieren, um sicherzustellen, dass die Daten korrekt organisiert werden
  group_by(Type,Region) %>%
  # Pivot-Wide-Transformation
  pivot_wider(
    names_from = c( Type,Region),  # "type" wird verwendet, um die Spalten zu erstellen
    values_from = c(mean,sd,min,max,n),  # "td" und "Mod" werden in den Zellen eingefügt
    names_sep = "_"  # Spaltenname wird kombiniert mit einem Unterstrich, z.B. td_f_type1
  ) %>%
  # Spalte sex beibehalten
  ungroup()
genus_order <- c("Ampulex","Sceliphron", "Ammophila","Odontosphex","Heliocausus","Bembix","Bicyrtes","Microbembex", "Stictiella","Liris","Tachysphex","Tachytes")
pb$Genus <- factor(pb$Genus, levels = genus_order)
pb <- pb %>%
  arrange(Genus)
write.table(pb,paste0(wd,"/RBL_Summary_pb.txt"),sep = "\t", quote = FALSE, row.names = FALSE )
####------------------------------------------------------------------####
Rblength_mean_sd$mean <- as.numeric(Rblength_mean_sd$mean)
Rblength_mean_sd$sd <- as.numeric(Rblength_mean_sd$sd)


genus_order <- c("Ampulex","Sceliphron", "Ammophila","Odontosphex","Heliocausus","Bembix","Bicyrtes","Microbembex", "Stictiella","Liris","Tachysphex","Tachytes")
genus_order <- rev(genus_order)
Rblength_mean_sd$Genus <- factor(Rblength_mean_sd$Genus, levels = genus_order)

f <- filter(Rblength_mean_sd, Sex =="f")
m <- filter(Rblength_mean_sd, Sex =="m")


LO <- filter(f, Type == "LO")
MO <- filter(f, Type == "MO")

LO_vt <- filter(LO, Region == "vt")
LO_do <- filter(LO, Region == "do")
MO_vt <- filter(MO, Region == "vt")
MO_do <- filter(MO, Region == "do")

female_length<- ggplot(data = Rblength_mean_sd, aes(x = Genus)) +
  scale_x_discrete(drop = FALSE) + 
  scale_y_continuous(limits = c(0, 50)) +
  
  geom_bar(data = LO_do, aes(x = Genus, y = mean), stat = "identity", fill = "#4F4F4F", position = position_nudge(x = -0.325, y = 0), width = 0.2) +  
  geom_bar(data = MO_do, aes(x = Genus, y = mean), stat = "identity", fill = "#4F4F4F", position = position_nudge(x = 0.125, y = 0), width = 0.2) +  
  geom_bar(data = MO_vt, aes(x = Genus, y = mean), stat = "identity", fill = "#A9A9A9", position = position_nudge(x = 0.325, y = 0), width = 0.2) +
  geom_bar(data = LO_vt, aes(x = Genus, y = mean), stat = "identity", fill = "#A9A9A9", position = position_nudge(x = -0.125, y = 0), width = 0.2) +  
  
  
  geom_errorbar(data = LO_do, aes(x = Genus, ymin = mean-sd, ymax = mean+sd), position = position_nudge(x = -0.325, y = 0), width = 0)+
  geom_errorbar(data = MO_do, aes(x = Genus, ymin = mean-sd, ymax = mean+sd), position = position_nudge(x = 0.125, y = 0), width = 0) +  
  geom_errorbar(data = MO_vt, aes(x = Genus, ymin = mean-sd, ymax = mean+sd), position = position_nudge(x = 0.325, y = 0), width = 0) +
  geom_errorbar(data = LO_vt, aes(x = Genus, ymin = mean-sd, ymax = mean+sd), position = position_nudge(x = -0.125, y = 0), width = 0) +  
  
  # Textbeschriftung für MO
  geom_text(data = MO_do, aes(x = Genus, y = 0.002, label = "MO"), hjust = 1.2, size = 5, color = "black", position = position_nudge(0.225)) +  
  # Textbeschriftung für LO
  geom_text(data = LO_do, aes(x = Genus, y = 0.002, label = "LO"), hjust = 1.2, size = 5, color = "black", position = position_nudge(-0.225)) +  
  
  coord_flip() +  # Balken waagerecht drehen
  labs(y = "mean rhabdom length (lateral) in µm", x = "") +
  theme(
    axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
    axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
    panel.background = element_rect(fill = "white"),  # Hintergrund des Panels auf Weiß setzen
    plot.background = element_rect(fill = "white"),  # Hintergrund des Plots auf Weiß setzen
    panel.grid.major = element_blank(),  # Hauptgitterlinien entfernen
    panel.grid.minor = element_blank(),  # Nebengitterlinien entfernen
    axis.title = element_text(color = "black"),  # Achsentitel in Schwarz
    axis.text = element_text(color = "black", size = 15),  # Achsentext in Schwarz
    axis.ticks = element_line(color = "black"),  # Achsenmarkierungen in Schwarz
    legend.background = element_rect(fill = "white"),  # Hintergrund der Legende auf Weiß setzen
    legend.text = element_text(color = "black", size = 20),  # Legendentext in Schwarz
    legend.title = element_text(color = "black", size = 20),  # Legendentitel in Schwarz
    plot.title = element_text(color = "black", hjust = 0.5, size = 30),  # Plot-Titel in Schwarz
    plot.subtitle = element_text(color = "black", hjust = 0.5),  # Untertitel in Schwarz
    plot.caption = element_text(color = "black"),  # Bildunterschrift in Schwarz
    axis.line = element_line(color = "black")  # Achsenlinien in Schwarz
  ) + 
  ggtitle("female")
ggsave(paste0(wd,"/","female_length",".png"), width = 12, height = 18, plot = female_length, dpi = 150)
  

LO <- filter(m, Type == "LO")
MO <- filter(m, Type == "MO")

LO_vt <- filter(LO, Region == "vt")
LO_do <- filter(LO, Region == "do")
MO_vt <- filter(MO, Region == "vt")
MO_do <- filter(MO, Region == "do")

male_length<- ggplot(data = Rblength_mean_sd, aes(x = Genus)) +
  scale_x_discrete(drop = FALSE) + 
  scale_y_continuous(limits = c(0, 50)) +
  
  geom_bar(data = LO_do, aes(x = Genus, y = mean), stat = "identity", fill = "#4F4F4F", position = position_nudge(x = -0.325, y = 0), width = 0.2) +  
  geom_bar(data = MO_do, aes(x = Genus, y = mean), stat = "identity", fill = "#4F4F4F", position = position_nudge(x = 0.125, y = 0), width = 0.2) +  
  geom_bar(data = MO_vt, aes(x = Genus, y = mean), stat = "identity", fill = "#A9A9A9", position = position_nudge(x = 0.325, y = 0), width = 0.2) +
  geom_bar(data = LO_vt, aes(x = Genus, y = mean), stat = "identity", fill = "#A9A9A9", position = position_nudge(x = -0.125, y = 0), width = 0.2) +  
  geom_errorbar(data = LO_do, aes(x = Genus, ymin = mean-sd, ymax = mean+sd), position = position_nudge(x = -0.325, y = 0), width = 0)+
  geom_errorbar(data = MO_do, aes(x = Genus, ymin = mean-sd, ymax = mean+sd), position = position_nudge(x = 0.125, y = 0), width = 0) +  
  geom_errorbar(data = MO_vt, aes(x = Genus, ymin = mean-sd, ymax = mean+sd), position = position_nudge(x = 0.325, y = 0), width = 0) +
  geom_errorbar(data = LO_vt, aes(x = Genus, ymin = mean-sd, ymax = mean+sd), position = position_nudge(x = -0.125, y = 0), width = 0) +  
  
  # Textbeschriftung für MO
  geom_text(data = MO_do, aes(x = Genus, y = 0.002, label = "MO"), hjust = 1.2, size = 5, color = "black", position = position_nudge(0.225)) +  
  # Textbeschriftung für LO
  geom_text(data = LO_do, aes(x = Genus, y = 0.002, label = "LO"), hjust = 1.2, size = 5, color = "black", position = position_nudge(-0.225)) +  
  
  coord_flip() +  # Balken waagerecht drehen
  labs(y = "mean rhabdom length (lateral) in µm", x = "") +
  theme(
    axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
    axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
    panel.background = element_rect(fill = "white"),  # Hintergrund des Panels auf Weiß setzen
    plot.background = element_rect(fill = "white"),  # Hintergrund des Plots auf Weiß setzen
    panel.grid.major = element_blank(),  # Hauptgitterlinien entfernen
    panel.grid.minor = element_blank(),  # Nebengitterlinien entfernen
    axis.title = element_text(color = "black"),  # Achsentitel in Schwarz
    axis.text = element_text(color = "black", size = 15),  # Achsentext in Schwarz
    axis.ticks = element_line(color = "black"),  # Achsenmarkierungen in Schwarz
    legend.background = element_rect(fill = "white"),  # Hintergrund der Legende auf Weiß setzen
    legend.text = element_text(color = "black", size = 20),  # Legendentext in Schwarz
    legend.title = element_text(color = "black", size = 20),  # Legendentitel in Schwarz
    plot.title = element_text(color = "black", hjust = 0.5, size = 30),  # Plot-Titel in Schwarz
    plot.subtitle = element_text(color = "black", hjust = 0.5),  # Untertitel in Schwarz
    plot.caption = element_text(color = "black"),  # Bildunterschrift in Schwarzp
    axis.line = element_line(color = "black")  # Achsenlinien in Schwarz
  ) + 
  ggtitle("male")
ggsave(paste0(wd,"/","male_length",".png"), width = 12, height = 18, plot = male_length, dpi = 150)


#### visuallisation####
Rblength_mean_sd$mean <- as.numeric(Rblength_mean_sd$mean)
Rblength_mean_sd$sd <- as.numeric(Rblength_mean_sd$sd)
#difference Ventral/dorsal

mean_by_genus_region <- Rblength_mean_sd %>%
  group_by(Genus, Region, Type) %>%        
  summarize(mean_value = mean(mean)) 
do <- filter(mean_by_genus_region, Region == "do")
vt <- filter(mean_by_genus_region, Region == "vt")
             

LO <- filter(mean_by_genus_region, Type == "LO")

ggplot(data = LO, aes(x = reorder(Genus, mean_value), y = mean_value,fill = Region)) +
  geom_jitter( width = 0.2, shape = 21, size = 5, color = "black")+
  theme(
    axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
    axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
    panel.background = element_rect(fill = "white"),  # Hintergrund des Panels auf Weiß setzen
    plot.background = element_rect(fill = "white"),  # Hintergrund des Plots auf Weiß setzen
    panel.grid.major = element_blank(),  # Hauptgitterlinien entfernen
    panel.grid.minor = element_blank(),  # Nebengitterlinien entfernen
    axis.title = element_text(color = "black"),  # Achsentitel in Schwarz
    axis.text = element_text(color = "black", size = 30),  # Achsentext in Schwarz
    axis.ticks = element_line(color = "black"),  # Achsenmarkierungen in Schwarz
    legend.background = element_rect(fill = "white"),  # Hintergrund der Legende auf Weiß setzen
    legend.text = element_text(color = "black", size = 20),  # Legendentext in Schwarz
    legend.title = element_text(color = "black", size = 20),  # Legendentitel in Schwarz
    plot.title = element_text(color = "black", hjust = 0.5),  # Plot-Titel in Schwarz
    plot.subtitle = element_text(color = "black", hjust = 0.5),  # Untertitel in Schwarz
    plot.caption = element_text(color = "black"),  # Bildunterschrift in Schwarz
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1)# Achsenlinien in Schwarz
  ) + 
  labs(x = "genus", y = "mean rhabdom length in µm")+
  scale_fill_manual(values = c("vt" = "#FFFFB3", "do" = "#B2D8B2"))
  

MO <- filter(mean_by_genus_region, Type == "MO")

ggplot(data = MO, aes(x = reorder(Genus, mean_value), y = mean_value,fill = Region)) +
  geom_jitter( width = 0.2, shape = 21, size = 5, color = "black")+
  theme(
    axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
    axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
    panel.background = element_rect(fill = "white"),  # Hintergrund des Panels auf Weiß setzen
    plot.background = element_rect(fill = "white"),  # Hintergrund des Plots auf Weiß setzen
    panel.grid.major = element_blank(),  # Hauptgitterlinien entfernen
    panel.grid.minor = element_blank(),  # Nebengitterlinien entfernen
    axis.title = element_text(color = "black"),  # Achsentitel in Schwarz
    axis.text = element_text(color = "black", size = 30),  # Achsentext in Schwarz
    axis.ticks = element_line(color = "black"),  # Achsenmarkierungen in Schwarz
    legend.background = element_rect(fill = "white"),  # Hintergrund der Legende auf Weiß setzen
    legend.text = element_text(color = "black", size = 20),  # Legendentext in Schwarz
    legend.title = element_text(color = "black", size = 20),  # Legendentitel in Schwarz
    plot.title = element_text(color = "black", hjust = 0.5),  # Plot-Titel in Schwarz
    plot.subtitle = element_text(color = "black", hjust = 0.5),  # Untertitel in Schwarz
    plot.caption = element_text(color = "black"),  # Bildunterschrift in Schwarz
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1)# Achsenlinien in Schwarz
  ) + 
  labs(x = "genus", y = "mean rhabdom length in µm")+
  scale_fill_manual(values = c("vt" = "#FFFFB3", "do" = "#B2D8B2"))
####Seperated by sex and Type####    

LOf <- filter(Rblength_mean_sd, Type == "LO" & Sex == "f")
LOf_g <- ggplot(data = LOf, aes(x = Genus, y = mean,fill = Region)) +
  geom_point( width = 0.2, shape = 21, size = 5, color = "black")+
  theme(
    axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
    axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
    panel.background = element_rect(fill = "white"),  # Hintergrund des Panels auf Weiß setzen
    plot.background = element_rect(fill = "white"),  # Hintergrund des Plots auf Weiß setzen
    panel.grid.major = element_blank(),  # Hauptgitterlinien entfernen
    panel.grid.minor = element_blank(),  # Nebengitterlinien entfernen
    axis.title = element_text(color = "black"),  # Achsentitel in Schwarz
    axis.text = element_text(color = "black", size = 15),  # Achsentext in Schwarz
    axis.ticks = element_line(color = "black"),  # Achsenmarkierungen in Schwarz
    legend.background = element_rect(fill = "white"),  # Hintergrund der Legende auf Weiß setzen
    legend.text = element_text(color = "black", size = 20),  # Legendentext in Schwarz
    legend.title = element_text(color = "black", size = 20),  # Legendentitel in Schwarz
    plot.title = element_text(color = "black", hjust = 0.5, size = 30),  # Plot-Titel in Schwarz, Größe 30
    plot.subtitle = element_text(color = "black", hjust = 0.5),  # Untertitel in Schwarz
    plot.caption = element_text(color = "black"),  # Bildunterschrift in Schwarz
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1)  # X-Achsen-Beschriftungen schräg
  )+
  labs(x = "genus", y = "mean rhabdom length in µm")+
  scale_fill_manual(values = c("vt" = "#FFFFB3", "do" = "#B2D8B2"))+
  ggtitle("LO female")
ggsave(paste0(wd,"/","LO female",".png"), width = 18, height = 12, plot = LOf_g, dpi = 300)



LOm <- filter(Rblength_mean_sd, Type == "LO"& Sex == "m")
LOm_g <- ggplot(data = LOm, aes(x = Genus, y = mean,fill = Region)) +
  geom_point( width = 0.2, shape = 21, size = 5, color = "black")+
  theme(
    axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
    axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
    panel.background = element_rect(fill = "white"),  # Hintergrund des Panels auf Weiß setzen
    plot.background = element_rect(fill = "white"),  # Hintergrund des Plots auf Weiß setzen
    panel.grid.major = element_blank(),  # Hauptgitterlinien entfernen
    panel.grid.minor = element_blank(),  # Nebengitterlinien entfernen
    axis.title = element_text(color = "black"),  # Achsentitel in Schwarz
    axis.text = element_text(color = "black", size = 15),  # Achsentext in Schwarz
    axis.ticks = element_line(color = "black"),  # Achsenmarkierungen in Schwarz
    legend.background = element_rect(fill = "white"),  # Hintergrund der Legende auf Weiß setzen
    legend.text = element_text(color = "black", size = 20),  # Legendentext in Schwarz
    legend.title = element_text(color = "black", size = 20),  # Legendentitel in Schwarz
    plot.title = element_text(color = "black", hjust = 0.5, size = 30),  # Plot-Titel in Schwarz, Größe 30
    plot.subtitle = element_text(color = "black", hjust = 0.5),  # Untertitel in Schwarz
    plot.caption = element_text(color = "black"),  # Bildunterschrift in Schwarz
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1)  # X-Achsen-Beschriftungen schräg
  )+ 
  labs(x = "genus", y = "mean rhabdom length in µm")+
  scale_fill_manual(values = c("vt" = "#FFFFB3", "do" = "#B2D8B2"))+
  ggtitle("LO male")
ggsave(paste0(wd,"/","LO male",".png"), width = 18, height = 12, plot = LOm_g, dpi = 300)

MOf <- filter(Rblength_mean_sd, Type == "MO"& Sex == "f")
MOf_g <- ggplot(data = MOf, aes(x = Genus, y = mean,fill = Region)) +
  geom_point( width = 0.2, shape = 21, size = 5, color = "black")+
  theme(
    axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
    axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
    panel.background = element_rect(fill = "white"),  # Hintergrund des Panels auf Weiß setzen
    plot.background = element_rect(fill = "white"),  # Hintergrund des Plots auf Weiß setzen
    panel.grid.major = element_blank(),  # Hauptgitterlinien entfernen
    panel.grid.minor = element_blank(),  # Nebengitterlinien entfernen
    axis.title = element_text(color = "black"),  # Achsentitel in Schwarz
    axis.text = element_text(color = "black", size = 15),  # Achsentext in Schwarz
    axis.ticks = element_line(color = "black"),  # Achsenmarkierungen in Schwarz
    legend.background = element_rect(fill = "white"),  # Hintergrund der Legende auf Weiß setzen
    legend.text = element_text(color = "black", size = 20),  # Legendentext in Schwarz
    legend.title = element_text(color = "black", size = 20),  # Legendentitel in Schwarz
    plot.title = element_text(color = "black", hjust = 0.5, size = 30),  # Plot-Titel in Schwarz, Größe 30
    plot.subtitle = element_text(color = "black", hjust = 0.5),  # Untertitel in Schwarz
    plot.caption = element_text(color = "black"),  # Bildunterschrift in Schwarz
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1)  # X-Achsen-Beschriftungen schräg
  )+
  labs(x = "genus", y = "mean rhabdom length in µm")+
  scale_fill_manual(values = c("vt" = "#FFFFB3", "do" = "#B2D8B2"))+
  ggtitle("MO female")
ggsave(paste0(wd,"/","MO female",".png"), width = 18, height = 12, plot = MOf_g, dpi = 300)

MOm <- filter(Rblength_mean_sd, Type == "MO"& Sex == "m")
MOm_g <- ggplot(data = MOm, aes(x = Genus, y = mean,fill = Region)) +
  geom_point( width = 0.2, shape = 21, size = 5, color = "black")+
  theme(
    axis.title.x = element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
    axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
    panel.background = element_rect(fill = "white"),  # Hintergrund des Panels auf Weiß setzen
    plot.background = element_rect(fill = "white"),  # Hintergrund des Plots auf Weiß setzen
    panel.grid.major = element_blank(),  # Hauptgitterlinien entfernen
    panel.grid.minor = element_blank(),  # Nebengitterlinien entfernen
    axis.title = element_text(color = "black"),  # Achsentitel in Schwarz
    axis.text = element_text(color = "black", size = 15),  # Achsentext in Schwarz
    axis.ticks = element_line(color = "black"),  # Achsenmarkierungen in Schwarz
    legend.background = element_rect(fill = "white"),  # Hintergrund der Legende auf Weiß setzen
    legend.text = element_text(color = "black", size = 20),  # Legendentext in Schwarz
    legend.title = element_text(color = "black", size = 20),  # Legendentitel in Schwarz
    plot.title = element_text(color = "black", hjust = 0.5, size = 30),  # Plot-Titel in Schwarz, Größe 30
    plot.subtitle = element_text(color = "black", hjust = 0.5),  # Untertitel in Schwarz
    plot.caption = element_text(color = "black"),  # Bildunterschrift in Schwarz
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1)  # X-Achsen-Beschriftungen schräg
  )+
  labs(x = "genus", y = "mean rhabdom length in µm")+
  scale_fill_manual(values = c("vt" = "#FFFFB3", "do" = "#B2D8B2"))+
  ggtitle("MO male")
ggsave(paste0(wd,"/","MO male",".png"), width = 18, height = 12, plot = MOm_g, dpi = 300)
