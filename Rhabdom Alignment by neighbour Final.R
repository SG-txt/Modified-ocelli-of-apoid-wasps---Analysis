library("dplyr")
library("tiff")
library("ggplot2")
library ("viridis")
library("RColorBrewer")
library(RANN)
library("skimr")
library ("qqplotr")
library ("rstatix")

{  

setwd("C:/Users/Stefan/Desktop/Dr/Analysen/Pol.Sens. R-Script/PolTest")
wd<-getwd()
if (dir.exists(paste0(wd,"/save/NN-Analysis"))){
  print("folder already exists")
} else {
  dir.create(paste0(wd,"/save/NN-Analysis"))}
if (dir.exists(paste0(wd,"/controlnn/"))){
  print("folder already exists")
} else {
  dir.create(paste0(wd,"/controlnn/"))}

# Liste aller Dateien mit .txt format im Zielordner
filelist <- list.files(pattern = "*.txt")

#How many neigbours should be used (always intended number +1 for first line)
nnn <- 20
#Tabelle für Gesamtergebnisse
Dev_master <- data.frame()

Zähler <-1

for (file in filelist) {
  
  
  file_name <- filelist [Zähler]
  RAW_temp<- read.table(file_name, header = TRUE, skip = 7)
  
  
  
  # species sollte in der Ergebnisstabelle als Spalte auftauchen, ebenso ID
  table <- read.table(file_name, header = FALSE, nrows = 1, sep = ",")
  genus <- table [,1]
  species <- table [,2]
  ID <- table [,3]
  type <- table [,4]
  sblength <- table [,5]
  sex <- table [,6]
  mod <- table [,7]
  complete <-table [,8]
  
  corner1<- read.table(file_name, header = FALSE, nrows = 1, skip = 1)
  corner2<- read.table(file_name, header = FALSE, nrows = 1, skip = 2)
  corner3<- read.table(file_name, header = FALSE, nrows = 1, skip = 3)
  corner4<- read.table(file_name, header = FALSE, nrows = 1, skip = 4)
  
  sbpoint1<- read.table(file_name, header = FALSE, nrows = 1, skip = 5)
  sbpoint2<- read.table(file_name, header = FALSE, nrows = 1, skip = 6)
  
  #calculate factor to convert distance between points to micrometer 
  sbfactor <- sblength/sqrt((sbpoint2[1,1]-sbpoint1[1,1])^2+(sbpoint2[1,2]-sbpoint1[1,2])^2)
  

  #create RAW and populate with data from RAW_temp 
  RAW2 <- data.frame(x = RAW_temp$X, y = RAW_temp$Y, time = RAW_temp$time)
  if (nrow(RAW2) %% 5 != 0) {
    stop("Number of points is not evenly dividable by 5")
  }
  ####Shorten RAW-table to only include the first and last point of each five Rhabdom pair####
  filtered_df <- data.frame()
  for (i in seq(1, nrow(RAW2), by = 5)) {
    
    filtered_df <- rbind(filtered_df, RAW2[i, ])
    filtered_df <- rbind(filtered_df, RAW2[i+4, ])
    
  }

  #x and y coordinates of each of the two end points of a rhabdom have to be in the same row as x1,y1,x2,y2
  transposed_df <- data.frame(x1 = filtered_df[1,1], y1 = filtered_df[1,2], x2 = filtered_df[2,1], y2 = filtered_df[2,2])
  
  i <- 3
  while (i<nrow(filtered_df)) {
    temp_df2 <- data.frame(x1 = filtered_df[i,1], y1 = filtered_df[i,2], x2 = filtered_df[i+1,1], y2 = filtered_df[i+1,2])
    transposed_df <- rbind(transposed_df,temp_df2)
    i <- i+2
  }
  
  #X1/2 SHOULD NEVER BE THE EXACT SAME!!!!!!!!!!!!!!!!!!!!!!!!
  #CHECK FOR SAME X1/2 AND ADD A MINIMAL DEVIATION TO X2 IN THIS CASE 

  
 for (s in 1:nrow(transposed_df)) {
   
   if (transposed_df[s,1]==transposed_df[s,3]){c(transposed_df[s,3] <- transposed_df[s,3]+0.001,print(paste0("Warning:TRUE",s," -> x1/2 with identical value, 0.001 has been added to x2" )))}else{}
   
 }
  
  
   # calculate midpoints for each rhabdom and unite with coordinates 
  midpoints <- data.frame( xmid = ((transposed_df$x1 + transposed_df$x2)/2), ymid = ((transposed_df$y1 + transposed_df$y2/2)))
    Comb <- cbind(transposed_df,midpoints)
    
#calculate nearest neighbour from midpoints 
    nn_results <- nn2(Comb[,c(5,6)], k=nnn)
    nn_id <- nn_results$nn.idx
    
    #ggplot to access how often every ID was used in the nncomparison 
    # idealy I would like to have an equal use of each ID?!
    häufigkeit <- table(nn_id)
    Hdf <- as.data.frame(häufigkeit)
    Hdf <- cbind(Hdf, midpoints)
    control <- ggplot(data = Hdf,  mapping = aes(x=xmid, y=ymid)) +
      annotate("text",label = paste0( "mean Freq: ", mean(Hdf$Freq),"  -  STD: ",sd(Hdf$Freq),"\n","mean Freq/",nnn,"(#nn): ", mean(Hdf$Freq)/nnn, "  -  STD: ", sd(Hdf$Freq)/nnn), x=2,y=2)+
      geom_point(aes(colour = Freq))
    ggsave(paste0(wd,"/controlnn/",ID,genus,type,sex,"control.png"), width = 9,height = 6, plot = control, dpi = 100) 
    controlH <- ggplot(data = Hdf,  mapping = aes(y=Freq)) +
      geom_histogram()
    ggsave(paste0(wd,"/controlnn/",ID,genus,type,sex,"controlH.png"), width = 9,height = 6, plot = controlH, dpi = 100) 
    controlS <- ggplot(data=Hdf, mapping = aes(x=nn_id,y=Freq))+
      geom_point()
    ggsave(paste0(wd,"/controlnn/",ID,genus,type,sex,"controlS.png"), width = 9,height = 6, plot = controlS, dpi = 100)
# analyse varianz of Freq/RHabdom used to see how effect varies between #nn? 
    
# calculate the angel between the first line of each matrix and its five neighbours 
  # referenz matrixes:   nn_id[ matrix nr 1, first value within matrix 1] 
    #-> We have as many matrixes as rhabdomers and as many values within each matrix as chosen nnn
    
   
    
  dev_table <- data.frame()
   l <- 1 
  for (l in 1:nrow(Comb)) {
    Incl_table <- data.frame()
    # calculate incline for first rhabdom y2-y1/x2-x1
    Incl0 <- (Comb[nn_id[l,1],4]-Comb[nn_id[l,1],2])/(Comb[nn_id[l,1],3]-Comb[nn_id[l,1],1])
    Incl_table <- rbind(Incl_table,Incl0)
    for (j in 2:nnn) {
    #calculate incline for other rhabdoms 
      Inclnnn <- (Comb[nn_id[l,j],4]-Comb[nn_id[l,j],2])/(Comb[nn_id[l,j],3]-Comb[nn_id[l,j],1]) 
      Incl_table <- rbind(Incl_table,Inclnnn)
        }
    # Angels between inlcine of R1 and all other rhabdoms have to be calculated 
        deg_table <- data.frame()
    
    for (m in 2:nnn) {
      angel_rad <- atan((Incl_table[m,1]-Incl_table[1,1])/(1+Incl_table[m,1]*Incl_table[1,1]))
      angel_deg <- angel_rad * 180 / pi
      deg_table <- rbind(deg_table,angel_deg)
    }
    
    # mean of angels requires the absolute value, this should avoid that negative and positive inclines add up to zero, instead of their combined value
    deg_table <- abs(deg_table)
    mean_dev <- mean(deg_table[,1])
    dev_table <- rbind(dev_table,c(mean_dev,l,ID,type,genus,species,sex,mod,complete))
     
    
  }
  
   # assign ID name to dev_table
   #save as .txt
   names(dev_table) <- c("DegA","Nr","ID","type","genus","species","sex","mod","complete")
   dev_table <- cbind(dev_table,Comb)
   #temp_name <- paste0("Dev_table_",ID)  
   filename <- paste0(wd,"/save/NN-Analysis/NN_DEV",ID,type, ".txt")
   write.table(dev_table, file = filename, sep = "\t", quote = FALSE, row.names = FALSE)
   
   Dev_master <- rbind(Dev_master,dev_table)
   Dev_master$DegA <- as.numeric(Dev_master$DegA)
   
   
   filename_m <- paste0(wd,"/save/Alg_master.txt")
   write.table(Dev_master, file = filename_m, sep = "\t", quote = FALSE, row.names = FALSE)
 Zähler <- Zähler+1
}  
}

#x1/2with identical value, 0.001 added to x2 in:
#three cases
# 45,19,116