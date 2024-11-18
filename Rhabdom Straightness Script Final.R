
#START
#prepare txt-files: Insert one line at the very top in the following format: Gattung, Art, ID
#The next four lines are the coordinates of the four corners, so only two numbers
#another two lines of coordinates for the scalebar, not yet used, might be NA
#IMPORTANT: The "counter" for the number of points has no header in 8 two.
#Thereby it will not be added as a column in the data frame. Add the header "time" manually to the .txt file over the first column (total column nr of data frame should be seven)
#TXT-FILENAME IS NOT IMPORTANT
#The image the points where put in should be saved as a .tif file in the workdirectory
#IMAGE NAME IS IMPORTANT, SHOULD BE: "ID"points
#IMAGE SHOULD HAVE 400 DPI - or adjust DPI in the plot section: map of aligment plot
####Library prep####
{if (!requireNamespace(c("tiff","dplyr","ggplot2","viridis","RColorBrewer"), quietly = TRUE)) {
  install.packages(c("tiff","dplyr","ggplot2","viridis","RColorBrewer"))}
  
  library("dplyr")
  library("tiff")
  library("ggplot2")
  library ("viridis")
  library("RColorBrewer")
  
}
####Loop containing all calculations and plots for one specimen####
setwd("C:/Users/Stefan/Desktop/Dr/Analysen/Pol.Sens. R-Script/PolTest")
wd<-getwd()


{
  Str_master <- data.frame()
  ####SetWD, disect input .txt, set-up and control folder structure
  
  # Liste aller Dateien mit .txt format im Zielordner
  filelist <- list.files(pattern = "*.txt")
  #filelist <- filelist[4]
  #DataFrame für Endergebnisse
  Results_Species <-data.frame(Column1= character(),Column2=numeric(), Column3=numeric(), Column4= numeric(), Column5=numeric(), COlumn6=numeric(), Column7=numeric())
  
  
  
  Zähler <- 1
  for (file in filelist) {
    
    
    file_name <- filelist [Zähler]
    #if (dir.exists(paste0(wd,file_name)) = TRUE){
    # print("already processed")
    #} else {
    RAW_temp<- read.table(file_name, header = TRUE, skip = 7)
    
    # species sollte in der Ergebnisstabelle als Spalte auftauchen, ebenso ID
    table <- read.table(file_name, header = FALSE, nrows = 1, sep = ",")
    genus <- table [,1]
    species <- table [,2]
    ID <- table [,3]
    type <- table [,4]
    sblength <- table [,5]
    sex <- table [,6]
    mod <- table[,7]
    complete <-table [,8]
    
    corner1<- read.table(file_name, header = FALSE, nrows = 1, skip = 1)
    corner2<- read.table(file_name, header = FALSE, nrows = 1, skip = 2)
    corner3<- read.table(file_name, header = FALSE, nrows = 1, skip = 3)
    corner4<- read.table(file_name, header = FALSE, nrows = 1, skip = 4)
    
    sbpoint1<- read.table(file_name, header = FALSE, nrows = 1, skip = 5)
    sbpoint2<- read.table(file_name, header = FALSE, nrows = 1, skip = 6)
    
    #calculate factor to convert distance between points to micrometer 
    sbfactor <- sblength/sqrt((sbpoint2[1,1]-sbpoint1[1,1])^2+(sbpoint2[1,2]-sbpoint1[1,2])^2)
    
    
    #if necessary create the folder structure and respective folders
    if (dir.exists(paste0(wd,"/save/"))){
      print("folder already exists")
    } else {
      dir.create(paste0(wd,"/save/"))}
    if (dir.exists(paste0(wd,"/save/",ID,type))){
      
      print("folder already exists")
    } else {
      dir.create(paste0(wd,"/save/",ID,type))}
    if (dir.exists(paste0(wd, "/",ID,type))){
      print("folder already exists")
    } else {
      dir.create(paste0(wd, "/",ID,type))}
    
    #create RAW and populate with data from RAW_temp 
    RAW <- data.frame(x = RAW_temp$X, y = RAW_temp$Y, time = RAW_temp$time)
   
    ####Zerlegen in Einzeltabellen für Rhabdomere ####
    # Anzahl der Einträge pro kleiner Tabelle
    split_size <- 5
    
    # Kontrolle: RAW muss glatt durch split_size teilbar sein, sonst skript abgebrochen 
    if (nrow(RAW) %% split_size != 0) {
      stop("Number of points is not evenly dividable by the split_size")
    }
    
    # Schleife zum Zerlegen der großen Tabelle
    for (i in seq(1, nrow(RAW), split_size)) {
      # Erstellen einer kleinen Tabelle
      small_table <- RAW[i:(i+split_size-1), ]
      # Benennen der kleinen Tabelle
      small_table_name <- paste0("Rhabdomer", (i-1)/split_size + 1)
      assign(small_table_name, small_table)
    }
    
    #### Ergebnistabellen anlegen ####
    # Dataframe für Ergebnisse der Rhabdom straightnes 
    Results <- data.frame(row.names = c("Min.","1st_Qu.", "Median", "Mean", "3rd_Qu.", "Max."))
    # Dataframe für die Gesamtausrichtung der Rhabdomere einer Gattung nach einzelnen Rhabdomeren 
    Rhabdom_Alignment <- data.frame(Colum1 = character(), Column2 = numeric())
    #data frame erstellen der alle Abweichungswinkel mit korrektem Vorzeichen enthält
    RAW_Angles <- data.frame()
    ####Overall loop, containing all calculations for one specimen####
    y<-1
    while (y<=(length(RAW$x)/split_size)) { 
      temp_name <- paste0("Rhabdomer",y)
      Temp_Tabel <- get(temp_name)
      incl0 <- (Temp_Tabel[5,2]-(Temp_Tabel[1,2]+0.000001))/(Temp_Tabel[5,1]-(Temp_Tabel[1,1]+0.000001))
      incl1 <- (Temp_Tabel[2,2]-(Temp_Tabel[1,2]+0.000001))/(Temp_Tabel[2,1]-(Temp_Tabel[1,1]+0.000001))
      incl2 <- (Temp_Tabel[3,2]-(Temp_Tabel[2,2]+0.000001))/(Temp_Tabel[3,1]-(Temp_Tabel[2,1]+0.000001))
      incl3 <- (Temp_Tabel[4,2]-(Temp_Tabel[3,2]+0.000001))/(Temp_Tabel[4,1]-(Temp_Tabel[3,1]+0.000001))
      incl4 <- (Temp_Tabel[5,2]-(Temp_Tabel[4,2]+0.000001))/(Temp_Tabel[5,1]-(Temp_Tabel[4,1]+0.000001))
      #speichern des Rhabdom_Straightn table 
      #temp_name6 <- paste0("Straightn_per_Rh_", ID,"Rh", y)
      #assign(temp_name6,c(incl0,incl1,incl2,incl3,incl4))
      #zusätzlich speichern als .txt datei
      #filename2 <- paste0("C:/Users/Stefa/Desktop/Test Ocelli Straightness/save/",temp_name6, ".txt")
      #write.table(get(temp_name6), file = filename2, sep = "\t", quote = FALSE, row.names = TRUE)
      
      
      #####Für Rhabdomer-overall alignment ####
      angle_overall<- atan((incl0-0)/(1+incl0*0))*180/pi# hier nicht den betrag verwenden! kein abs()! Denn -80 Grad ist nicht gleich ausgerichtet wie +80 Grad!
     
      temp_name2 <- paste0("overall",y)
      assign(temp_name2,angle_overall)
      Rhabdom_Alignment <- rbind(Rhabdom_Alignment, get(paste0("overall",y)))
      colnames(Rhabdom_Alignment) <- "Dev_from_horizon_in_degree"
      #speichern des Rhabdom_Alignment table 
      temp_name4 <- paste0("Alignment_per_Rh_", ID)
      assign(temp_name4,Rhabdom_Alignment)
      #zusätzlich speichern als .txt datei
      filename1 <- paste0(wd,"/save/",ID,type,"/",temp_name4, ".txt")
      write.table(get(temp_name4), file = filename1, sep = "\t", quote = FALSE, row.names = FALSE)
      
      
      
      #####nächster loop berechnet aus allen inclines die jeweiligen Winkel####
      x<-1
      while (x<split_size) {
        var_name <- paste0("incl",x)
        current_varriabel <- get(var_name)
        angle_rad <- atan((current_varriabel-incl0)/(1+current_varriabel*incl0)) 
        name <- paste0("angle",x)
        assign(name,angle_rad * 180 / pi)
        x <- x+1
        #speichern der angle_rad zur Kontrolle
        #temp_name7 <- paste0("Angle_per_Rh_", ID,"Rh", y)
        #assign(temp_name7,c(angle1,angle2,angle3,angle4))
        #zusätzlich speichern als .txt datei
        #filename3 <- paste0("C:/Users/Stefa/Desktop/Test Ocelli Straightness/save/",temp_name7, ".txt")
        #write.table(get(temp_name7), file = filename3, sep = "\t", quote = FALSE, row.names = TRUE)
      }
      #calculate Rhabdomer length
      distance1 <- sqrt((Temp_Tabel[2,1]-Temp_Tabel[1,1])^2+(Temp_Tabel[2,2]-Temp_Tabel[1,2])^2)
      distance2 <- sqrt((Temp_Tabel[3,1]-Temp_Tabel[2,1])^2+(Temp_Tabel[3,2]-Temp_Tabel[2,2])^2)
      distance3 <- sqrt((Temp_Tabel[4,1]-Temp_Tabel[3,1])^2+(Temp_Tabel[4,2]-Temp_Tabel[3,2])^2)
      distance4 <- sqrt((Temp_Tabel[5,1]-Temp_Tabel[4,1])^2+(Temp_Tabel[5,2]-Temp_Tabel[4,2])^2)
      
      Rhabdom_length <- sum(c(distance1,distance2,distance3,distance4))*sbfactor
      
      # assemble Result for this Rhabdomer, including Rhabdom length
      vector <- c(abs(angle1), abs(angle2), abs(angle3), abs(angle4))#hier abs() für Betrag -> Alle abweichungen müssen aufsummiert werde, +/- ist "egal"
      name_result <- paste0("R_straightn_",y) 
      Results <- rbind(Results, assign(name_result,c(unname(summary(vector)),Rhabdom_length)))
      colnames(Results) <- c("Min.","1st_Qu.", "Median", "Mean", "3rd_Qu.", "Max.","Length")
      
      
      
      
      row <- c(as.numeric(mean(vector)),y,Rhabdom_length,ID,type,genus,species,sex,mod )
      Str_master <- rbind(Str_master,row)
  
      
      #RAW Angles dieser iteration zum Plotten speichern 
      RAW_Angles <- rbind (RAW_Angles,angle1)
      RAW_Angles <- rbind (RAW_Angles,angle2)
      RAW_Angles <-rbind (RAW_Angles,angle3)
      RAW_Angles <-rbind (RAW_Angles,angle4)
      #RAW  Anlges Tabelle als .txt speichern
      
      temp_name5 <- paste0("RAW_Angles", ID)
      assign(temp_name5,RAW_Angles)
      filename4 <- paste0(wd,"/save/",ID,type,"/","RAW_Angles",ID,type, ".txt")
      colnames(RAW_Angles)[1] <- "dev"
      write.table(get(temp_name5), file = filename4, sep = "\t", quote = FALSE, row.names = FALSE)
      
     
      
      #speichern der Results table zur Kontrolle oder für weitere Rechnungen
      temp_name3 <- paste0("angle_Summary_per_Rh_", ID)
      assign(temp_name3,Results)
      #zusätzlich speichern als .txt datei
      filename2 <- paste0(wd,"/save/",ID,type,"/",temp_name3, ".txt")
      write.table(get(temp_name3), file = filename2, sep = "\t", quote = FALSE, row.names = FALSE)
      
      
      
      
      
      y <- y+1
    }
    
    # Hier werden Daten in den Results_species table gespeichert, hier kann also am einfachsten Statistik implementiert werden,
    #die mit jeder Itteration ausgeführt werden soll. 
    new_row <- data.frame(Column1 = genus, Column2 = species, column3 = ID, Column4 = mean(Results$Mean), Column5 = sd(Results$Mean), Column6 = mean(Rhabdom_Alignment$Dev_from_horizon_in_degree), Column7 = sd(Rhabdom_Alignment$Dev_from_horizon_in_degree))
    names(new_row) <- c("genus","species", "ID","mean_angle_rhabdomstrn","sd","mean_deviation_from_horizon","sd_deviation_from_horizon")
    Results_Species <- rbind(Results_Species, new_row[1,])
    #Spalten umbennen
    names(Results_Species) <- c("genus","species", "ID","mean_angle_rhabdomstrn","sd","mean_deviation_from_horizon","sd_deviation_from_horizon")
    #zusätzlich speichern als .txt datei
    filename3 <- paste0(wd,"/save/Results_Species.txt")
    write.table(Results_Species, file = filename3, sep = "\t", quote = FALSE, row.names = FALSE)
    
    ####All Plots for one specimen####
    #make a new data frame containing the coordinates to display each rhabdom as well as the deviation from horizon
    # for manual analysis simply insert correct Alignment data here. z.B. Alignment_per_RhG96
    {
      Alignment_chart <- get(paste0("Alignment_per_Rh_",ID))
      
      #new data frame for visualisation
      Visualisation_clgr <- data.frame(new_column = numeric(length(RAW$x)))
      Visualisation_clgr$x <- RAW$x
      Visualisation_clgr$y <- RAW$y
      Visualisation_clgr <- Visualisation_clgr[,!seq_along(RAW) %in% 1]
      # only keep the first and last point of every rhabdomer
      filtered_df <- data.frame()
      for (i in seq(1, nrow(Visualisation_clgr), by = 5)) {
        
        filtered_df <- rbind(filtered_df, Visualisation_clgr[i, ])
        filtered_df <- rbind(filtered_df, Visualisation_clgr[i+4, ])
        
      }
      Visualisation_clgr <- filtered_df
      # associate the correct deviation in degree to the two coordinates belonging to each rhabdom 
      # Erstellen des Ziel-Dataframes
      Visualisation_clgr <- data.frame(
        Dev = rep(abs(Alignment_chart$Dev_from_horizon_in_degree), each = 2),  # Jede Zeile aus df_origin wird dupliziert
        
        x = rep(Visualisation_clgr$x),
        y = rep(Visualisation_clgr$y)
      )
      #dataframe mit Horizont linie erstellen
      horizon <- data.frame(
        x = c(max(Visualisation_clgr$x),min(Visualisation_clgr$x)),
        y= c(mean(Visualisation_clgr$y),mean(Visualisation_clgr$y))
      )
      
      #Seitenverhältnis aus .tiff Datei auslesen 
      # Öffnen der TIFF-Datei
      img <- readTIFF(paste0(wd,"/",ID,type,"points.tif"))
      # Abmessungen und Seitenverhältnis ermitteln
      width <- dim(img)[2]
      height <- dim(img)[1]
      
      #visualisation_clgr vier weitere punkte mit den xy coordinaten aus den "corner" hinzufügen
      #Dev sollte NA sein, und NA in der colour scale als weiß definiert sein 
      
   
      
      colnames(corner1) <- c("x","y")
      colnames(corner2) <- c("x","y")
      colnames(corner3) <- c("x","y")
      colnames(corner4) <- c("x","y")
      Visualisation_clgr_2<- data.frame(Dev = c(0,0,0,0),x= c(corner1[1,1],corner2[1,1],corner3[1,1],corner4[1,1]),y= c(corner1[1,2],corner2[1,2],corner3[1,2],corner4[1,2]))
      #plot
      p <- ggplot(data = Visualisation_clgr,
                  mapping = aes(
                    x = x,
                    y = y,
                    colour = Dev
                  ))+
        #geom_point()+
        theme_classic()+
        theme(axis.line = element_blank(),
              axis.text.x = element_blank(),  # Achsenbeschriftungen ausblenden
              axis.title.x = element_blank(),  # X-Achsentitel ausblenden http://127.0.0.1:38167/graphics/plot_zoom_png?width=1911&height=823
              axis.text.y = element_blank(),  # Y-Achsenbeschriftungen ausblenden
              axis.title.y = element_blank(),
              axis.ticks = element_blank(),legend.key.size = unit(4, "lines"),
              panel.background = element_rect(fill = "black"),
              #legend.background = element_rect(fill = "black", color = "black"),
              #legend.text = element_text(color = "white", size = 60),
              #legend.key.height = unit(6, "cm"), legend.key.width = unit(3, "cm"),
              #legend.title = element_text(size = 60,margin = margin (t = -20, b = -20)),
              legend.position = "none",
              panel.border = element_rect(color = "black", fill = NA, linewidth = 4))+
        scale_x_continuous(expand=(c(0,0)))+ scale_y_continuous(expand=(c(0,0)))+
        #scale_color_manual(values = "Dev"="Deviation from Horizon in °")+
        labs(color="Devation in °")+
        geom_line(data = horizon, aes(x=x,y=y), colour = "grey",linewidth=4)
      
      i <- 1
      while(i<length(Visualisation_clgr_2$x)){
        
        v <-i+1
        p <- p+ geom_line(data = Visualisation_clgr_2[i:v, ],
                          aes(x=x,
                              y=y),
                          linewidth=1,
                          color = "black"
        )
        
        i<-i+2
      }
      
      
      i <- 1
      while(i<length(Visualisation_clgr$x)){
        
        v <-i+1
        p <- p+ geom_line(data = Visualisation_clgr[i:v, ],
                          aes(x=x,
                              y=y),
                          linewidth=6
        )
        
        i<-i+2
      }
      i<-1
      
      
      #p=p+scale_colour_gradientn(colors = rev(brewer.pal(8, "YlGnBu")), limits = c(0, 90), breaks = seq(0, 90, by = 10))
      p <- p+scale_color_viridis(na.value = "white", option="turbo", limits = c(0, 90),breaks = seq(0, 90, by = 10), begin = 0.1)
      p=p+scale_y_reverse()
      #p=p+scale_x_reverse()
      #p=p+theme(legend.position = "none")
      #Seitenverhältnisse auf Plot übetragen 
      #p=p + coord_fixed(ratio = width / height) #Gibt das Seitenverhältnis des Zeichenbereiches an, muss noch dynamisch vom Bild abgenommen werden!!!
      print(p)
      
      ggsave(paste0(wd,"/",ID,type,"/",ID,type,"overall_Alig_map.png"), height = height/400*5, width = width/400*5, plot = p, dpi = 150 , limitsize = FALSE) #Hier muss die entsprechende DPI eingesetzt werden, Input fotos sollten also immer 400 dpi haben
    }
    
    
    #verbinde visualisation_clgr mit angle_Summary_per_Rh_ID
    temp_Table2<-get(paste0("angle_Summary_per_Rh_",ID))
    temp_Table3<-Visualisation_clgr
    #Ecken entfernen
    temp_Table3 <- temp_Table3[1:(nrow(temp_Table3) - 4), ]
   
    #Mitte der Rhabdomere ermitteln und in Tabelle eintragen 
    temp_Table4 <- data.frame(x = numeric(0), Mean_DEV_RH = numeric(0),  y = numeric(0))
    
    o<-1
    j<-1
    while (o < length(temp_Table3$x)){
      
      new_row2 <- data.frame(x = (temp_Table3[o+1,2]-temp_Table3[o,2])/2+temp_Table3[o,2], Mean_DEV_RH = temp_Table2[j,3],y=1)
      temp_Table4<-rbind(temp_Table4,new_row2)
     o<-o+2
    j<-j+1
    }
    
     #Tabelle nach x-Positionen sortieren
    temp_Table4 <- temp_Table4[order(temp_Table4$x), ]
    
    #Spalte Z einfügen, sollte eine sequenz in Reihenfolge enthalten -> 1-length(temp_Table4)
    
    temp_Table4$z <- seq(1,length(temp_Table4$x))
    
  #plotten
    p5 <- ggplot(data = temp_Table4,
           mapping = aes(
             x = z,
             y = y,
             fill = Mean_DEV_RH))+
      geom_bar(stat = "identity")+
      scale_y_continuous(limits = c(0,1))+
      coord_fixed(ratio = 2)+
      scale_fill_viridis(option = "plasma", limits = c(0,90))+
      xlab("")+
      theme(legend.position = "none",
            panel.background = element_rect(fill = "black"),  # Weißer Hintergrund
            panel.grid = element_blank(),
            axis.title.x = element_text(size = 60),
            axis.text.x = element_blank(),  # X-Achsentitel ausblenden 
            axis.text.y = element_blank(),  # Y-Achsenbeschriftungen ausblenden
            axis.title.y = element_blank(),
            axis.ticks = element_blank())
      
    
    ggsave(paste0(wd,"/",ID,type,"/",ID,type,"strgh_within_map.png"), height = height/400*5, width = width/400*5, plot = p5, dpi = 150 , limitsize = FALSE) #Hier muss die entsprechende DPI eingesetzt werden, Input fotos sollten also immer 400 dpi haben
    
   
    #Balkendiagramm
    
    # für Balkendiagramm gruppieren
    
    #breaks <- c(seq(0, min (Alignment_per_Rh_G96$Dev_from_horizon_in_degree)-5, by = -5), seq(0, max(Alignment_per_Rh_G96$Dev_from_horizon_in_degree)+5, by = 5))
    breaks <-  seq(-95, 95, by = 10)
    Alignment_plot <- get(paste0("Alignment_per_Rh_", ID))
    Alignment_plot$group <- cut(Alignment_plot$Dev_from_horizon_in_degree, breaks = breaks, labels = FALSE)
    
    # im Blakendiagramm müssen immer die selben bins abgefragt und eventuell mit null "befüllt" werden
    
    p2 <-ggplot(data=Alignment_plot, mapping=aes())+
      geom_bar(data = Alignment_plot, aes(x=group), position = "identity")+
      scale_x_continuous(breaks = c(seq(0.5,8.5, by =2),9.5,10.5,seq(11.5,19.5, by=2)), labels = c(seq(-95,-15,by =20),-5,5,seq(15,95, by= 20)),expand = c(0,0))+
      scale_y_continuous(expand = c(0,0))+
      theme_classic()+
      theme(axis.text = element_text(size= 25), axis.title = element_text(size = 50),plot.margin = margin(30, 60, 30, 30))+
      xlab("Deviation from Horizon in °")+
      ylab("Nr. of Rhabdomers")+
      coord_cartesian(xlim = c(0, 20))+ 
      geom_text(size=18, aes(x = max(20), y = max(as.numeric(table(group)[as.character(names(which.max(table(group))))])), label = paste0("mean = ",round(mean(Dev_from_horizon_in_degree), digits = 2),"\n",
                                                                                                                                          "std = ",round(sd(Dev_from_horizon_in_degree), digits = 2),"\n",
                                                                                                                                          "n = ",length(Dev_from_horizon_in_degree)), hjust = 1, vjust = 1))
    #Kontrolliere den Durschnitt einer Gruppe, um sicher zu gehen das der Plot die richtigen Label hat!
    print(Alignment_plot %>%
            filter(group == 8) %>%
            summarise( mean(Dev_from_horizon_in_degree))) 
    #speichern
    ggsave(paste0(wd,"/",ID,type,"/",ID,type,"Dev._Horz.png"), width = 12, height = 12, plot = p2, dpi = 150 ) 
    
    #erstmal übersicht der Mean absolute standard deviation
    #new data frame for visualisation
    angle_Summary_plot <- get(paste0("angle_Summary_per_Rh_",ID))
    straight_table <- data.frame(new_column = numeric(length(angle_Summary_plot$Mean)))
    straight_table$mean <- angle_Summary_plot$Mean
    
    # für Balkendiagramm gruppieren
    breaks2 <-  seq(0, 60, by = 2.5)
    straight_table$group <- cut(straight_table$mean, breaks = breaks2, labels = FALSE)
    #remove unnecessary column
    straight_table <- straight_table[, -which(names(straight_table) == "new_column")]
    
    p3 <- ggplot(data=straight_table, mapping=aes())+
      geom_bar(data = straight_table, aes(x=group))+
      scale_x_continuous(breaks = seq(0.5,16.5, by = 2),labels = seq(0,40, by = 5), expand = c(0,0))+
      scale_y_continuous(expand = c(0,0))+
      theme_classic()+
      theme(axis.text = element_text(size= 25), axis.title = element_text(size = 50),plot.margin = margin(30, 60, 30, 30))+
      xlab("Mean segment deviation in °")+
      ylab("Nr. of Rhabdomers")+
      coord_cartesian(xlim = c(0, 15))+ 
      geom_text(size=18, aes(x = max(15), y = max(as.numeric(table(straight_table$group)[as.character(names(which.max(table(straight_table$group))))])), label = paste0("mean = ",round(mean(mean), digits = 2),"\n",
                                                                                                                                                                        "std = ",round(sd(mean), digits = 2),"\n",
                                                                                                                                                                        "n = ",length(mean)), hjust = 1, vjust = 1))
    #Kontrolliere den Durschnitt einer Gruppe, um sicher zu gehen das der Plot die richtigen Label hat!
    print(straight_table %>%
            filter(group == 2) %>%
            summarise( mean(mean))) 
    #speichern
    ggsave(paste0(wd,"/",ID,type,"/",ID,type,"Rhb_Str.png"), width = 12, height = 12, plot = p3, dpi = 150 ) 
    
    # für Balkendiagramm gruppieren 5-step
    breaks3 <-  seq(-60, 60, by = 5)
    RAW_Angles$group <- cut(RAW_Angles$dev, breaks = breaks3, labels = FALSE)
    #plot
    p4 <- ggplot(data= RAW_Angles)+
      geom_bar( data = RAW_Angles, mapping = aes(x=group))+
      scale_x_continuous(breaks = c(seq(0.5,11.5, by = 2),11.5,12.5,13.5,seq(14.5,25.5, by=2)),labels = c(seq(-60,-10, by=10),-5,0,5,seq(10,60, by=10)), expand = c(0,0))+
      scale_y_continuous(expand = c(0,0))+
      theme_classic()+
      theme(axis.text = element_text(size= 25), axis.title = element_text(size = 50),plot.margin = margin(30, 60, 30, 30))+
      xlab("Absolute segment deviation in °")+
      ylab("Nr. of Segments")+
      coord_cartesian(xlim = c(0, 25)) + 
      geom_text(size=18, aes(x = max(24), y = max(as.numeric(table(group)[as.character(names(which.max(table(group))))])), label = paste0("mean = ",round(mean(dev), digits = 2),"\n",
                                                                                                                                          "std = ",round(sd(dev), digits = 2),"\n",
                                                                                                                                          "n = ",length(dev)), hjust = 1, vjust = 1))
    #Kontrolliere den Durschnitt einer Gruppe, um sicher zu gehen das der Plot die richtigen Label hat!
    print(RAW_Angles %>%
            filter(group == 13) %>%
            summarise( mean(dev)))
    #speichern
    ggsave(paste0(wd,"/",ID,type,"/",ID,type,"Abs_Seg_Dev.png"), width = 12, height = 12, plot = p4, dpi = 150 ) 
    
    #Balkendiagramm - länge der Rhabdomere 
    
    temp_Table5 <- get(paste0("angle_Summary_per_Rh_",ID))
    breaks5 <-  seq(0, 20, by = 0.5)
    temp_Table5$group <- cut(temp_Table5$Length, breaks = breaks5, labels = FALSE)
    
    max(temp_Table5$group)
    
    p5 <- ggplot(data = temp_Table5)+
      geom_bar(data = temp_Table5, mapping = aes (x=group))+
      scale_x_continuous(breaks=c(0,seq(2.5,24.5,2)), labels = seq(0,12,1))+
      scale_y_continuous(expand = c(0,0))+
      theme_classic()+
      theme(axis.text = element_text(size= 25), axis.title = element_text(size = 50),plot.margin = margin(30, 60, 30, 30))+
      xlab("Length of Rhabdomers in µm")+
      ylab("Nr. of Rhabdomers")+
      #coord_cartesian(xlim = c(0, 25)) + 
      geom_text(size=18, aes(x = max(24), y = max(as.numeric(table(group)[as.character(names(which.max(table(group))))])), label = paste0("n = ",length(Length)), hjust = 1, vjust = 1))    
    #speichern
    ggsave(paste0(wd,"/",ID,type,"/",ID,type,"Rahbdom length.png"), width = 12, height = 12, plot = p5, dpi = 150 ) 
    
    min(temp_Table5$group)
    #Scatter plot Länge/Geradheit/
    p6 <- ggplot(data = temp_Table5,aes (x=Length,y=Mean))+
      geom_point(data = temp_Table5)+
      geom_smooth(method = "lm", se = FALSE)+
      theme_classic()+
      xlab("Length of Rhabdomers in µm")+
      ylab("Mean deviation of rhabdos in ° - straightnes")+
    theme(axis.text = element_text(size= 25), axis.title = element_text(size = 25),plot.margin = margin(30, 60, 30, 30))
    #speichern
    ggsave(paste0(wd,"/",ID,type,"/",ID,type,"LängeGeradheit.png"), width = 12, height = 12, plot = p6, dpi = 150 ) 
    
    
    
    
    ####House keeping####
    
  #Variablen aufräumen
temp_list <- ls(pattern = "plot")
rm(list = temp_list)
rm(RAW_Angles)
    
    
    ####House keeping####
    
    Zähler <- Zähler+1
 temp_list <- ls(pattern = "Rhabdomer")
 rm(list = temp_list)
 temp_list <- ls(pattern = "overall")
 rm(list = temp_list)
temp_list <- ls(pattern = "R_straightn_")
 rm(list = temp_list)
temp_list <- ls(pattern = "temp")
rm(list = temp_list)
    
    #nicht genutzte Daten sollten nach jedem Durchlauf eines Loops gelöscht werden, 
    #damit sie in weiteren Iterationen nicht genutzt werden wenn der Zweite Datensatz 
    #zu klein ist um alle zu überschreiben
  }
  
  colnames(Str_master) <-  c("DegS","nr","rhlength","ID","type","genus","species","sex","mod")
  filename_m <- paste0(wd,"/save/Str_master.txt")
  write.table(Str_master, file = filename_m, sep = "\t", quote = FALSE, row.names = FALSE)
  #Tada!
}
####Ende####
