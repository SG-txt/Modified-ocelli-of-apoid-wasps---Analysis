####Library and WD####
{
library ("plyr")
library("dplyr")
library("tiff")
library("ggplot2")
library ("viridis")
library(RANN)
library("skimr")
library ("qqplotr")
library ("rstatix")
library("effsize")
  library("tidyr")
setwd("C:/Users/Stefan/Desktop/Dr/Analysen/Pol.Sens. R-Script/PolTest")
wd <- getwd()
genus_order <- c("Ampulex","Sceliphron", "Ammophila","Bembix","Bicyrtes","Microbembex", "Stictiella","Liris","Tachysphex","Tachytes")
}
####------------------------------------------------------------------####
#### load master tables and species list####
Str <- read.table(paste0(wd,"/save/Str_master.txt"), header = TRUE)
Str$DegS <- as.numeric(Str$DegS)
if (dir.exists(paste0(wd,"/save/Str"))){
  print("folder already exists")
} else {
  dir.create(paste0(wd,"/save/Str"))}
Alg <- read.table(paste0(wd,"/save/Alg_master.txt"), header = TRUE)
Alg$DegA <- as.numeric(Alg$DegA)
if (dir.exists(paste0(wd,"/save/Alg"))){
  print("folder already exists")
} else {
  dir.create(paste0(wd,"/save/Alg"))}
genera_list <- unique(Str$genus)
Str_temp <- Str
Alg_temp <- Alg
####------------------------------------------------------------------####
####Rainclouds####
{
### This script creates an R function to generate raincloud plots, then simulates
### data for plots. If using for your own data, you only need lines 1-80.
### It relies largely on code previously written by David Robinson
### (https://gist.github.com/dgrtwo/eb7750e74997891d7c20)
### and the package ggplot2 by Hadley Wickham

  # Check if required packages are installed ----
packages <- c("cowplot", "readr", "ggplot2", "dplyr", "lavaan", "smooth", "Hmisc")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))
}

  # Load packages ----
library(ggplot2)

  # Defining the geom_flat_violin function ----
# Note: the below code modifies the
# existing github page by removing a parenthesis in line 50

"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
            data$width <- data$width %||%
              params$width %||% (resolution(data$x, FALSE) * 0.9)
            
            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            data %>%
              group_by(group) %>%
              mutate(
                ymin = min(y),
                ymax = max(y),
                xmin = x,
                xmax = x + width / 2
              )
          },
          
          draw_group = function(data, panel_scales, coord) {
            # Find the points for the line to go all the way around
            data <- transform(data,
                              xminv = x,
                              xmaxv = x + violinwidth * (xmax - x)
            )
            
            # Make sure it's sorted properly to draw the outline
            newdata <- rbind(
              plyr::arrange(transform(data, x = xminv), y),
              plyr::arrange(transform(data, x = xmaxv), -y)
            )
            
            # Close the polygon: set first and last point the same
            # Needed for coord_polar and such
            newdata <- rbind(newdata, newdata[1, ])
            
            ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
          },
          
          draw_key = draw_key_polygon,
          
          default_aes = aes(
            weight = 1, colour = "grey20", fill = "white", size = 0.5,
            alpha = NA, linetype = "solid"
          ),
          
          required_aes = c("x", "y")
  )
}
####------------------------------------------------------------------####
####Sample Size####
{
  nrN <- data.frame()
  genera_list <- unique(Str$genus)
  for (i in 1:length(genera_list)) {
    
    subset_temp <- filter (Str, genus == genera_list[i] & sex == "f" & type == "LO" )
    if (length(subset_temp$DegS)<1) {
      n1 <- "NA"
    }else{
      n1 <- length(subset_temp$DegS)
    }
    
    subset_temp <- filter (Str, genus == genera_list[i] & sex == "f" & (type == "MO" | type == "hMO"))
    if (length(subset_temp$DegS)<1) {
      n2 <- "NA"
    }else{
      n2 <- length(subset_temp$DegS)
    }
    
    subset_temp <- filter (Str, genus == genera_list[i] & sex == "m" & type == "LO" )
    if (length(subset_temp$DegS)<1) {
      n3 <- "NA"
    }else{
      n3 <- length(subset_temp$DegS)
    }
    
    subset_temp <- filter (Str, genus == genera_list[i] & sex == "m" & (type == "MO" | type == "hMO"))
    if (length(subset_temp$DegS)<1) {
      n4 <- "NA"
    }else{
      n4 <- length(subset_temp$DegS)
    }
    
    row <- c(n1,n2,n3,n4)
    
    nrN <- rbind(nrN,row)
  }
  colnames(nrN) <- c("fLO", "fMO", "mLO", "mMO")
  rownames(nrN) <- genera_list
  filename <- paste0(wd,"/save/SampleSize.txt")
  write.table(nrN, file = filename, sep = "\t", row.names = TRUE)
  
}
####------------------------------------------------------------------####
####Mean values####


df <- filter(Alg, ID == "G116", type == "LO")
mean(df$DegA)

comb <- cbind(Alg,Str)
comb <- comb[ , !duplicated(names(comb))] # remove double columns
comb <- comb %>% select(-c(nr,xmid,ymid,x1,x2,y1,y2)) # drop nr column with smal"n"
comb <- comb[, c("Nr", "ID", "genus","species","sex","type","mod", "DegA","DegS","rhlength", "complete")] #reorder columns


mean_master <- comb %>%
  group_by(genus, sex, type) %>%  # grouping
  summarize(
    mean_DegA = round(mean(DegA, na.rm = TRUE), 1),
    sd_DegA = round(sd(DegA, na.rm = TRUE), 1),
    mean_DegS = round(mean(DegS, na.rm = TRUE), 1),
    sd_DegS = round(sd(DegS, na.rm = TRUE), 1),
    mean_length = round(mean(rhlength, na.rm = TRUE), 1),
    sd_length = round(sd(rhlength, na.rm = TRUE), 1),
    # keep additional information from ID, complete, mod
    ID = first(ID),  # first -> keep first value in group (doesn't really matter, all the same)
    complete = first(complete),
    mod = first(mod)
  ) %>%
  ungroup()



filename <- paste0(wd,"/save/mean_master.txt")
write.table(mean_master, file = filename, sep = "\t", quote = FALSE, row.names = TRUE)


####for publication female####
{pb <- mean_master
pb_f <- filter(pb, sex == "f")
pb_f <- pb_f  %>% select(-sex)
pb_f <- pb_f  %>% select(-ID)
pb_f <- pb_f  %>% select(-complete)

pb_f <- pb_f %>%
  # Zuerst nach Genus und Sex gruppieren, um sicherzustellen, dass die Daten korrekt organisiert werden
  group_by(type) %>%
  # Pivot-Wide-Transformation
  pivot_wider(
    names_from = c(type),  # "type" wird verwendet, um die Spalten zu erstellen
    values_from = c(mean_DegA,sd_DegA,mean_DegS,sd_DegS,mean_length,sd_length,mod),  # "td" und "Mod" werden in den Zellen eingefügt
    names_sep = "_"  # Spaltenname wird kombiniert mit einem Unterstrich, z.B. td_f_type1
  ) %>%
  # Spalte sex beibehalten
  ungroup()

pb_f <- pb_f[,c("genus","mean_DegA_MO","sd_DegA_MO","mean_DegS_MO","sd_DegS_MO","mean_length_MO","sd_length_MO","mod_MO","mean_DegA_LO","sd_DegA_LO","mean_DegS_LO","sd_DegS_LO","mean_length_LO","sd_length_LO","mod_LO")]

genus_order <- c("Ampulex","Sceliphron", "Ammophila","Bembix","Microbembex", "Liris","Tachysphex","Tachytes")
pb_f$genus <- factor(pb_f$genus, levels = genus_order)
pb_f <- pb_f %>%
  arrange(genus)
write.table(pb_f,paste0(wd,"/save/A_f_Mean.txt"),sep = "\t", quote = FALSE, row.names = FALSE )}

####for publication male####
{pb <- mean_master
pb_m <- filter(pb, sex == "m")
pb_m <- pb_m  %>% select(-sex)
pb_m <- pb_m  %>% select(-ID)
pb_m <- pb_m  %>% select(-complete)

pb_m <- pb_m %>%
  # Zuerst nach Genus und Sex gruppieren, um sicherzustellen, dass die Daten korrekt organisiert werden
  group_by(type) %>%
  # Pivot-Wide-Transformation
  pivot_wider(
    names_from = c(type),  # "type" wird verwendet, um die Spalten zu erstellen
    values_from = c(mean_DegA,sd_DegA,mean_DegS,sd_DegS,mean_length,sd_length,mod),  # "td" und "Mod" werden in den Zellen eingefügt
    names_sep = "_"  # Spaltenname wird kombiniert mit einem Unterstrich, z.B. td_f_type1
  ) %>%
  # Spalte sex beibehalten
  ungroup()

pb_m <- pb_m[,c("genus","mean_DegA_MO","sd_DegA_MO","mean_DegS_MO","sd_DegS_MO","mean_length_MO","sd_length_MO","mod_MO","mean_DegA_LO","sd_DegA_LO","mean_DegS_LO","sd_DegS_LO","mean_length_LO","sd_length_LO","mod_LO")]

genus_order <- c("Ampulex","Sceliphron", "Ammophila","Bembix","Bicyrtes","Microbembex", "Stictiella","Liris","Tachysphex")
pb_m$genus <- factor(pb_m$genus, levels = genus_order)
pb_m <- pb_m %>%
  arrange(genus)
write.table(pb_m,paste0(wd,"/save/A_m_Mean.txt"),sep = "\t", quote = FALSE, row.names = FALSE )}

####------------------------------------------------------------------####
####Rhabdom straightnes####
####Differences between sexes?####
####Stats####
{
  result <- data.frame(matrix(ncol = 0, nrow = length(genera_list)))
  temp_df <- data.frame() 
  test_name <- "Str_fm"  
  genera_list <- unique(Str$genus)
  
  for (i in 1:length(genera_list)) {
    subset_temp <- filter (Str, genus == genera_list[i] & type != "cMO") 
    if (length(unique(subset_temp$sex)) ==2 ){
      #Test for normal distribution via qqplotr
      data <- as.data.frame(subset_temp$DegS)
      gg <- ggplot(subset_temp = data, mapping = aes(sample = subset_temp$DegS)) +
        stat_qq_band() +
        stat_qq_line() +
        stat_qq_point() +
        labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
      
      print(gg)
      
      #User Input to confirm normal distribution in qqplot
      question_qqplot <- function(question) {
        answer <- readline(prompt = question)
        while (!(answer %in% c("y", "n"))) {
          cat("Is the data normally distributed according to ggplot.\n answer: (y/n) ")
          answer <- readline(prompt = question)
        }
        return(answer)
      }
      
      answer <- question_qqplot("Is the data normally distributed according to ggplot.\n answer: (y/n) ")
      
      if (answer == "y") {
        cat("Continue with population varinace test and corresponding t-test\n")
        
        var_result <- var.test(DegS ~ sex, subset_temp, 
                               alternative = "two.sided")
        var_result <- var_result$p.value
        
        if (var_result>0.05) {
          cat("no significant difference in varianz -> >Student's t-test")
          
          a <-filter(subset_temp, type == "MO")
          b <-filter(subset_temp, type == "LO")
          
          summary_temp <- t.test(a$DegS,b$DegS, var.equal = TRUE)
          #pv <- paste0 (summary_temp$p.value,"s")
          #temp_df <- rbind(temp_df,pv)
          if (summary_temp$p.value<0.05) {
            
            effsize_temp <- cohens_d(subset_temp, DegS ~ genus, var.equal = TRUE)
            
            pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          } else {
            pv <- paste0 ("INSIG., " , " s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          }
          
        }else{
          cat("significant difference in varianz -> >Welch's t-test")
          
          a <-filter(subset_temp, type == "MO")
          b <-filter(subset_temp, type == "LO")
          summary_temp <- t.test(a$DegS,b$DegS, var.equal = FALSE)
          #pv <- paste0(summary_temp$p.value,"w")
          #temp_df <- rbind(temp_df,pv)
          
          if (summary_temp$p.value<0.05) {
            
            effsize_temp <- cohens_d(subset_temp, DegS ~ genus, var.equal = FALSE)
            
            pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", w-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          } else {
            pv <- paste0 ("INSIG., " , " w-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          }
          
        }
        
        
        
      } else {
        cat("Continue with Mann-Whitney-U test\n")
        
        
        summary_temp <- wilcox.test(DegS~sex, data = subset_temp, 
                                    exact = FALSE, 
                                    correct = FALSE, 
                                    conf.int = FALSE)
        W_value <- summary_temp$statistic
        #U-value 
        #subsets
        a <- filter(subset_temp, sex == "f")
        b <- filter(subset_temp, sex == "m")
        #ranked value
        a <- a$DegS
        b <- b$DegS
        #ranks
        ranks <- rank(c(a,b))
        #nr of ranks
        a <- length(a)
        b <- length(b)
        #rank sum
        R1 <- sum(ranks[1:a])
        R2 <- sum(ranks[(a+1):(a+b)])
        # U-Values formula
        U1 <- a * b + (a * (a + 1)) / 2 - R1
        U2 <- a * b + (b * (b + 1)) / 2 - R2
        #choose the smaller U value
        U_value <- min(U1, U2)
        
        
        if (summary_temp$p.value<0.05) {
          effsize_temp <- wilcox_effsize(DegS~sex, data = subset_temp, 
                                         exact = FALSE, 
                                         correct = FALSE, 
                                         conf.int = FALSE)
           pv <- paste0(toupper(effsize_temp$magnitude),", r = " , sub("^0+", "",sprintf("%.2f", effsize_temp$effsize)), ", U = " ,U_value,  ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
          temp_df <- rbind(temp_df,pv)
        } else {
          pv <-  paste0("INSIG., U = ",U_value, ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
          temp_df <- rbind(temp_df,pv)
        }
        
        
      }
    }else{
      pv <- "NA"
      temp_df <- rbind(temp_df,pv)
      
    }
    
  }
  rownames(temp_df) <- genera_list
  colnames(temp_df) <- test_name
  assign(test_name,temp_df)
  
  filename <- paste0(wd,"/save/Str/",test_name,".txt")
  write.table(get(test_name), file = filename, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
}
#MO only
{
  result <- data.frame(matrix(ncol = 0, nrow = length(genera_list)))
  temp_df <- data.frame() 
  test_name <- "Str_fm_MO"  
  genera_list <- unique(Str$genus)
  
  for (i in 1:length(genera_list)) {
    subset_temp <- filter (Str, genus == genera_list[i] & (type != "cMO" | type == "MO")) 
    if (length(unique(subset_temp$sex)) ==2 ){
      #Test for normal diStribution via qqplotr
      data <- as.data.frame(subset_temp$DegS)
      gg <- ggplot(subset_temp = data, mapping = aes(sample = subset_temp$DegS)) +
        stat_qq_band() +
        stat_qq_line() +
        stat_qq_point() +
        labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
      
      print(gg)
      
      #User Input to confirm normal diStribution in qqplot
      question_qqplot <- function(question) {
        answer <- readline(prompt = question)
        while (!(answer %in% c("y", "n"))) {
          cat("Is the data normally diStributed according to ggplot.\n answer: (y/n) ")
          answer <- readline(prompt = question)
        }
        return(answer)
      }
      
      answer <- question_qqplot("Is the data normally diStributed according to ggplot.\n answer: (y/n) ")
      
      if (answer == "y") {
        cat("Continue with population varinace test and corresponding t-test\n")
        
        var_result <- var.test(DegS ~ sex, subset_temp, 
                               alternative = "two.sided")
        var_result <- var_result$p.value
        
        if (var_result>0.05) {
          cat("no significant difference in varianz -> >Student's t-test")
          
          a <- filter(subset_temp, sex == "f")
          b <- filter(subset_temp, sex == "m")
          
          summary_temp <- t.test(a$DegS,b$DegS, var.equal = TRUE)
          #pv <- paste0 (summary_temp$p.value,"s")
          #temp_df <- rbind(temp_df,pv)
          if (summary_temp$p.value<0.05) {
            
            effsize_temp <- cohens_d(subset_temp, DegS ~ genus, var.equal = TRUE)
            
            pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          } else {
            pv <- paste0 ("INSIG., " , " s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          }
          
        }else{
          cat("significant difference in varianz -> >Welch's t-test")
          
          a <- filter(subset_temp, sex == "f")
          b <- filter(subset_temp, sex == "m")
          summary_temp <- t.test(a$DegS,b$DegS, var.equal = FALSE)
          #pv <- paste0(summary_temp$p.value,"w")
          #temp_df <- rbind(temp_df,pv)
          
          if (summary_temp$p.value<0.05) {
            
            effsize_temp <- cohens_d(subset_temp, DegS ~ genus, var.equal = FALSE)
            
            pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", w-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          } else {
            pv <- paste0 ("INSIG., " , " w-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          }
          
        }
        
        
      } else {
        cat("Continue with Mann-Whitney-U test\n")
        
        
        summary_temp <- wilcox.test(DegS~sex, data = subset_temp, 
                                    exact = FALSE, 
                                    correct = FALSE, 
                                    conf.int = FALSE)
        W_value <- summary_temp$statistic
        #U-value 
        #subsets
        a <- filter(subset_temp, sex == "f")
        b <- filter(subset_temp, sex == "m")
        #ranked value
        a <- a$DegS
        b <- b$DegS
        #ranks
        ranks <- rank(c(a,b))
        #nr of ranks
        a <- length(a)
        b <- length(b)
        #rank sum
        R1 <- sum(ranks[1:a])
        R2 <- sum(ranks[(a+1):(a+b)])
        # U-Values formula
        U1 <- a * b + (a * (a + 1)) / 2 - R1
        U2 <- a * b + (b * (b + 1)) / 2 - R2
        #choose the smaller U value
        U_value <- min(U1, U2)
        
        
        if (summary_temp$p.value<0.05) {
          effsize_temp <- wilcox_effsize(DegS~sex, data = subset_temp, 
                                         exact = FALSE, 
                                         correct = FALSE, 
                                         conf.int = FALSE)
           pv <- paste0(toupper(effsize_temp$magnitude),", r = " , sub("^0+", "",sprintf("%.2f", effsize_temp$effsize)), ", U = " ,U_value,  ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
          temp_df <- rbind(temp_df,pv)
        } else {
          pv <-  paste0("INSIG., U = ",U_value, ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
          temp_df <- rbind(temp_df,pv)
        }
        
        
      }
    }else{
      pv <- "NA"
      temp_df <- rbind(temp_df,pv)
      
    }
    
  }
  rownames(temp_df) <- genera_list
  colnames(temp_df) <- test_name
  assign(test_name,temp_df)
  
  filename <- paste0(wd,"/save/Str/",test_name,".txt")
  write.table(get(test_name), file = filename, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
}
#LO only
{
  result <- data.frame(matrix(ncol = 0, nrow = length(genera_list)))
  temp_df <- data.frame() 
  test_name <- "Str_fm_LO"  
  genera_list <- unique(Str$genus)
  
  for (i in 1:length(genera_list)) {
    subset_temp <- filter (Str, genus == genera_list[i] & ( type == "LO")) 
    if (length(unique(subset_temp$sex)) ==2 ){
      #Test for normal diStribution via qqplotr
      data <- as.data.frame(subset_temp$DegS)
      gg <- ggplot(subset_temp = data, mapping = aes(sample = subset_temp$DegS)) +
        stat_qq_band() +
        stat_qq_line() +
        stat_qq_point() +
        labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
      
      print(gg)
      
      #User Input to confirm normal diStribution in qqplot
      question_qqplot <- function(question) {
        answer <- readline(prompt = question)
        while (!(answer %in% c("y", "n"))) {
          cat("Is the data normally diStributed according to ggplot.\n answer: (y/n) ")
          answer <- readline(prompt = question)
        }
        return(answer)
      }
      
      answer <- question_qqplot("Is the data normally diStributed according to ggplot.\n answer: (y/n) ")
      
      if (answer == "y") {
        cat("Continue with population varinace test and corresponding t-test\n")
        
        var_result <- var.test(DegS ~ sex, subset_temp, 
                               alternative = "two.sided")
        var_result <- var_result$p.value
        
        if (var_result>0.05) {
          cat("no significant difference in varianz -> >Student's t-test")
          
          a <- filter(subset_temp, sex == "f")
          b <- filter(subset_temp, sex == "m")
          
          summary_temp <- t.test(a$DegS,b$DegS, var.equal = TRUE)
          #pv <- paste0 (summary_temp$p.value,"s")
          #temp_df <- rbind(temp_df,pv)
          if (summary_temp$p.value<0.05) {
            
            effsize_temp <- cohens_d(subset_temp, DegS ~ genus, var.equal = TRUE)
            
            pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          } else {
            pv <- paste0 ("INSIG., " , " s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          }
          
        }else{
          cat("significant difference in varianz -> >Welch's t-test")
          
          a <- filter(subset_temp, sex == "f")
          b <- filter(subset_temp, sex == "m")
          summary_temp <- t.test(a$DegS,b$DegS, var.equal = FALSE)
          #pv <- paste0(summary_temp$p.value,"w")
          #temp_df <- rbind(temp_df,pv)
          
          if (summary_temp$p.value<0.05) {
            
            effsize_temp <- cohens_d(subset_temp, DegS ~ genus, var.equal = FALSE)
            
            pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", w-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          } else {
            pv <- paste0 ("INSIG., " , " w-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          }
          
        }
        
        
      } else {
        cat("Continue with Mann-Whitney-U test\n")
        
        
        summary_temp <- wilcox.test(DegS~sex, data = subset_temp, 
                                    exact = FALSE, 
                                    correct = FALSE, 
                                    conf.int = FALSE)
        W_value <- summary_temp$statistic
        #U-value 
        #subsets
        a <- filter(subset_temp, sex == "f")
        b <- filter(subset_temp, sex == "m")
        #ranked value
        a <- a$DegS
        b <- b$DegS
        #ranks
        ranks <- rank(c(a,b))
        #nr of ranks
        a <- length(a)
        b <- length(b)
        #rank sum
        R1 <- sum(ranks[1:a])
        R2 <- sum(ranks[(a+1):(a+b)])
        # U-Values formula
        U1 <- a * b + (a * (a + 1)) / 2 - R1
        U2 <- a * b + (b * (b + 1)) / 2 - R2
        #choose the smaller U value
        U_value <- min(U1, U2)
        
        
        if (summary_temp$p.value<0.05) {
          effsize_temp <- wilcox_effsize(DegS~sex, data = subset_temp, 
                                         exact = FALSE, 
                                         correct = FALSE, 
                                         conf.int = FALSE)
           pv <- paste0(toupper(effsize_temp$magnitude),", r = " , sub("^0+", "",sprintf("%.2f", effsize_temp$effsize)), ", U = " ,U_value,  ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
          temp_df <- rbind(temp_df,pv)
        } else {
          pv <-  paste0("INSIG., U = ",U_value, ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
          temp_df <- rbind(temp_df,pv)
        }
        
        
      }
    }else{
      pv <- "NA"
      temp_df <- rbind(temp_df,pv)
      
    }
    
  }
  rownames(temp_df) <- genera_list
  colnames(temp_df) <- test_name
  assign(test_name,temp_df)
  
  filename <- paste0(wd,"/save/Str/",test_name,".txt")
  write.table(get(test_name), file = filename, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
}
####Differences between LO/MO?####
####Stats####
#Str_MO/LO_f
{
  result <- data.frame(matrix(ncol = 0, nrow = length(genera_list)))
  temp_df <- data.frame() 
  test_name <- "Str_MOLO_f"  
  genera_list <- unique(Str$genus)
  
  for (i in 1:length(genera_list)) {
    subset_temp <- filter (Str, genus == genera_list[i] & sex == "f" & type != "cMO") 
    if (length(unique(subset_temp$type)) == 2 ){
      #Test for normal distribution via qqplotr
      data <- as.data.frame(subset_temp$DegS)
      gg <- ggplot(subset_temp = data, mapping = aes(sample = subset_temp$DegS)) +
        stat_qq_band() +
        stat_qq_line() +
        stat_qq_point() +
        labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
      
      print(gg)
      
      #User Input to confirm normal distribution in qqplot
      question_qqplot <- function(question) {
        answer <- readline(prompt = question)
        while (!(answer %in% c("y", "n"))) {
          cat("Is the data normally distributed according to ggplot.\n answer: (y/n) ")
          answer <- readline(prompt = question)
        }
        return(answer)
      }
      
      answer <- question_qqplot("Is the data normally distributed according to ggplot.\n answer: (y/n) ")
      
      if (answer == "y") {
        cat("Continue with population varinace test and corresponding t-test\n")
        
        var_result <- var.test(DegS ~ type, subset_temp, 
                               alternative = "two.sided")
        var_result <- var_result$p.value
        
        if (var_result>0.05) {
          cat("no significant difference in varianz -> >Student's t-test")
          
          a <-filter(subset_temp, type == "MO")
          b <-filter(subset_temp, type == "LO")
          
          summary_temp <- t.test(a$DegS,b$DegS, var.equal = TRUE)
          #pv <- paste0 (summary_temp$p.value,"s")
          #temp_df <- rbind(temp_df,pv)
          if (summary_temp$p.value<0.05) {
            
            effsize_temp <- cohens_d(subset_temp, DegS ~ genus, var.equal = TRUE)
            
            pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          } else {
            pv <- paste0 ("INSIG., " , " s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          }
          
        }else{
          cat("significant difference in varianz -> >Welch's t-test")
          
          a <-filter(subset_temp, type == "MO")
          b <-filter(subset_temp, type == "LO")
          summary_temp <- t.test(a$DegS,b$DegS, var.equal = FALSE)
          #pv <- paste0(summary_temp$p.value,"w")
          #temp_df <- rbind(temp_df,pv)
          
          if (summary_temp$p.value<0.05) {
            
            effsize_temp <- cohens_d(subset_temp, DegS ~ genus, var.equal = FALSE)
            
            pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", w-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          } else {
            pv <- paste0 ("INSIG., " , " w-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          }
          
        }
        
        
      } else {
        cat("Continue with Mann-Whitney-U test\n")
        
        
        summary_temp <- wilcox.test(DegS~type, data = subset_temp, 
                                    exact = FALSE, 
                                    correct = FALSE, 
                                    conf.int = FALSE)
        W_value <- summary_temp$statistic
        #U-value 
        #subsets
        a <-filter(subset_temp, type == "MO")
        b <-filter(subset_temp, type == "LO")
        #ranked value
        a <- a$DegS
        b <- b$DegS
        #ranks
        ranks <- rank(c(a,b))
        #nr of ranks
        a <- length(a)
        b <- length(b)
        #rank sum
        R1 <- sum(ranks[1:a])
        R2 <- sum(ranks[(a+1):(a+b)])
        # U-Values formula
        U1 <- a * b + (a * (a + 1)) / 2 - R1
        U2 <- a * b + (b * (b + 1)) / 2 - R2
        #choose the smaller U value
        U_value <- min(U1, U2)
        
        
        if (summary_temp$p.value<0.05) {
          effsize_temp <- wilcox_effsize(DegS~type, data = subset_temp, 
                                         exact = FALSE, 
                                         correct = FALSE, 
                                         conf.int = FALSE)
           pv <- paste0(toupper(effsize_temp$magnitude),", r = " , sub("^0+", "",sprintf("%.2f", effsize_temp$effsize)), ", U = " ,U_value,  ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
          temp_df <- rbind(temp_df,pv)
        } else {
          pv <-  paste0("INSIG., U = ",U_value, ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
          temp_df <- rbind(temp_df,pv)
        }
        
        
      }
    }else{
      pv <- "NA"
      temp_df <- rbind(temp_df,pv)
      
    }
    
  }
  rownames(temp_df) <- genera_list
  colnames(temp_df) <- test_name
  assign(test_name,temp_df)
  
  filename <- paste0(wd,"/save/Str/",test_name,".txt")
  write.table(get(test_name), file = filename, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
}
#Str_MO/LO_m
{
  result <- data.frame(matrix(ncol = 0, nrow = length(genera_list)))
  temp_df <- data.frame() 
  test_name <- "Str_MOLO_m"  
  genera_list <- unique(Str$genus)
  
  for (i in 1:length(genera_list)) {
    subset_temp <- filter (Str, genus == genera_list[i] & sex == "m" & type != "cMO") 
    if (length(unique(subset_temp$type)) ==2 ){
      #Test for normal distribution via qqplotr
      data <- as.data.frame(subset_temp$DegS)
      gg <- ggplot(subset_temp = data, mapping = aes(sample = subset_temp$DegS)) +
        stat_qq_band() +
        stat_qq_line() +
        stat_qq_point() +
        labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
      
      print(gg)
      
      #User Input to confirm normal distribution in qqplot
      question_qqplot <- function(question) {
        answer <- readline(prompt = question)
        while (!(answer %in% c("y", "n"))) {
          cat("Is the data normally distributed according to ggplot.\n answer: (y/n) ")
          answer <- readline(prompt = question)
        }
        return(answer)
      }
      
      answer <- question_qqplot("Is the data normally distributed according to ggplot.\n answer: (y/n) ")
      
      if (answer == "y") {
        cat("Continue with population varinace test and corresponding t-test\n")
        
        var_result <- var.test(DegS ~ type, subset_temp, 
                               alternative = "two.sided")
        var_result <- var_result$p.value
        
        if (var_result>0.05) {
          cat("no significant difference in varianz -> >Student's t-test")
          
          a <-filter(subset_temp, type == "MO")
          b <-filter(subset_temp, type == "LO")
          
          summary_temp <- t.test(a$DegS,b$DegS, var.equal = TRUE)
          #pv <- paste0 (summary_temp$p.value,"s")
          #temp_df <- rbind(temp_df,pv)
          if (summary_temp$p.value<0.05) {
            
            effsize_temp <- cohens_d(subset_temp, DegS ~ genus, var.equal = TRUE)
            
            pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          } else {
            pv <- paste0 ("INSIG., " , " s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          }
          
        }else{
          cat("significant difference in varianz -> >Welch's t-test")
          
          a <-filter(subset_temp, type == "MO")
          b <-filter(subset_temp, type == "LO")
          summary_temp <- t.test(a$DegS,b$DegS, var.equal = FALSE)
          #pv <- paste0(summary_temp$p.value,"w")
          #temp_df <- rbind(temp_df,pv)
          
          if (summary_temp$p.value<0.05) {
            
            effsize_temp <- cohens_d(subset_temp, DegS ~ genus, var.equal = FALSE)
            
            pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", w-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          } else {
            pv <- paste0 ("INSIG., " , " w-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          }
          
        }
        
        
      } else {
        cat("Continue with Mann-Whitney-U test\n")
        
        
        summary_temp <- wilcox.test(DegS~type, data = subset_temp, 
                                    exact = FALSE, 
                                    correct = FALSE, 
                                    conf.int = FALSE)
        W_value <- summary_temp$statistic
        #U-value 
        #subsets
        a <-filter(subset_temp, type == "MO")
        b <-filter(subset_temp, type == "LO")
        #ranked value
        a <- a$DegS
        b <- b$DegS
        #ranks
        ranks <- rank(c(a,b))
        #nr of ranks
        a <- length(a)
        b <- length(b)
        #rank sum
        R1 <- sum(ranks[1:a])
        R2 <- sum(ranks[(a+1):(a+b)])
        # U-Values formula
        U1 <- a * b + (a * (a + 1)) / 2 - R1
        U2 <- a * b + (b * (b + 1)) / 2 - R2
        #choose the smaller U value
        U_value <- min(U1, U2)
        
        
        if (summary_temp$p.value<0.05) {
          effsize_temp <- wilcox_effsize(DegS~type, data = subset_temp, 
                                         exact = FALSE, 
                                         correct = FALSE, 
                                         conf.int = FALSE)
           pv <- paste0(toupper(effsize_temp$magnitude),", r = " , sub("^0+", "",sprintf("%.2f", effsize_temp$effsize)), ", U = " ,U_value,  ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
          temp_df <- rbind(temp_df,pv)
        } else {
          pv <-  paste0("INSIG., U = ",U_value, ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
          temp_df <- rbind(temp_df,pv)
        }
        
        
      }
    }else{
      pv <- "NA"
      temp_df <- rbind(temp_df,pv)
      
    }
    
  }
  rownames(temp_df) <- genera_list
  colnames(temp_df) <- test_name
  assign(test_name,temp_df)
  
  filename <- paste0(wd,"/save/Str/",test_name,".txt")
  write.table(get(test_name), file = filename, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
}
####Differences betweem Mod/Unmod####
####Stats####

#MO males
{
  Str_temp <- Str
  Str_temp$DegS <-as.numeric(Str_temp$DegS)
  genera_list <- unique(Str$genus)
  
  #Name the test
  test_name <- "Str_MO_m"
  #don't forget to adjust subset filter settings! 
  #currently 13 lines down from this one!!!
  
  #start analysis:
  {
    result <- data.frame(matrix(ncol = 0, nrow = length(genera_list)))
    
    #i<- 1
    for (i in 1:length(genera_list)) {
      
      temp_df <- data.frame()
      #temp_df <- colnames(temp_df) <- genera_list[i]
      #j<- 2
      for (j in 1:length(genera_list)) {
        subset_temp <- filter(Str_temp, (genus == genera_list[i] | genus == genera_list[j]) & type == "MO" & sex == "m")
        
        if (length(unique(subset_temp$genus))==2){
          #Test for normal distribution via qqplotr
          data <- as.data.frame(subset_temp$DegS)
          gg <- ggplot(subset_temp = data, mapping = aes(sample = subset_temp$DegS)) +
            stat_qq_band() +
            stat_qq_line() +
            stat_qq_point() +
            labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
          
          print(gg)
          
          #User Input to confirm normal distribution in qqplot
          question_qqplot <- function(question) {
            answer <- readline(prompt = question)
            while (!(answer %in% c("y", "n"))) {
              cat("Is the data normally distributed according to ggplot.\n answer: (y/n) ")
              answer <- readline(prompt = question)
            }
            return(answer)
          }
          
          answer <- question_qqplot("Is the data normally distributed according to ggplot.\n answer: (y/n) ")
          
          if (answer == "y") {
            cat("Continue with population varinace test and corresponding t-test\n")
            
            var_result <- var.test(DegS ~ genus, subset_temp, 
                                   alternative = "two.sided")
            var_result <- var_result$p.value
            
            if (var_result>0.05) {
              cat("no significant difference in varianz -> >Student's t-test")
              
              a <-filter(subset_temp, genus == genera_list[i])
              b <-filter(subset_temp, genus == genera_list[j])
              
              summary_temp <- t.test(a$DegS,b$DegS, var.equal = TRUE)
              #pv <- paste0 (summary_temp$p.value,"s")
              #temp_df <- rbind(temp_df,pv)
              if (summary_temp$p.value<0.05) {
                
                effsize_temp <- cohens_d(subset_temp, DegS ~ genus, var.equal = TRUE)
                
                            pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
                temp_df <- rbind(temp_df,pv)
              } else {
                pv <- paste0 ("INSIG., " , " s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
                temp_df <- rbind(temp_df,pv)
              }
              
            }else{
              cat("significant difference in varianz -> >Welch's t-test")
              
              a <-filter(subset_temp, genus == genera_list[i])
              b <-filter(subset_temp, genus == genera_list[j])
              summary_temp <- t.test(a$DegS,b$DegS, var.equal = FALSE)
              #pv <- paste0(summary_temp$p.value,"w")
              #temp_df <- rbind(temp_df,pv)
              
              if (summary_temp$p.value<0.05) {
                
                effsize_temp <- cohens_d(subset_temp, DegS ~ genus, var.equal = FALSE)
                
               pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", w-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
                temp_df <- rbind(temp_df,pv)
              } else {
                pv <- paste0 ("INSIG., " , " w-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value))) 
                temp_df <- rbind(temp_df,pv)
              }
              
            }
            
            
          } else {
            cat("Continue with Mann-Whitney-U test\n")
            
            
            summary_temp <- wilcox.test(DegS~genus, data = subset_temp, 
                                        exact = FALSE, 
                                        correct = FALSE, 
                                        conf.int = FALSE)
            W_value <- summary_temp$statistic
            #U-value 
            #subsets
            a <-filter(subset_temp, genus == genera_list[i])
            b <-filter(subset_temp, genus == genera_list[j])
            #ranked value
            a <- a$DegS
            b <- b$DegS
            #ranks
            ranks <- rank(c(a,b))
            #nr of ranks
            a <- length(a)
            b <- length(b)
            #rank sum
            R1 <- sum(ranks[1:a])
            R2 <- sum(ranks[(a+1):(a+b)])
            # U-Values formula
            U1 <- a * b + (a * (a + 1)) / 2 - R1
            U2 <- a * b + (b * (b + 1)) / 2 - R2
            #choose the smaller U value
            U_value <- min(U1, U2)
            if (summary_temp$p.value<0.05) {
              effsize_temp <- wilcox_effsize(DegS~genus, data = subset_temp, 
                                             exact = FALSE, 
                                             correct = FALSE, 
                                             conf.int = FALSE)
               pv <- paste0(toupper(effsize_temp$magnitude),", r = " , sub("^0+", "",sprintf("%.2f", effsize_temp$effsize)), ", U = " ,U_value,  ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
              temp_df <- rbind(temp_df,pv)
            } else {
              pv <-  paste0("INSIG., U = ",U_value, ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
              temp_df <- rbind(temp_df,pv)
            }
            
            
          }
        }else{
          pv <- "NA"
          temp_df <- rbind(temp_df,pv)
        }
      }
      result <- cbind(result,temp_df)
    }
    rownames(result) <- genera_list
    colnames(result) <- genera_list
    
  }
  
  assign(test_name,result)
  filename <- paste0(wd,"/save/Str/",test_name,".txt")
  write.table(get(test_name), file = filename, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
}
#MO females 
{
  Str_temp <- (Str)
  Str_temp$DegS <-as.numeric(Str_temp$DegS)
  genera_list <- unique(Str$genus)
  
  #Name the test
  test_name <- "Str_MO_f"
  #don't forget to adjust subset filter settings! 
  #currently 13 lines down from this one!!!
  
  #start analysis:
  {
    result <- data.frame(matrix(ncol = 0, nrow = length(genera_list)))
    
    #i<- 1
    for (i in 1:length(genera_list)) {
      
      temp_df <- data.frame()
      #temp_df <- colnames(temp_df) <- genera_list[i]
      #j<- 2
      for (j in 1:length(genera_list)) {
        subset_temp <- filter(Str_temp, (genus == genera_list[i] | genus == genera_list[j]) & (type == "MO") & sex == "f")
        
        if (length(unique(subset_temp$genus))==2){
          #Test for normal distribution via qqplotr
          data <- as.data.frame(subset_temp$DegS)
          gg <- ggplot(subset_temp = data, mapping = aes(sample = subset_temp$DegS)) +
            stat_qq_band() +
            stat_qq_line() +
            stat_qq_point() +
            labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
          
          print(gg)
          
          #User Input to confirm normal distribution in qqplot
          question_qqplot <- function(question) {
            answer <- readline(prompt = question)
            while (!(answer %in% c("y", "n"))) {
              cat("Is the data normally distributed according to ggplot.\n answer: (y/n) ")
              answer <- readline(prompt = question)
            }
            return(answer)
          }
          
          answer <- question_qqplot("Is the data normally distributed according to ggplot.\n answer: (y/n) ")
          
          if (answer == "y") {
            cat("Continue with population varinace test and corresponding t-test\n")
            
            var_result <- var.test(DegS ~ genus, subset_temp, 
                                   alternative = "two.sided")
            var_result <- var_result$p.value
            
            if (var_result>0.05) {
              cat("no significant difference in varianz -> >Student's t-test")
              
              a <-filter(subset_temp, genus == genera_list[i])
              b <-filter(subset_temp, genus == genera_list[j])
              
              summary_temp <- t.test(a$DegS,b$DegS, var.equal = TRUE)
              #pv <- paste0 (summary_temp$p.value,"s")
              #temp_df <- rbind(temp_df,pv)
              if (summary_temp$p.value<0.05) {
                
                effsize_temp <- cohens_d(subset_temp, DegS ~ genus, var.equal = TRUE)
                
                            pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
                temp_df <- rbind(temp_df,pv)
              } else {
                pv <- paste0 ("INSIG., " , " s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
                temp_df <- rbind(temp_df,pv)
              }
              
            }else{
              cat("significant difference in varianz -> >Welch's t-test")
              
              a <-filter(subset_temp, genus == genera_list[i])
              b <-filter(subset_temp, genus == genera_list[j])
              summary_temp <- t.test(a$DegS,b$DegS, var.equal = FALSE)
              #pv <- paste0(summary_temp$p.value,"w")
              #temp_df <- rbind(temp_df,pv)
              
              if (summary_temp$p.value<0.05) {
                
                effsize_temp <- cohens_d(subset_temp, DegS ~ genus, var.equal = FALSE)
                
               pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", w-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
                temp_df <- rbind(temp_df,pv)
              } else {
                pv <- paste0(round(summary_temp$p.value,4), ", INSIG, t: ",round(summary_temp$statistic,2), ", df: ",round(summary_temp$parameter,4), ",w")
              }
              
            }
            
            
          } else {
            cat("Continue with Mann-Whitney-U test\n")
            
            
            summary_temp <- wilcox.test(DegS~genus, data = subset_temp, 
                                        exact = FALSE, 
                                        correct = FALSE, 
                                        conf.int = FALSE)
            W_value <- summary_temp$statistic
            #U-value 
            #subsets
            a <-filter(subset_temp, genus == genera_list[i])
            b <-filter(subset_temp, genus == genera_list[j])
            #ranked value
            a <- a$DegS
            b <- b$DegS
            #ranks
            ranks <- rank(c(a,b))
            #nr of ranks
            a <- length(a)
            b <- length(b)
            #rank sum
            R1 <- sum(ranks[1:a])
            R2 <- sum(ranks[(a+1):(a+b)])
            # U-Values formula
            U1 <- a * b + (a * (a + 1)) / 2 - R1
            U2 <- a * b + (b * (b + 1)) / 2 - R2
            #choose the smaller U value
            U_value <- min(U1, U2)
            
            if (summary_temp$p.value<0.05) {
              effsize_temp <- wilcox_effsize(DegS~genus, data = subset_temp, 
                                             exact = FALSE, 
                                             correct = FALSE, 
                                             conf.int = FALSE)
               pv <- paste0(toupper(effsize_temp$magnitude),", r = " , sub("^0+", "",sprintf("%.2f", effsize_temp$effsize)), ", U = " ,U_value,  ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
              temp_df <- rbind(temp_df,pv)
            } else {
              pv <-  paste0("INSIG., U = ",U_value, ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
            }
            
            
          }
        }else{
          pv <- "NA"
          temp_df <- rbind(temp_df,pv)
        }
      }
      result <- cbind(result,temp_df)
    }
    rownames(result) <- genera_list
    colnames(result) <- genera_list
    
  }
  
  assign(test_name,result)
  filename <- paste0(wd,"/save/Str/",test_name,".txt")
  write.table(get(test_name), file = filename, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
}
#LO males
{
  Str_temp <- (Str)
  Str_temp$DegS <-as.numeric(Str_temp$DegS)
  genera_list <- unique(Str$genus)
  
  #Name the test
  test_name <- "Str_LO_m"
  #don't forget to adjust subset filter settings! 
  #currently 13 lines down from this one!!!
  
  #start analysis:
  {
    result <- data.frame(matrix(ncol = 0, nrow = length(genera_list)))
    
    #i<- 1
    for (i in 1:length(genera_list)) {
      
      temp_df <- data.frame()
      #temp_df <- colnames(temp_df) <- genera_list[i]
      #j<- 2
      for (j in 1:length(genera_list)) {
        subset_temp <- filter(Str_temp, (genus == genera_list[i] | genus == genera_list[j]) & type == "LO" & sex == "m")
        
        if (length(unique(subset_temp$genus))==2){
          #Test for normal distribution via qqplotr
          data <- as.data.frame(subset_temp$DegS)
          gg <- ggplot(subset_temp = data, mapping = aes(sample = subset_temp$DegS)) +
            stat_qq_band() +
            stat_qq_line() +
            stat_qq_point() +
            labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
          
          print(gg)
          
          #User Input to confirm normal distribution in qqplot
          question_qqplot <- function(question) {
            answer <- readline(prompt = question)
            while (!(answer %in% c("y", "n"))) {
              cat("Is the data normally distributed according to ggplot.\n answer: (y/n) ")
              answer <- readline(prompt = question)
            }
            return(answer)
          }
          
          answer <- question_qqplot("Is the data normally distributed according to ggplot.\n answer: (y/n) ")
          
          if (answer == "y") {
            cat("Continue with population varinace test and corresponding t-test\n")
            
            var_result <- var.test(DegS ~ genus, subset_temp, 
                                   alternative = "two.sided")
            var_result <- var_result$p.value
            
            if (var_result>0.05) {
              cat("no significant difference in varianz -> >Student's t-test")
              
              a <-filter(subset_temp, genus == genera_list[i])
              b <-filter(subset_temp, genus == genera_list[j])
              
              summary_temp <- t.test(a$DegS,b$DegS, var.equal = TRUE)
              #pv <- paste0 (summary_temp$p.value,"s")
              #temp_df <- rbind(temp_df,pv)
              if (summary_temp$p.value<0.05) {
                
                effsize_temp <- cohens_d(subset_temp, DegS ~ genus, var.equal = TRUE)
                
                            pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
                temp_df <- rbind(temp_df,pv)
              } else {
                pv <- paste0 ("INSIG., " , " s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
                temp_df <- rbind(temp_df,pv)
              }
              
            }else{
              cat("significant difference in varianz -> >Welch's t-test")
              
              a <-filter(subset_temp, genus == genera_list[i])
              b <-filter(subset_temp, genus == genera_list[j])
              summary_temp <- t.test(a$DegS,b$DegS, var.equal = FALSE)
              #pv <- paste0(summary_temp$p.value,"w")
              #temp_df <- rbind(temp_df,pv)
              
              if (summary_temp$p.value<0.05) {
                
                effsize_temp <- cohens_d(subset_temp, DegS ~ genus, var.equal = FALSE)
                
               pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", w-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
                temp_df <- rbind(temp_df,pv)
              } else {
                pv <- paste0 ("INSIG., " , " w-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value))) 
                temp_df <- rbind(temp_df,pv)
              }
              
            }
            
            
          } else {
            cat("Continue with Mann-Whitney-U test\n")
            
            
            summary_temp <- wilcox.test(DegS~genus, data = subset_temp, 
                                        exact = FALSE, 
                                        correct = FALSE, 
                                        conf.int = FALSE)
            W_value <- summary_temp$statistic
            #U-value 
            #subsets
            a <-filter(subset_temp, genus == genera_list[i])
            b <-filter(subset_temp, genus == genera_list[j])
            #ranked value
            a <- a$DegS
            b <- b$DegS
            #ranks
            ranks <- rank(c(a,b))
            #nr of ranks
            a <- length(a)
            b <- length(b)
            #rank sum
            R1 <- sum(ranks[1:a])
            R2 <- sum(ranks[(a+1):(a+b)])
            # U-Values formula
            U1 <- a * b + (a * (a + 1)) / 2 - R1
            U2 <- a * b + (b * (b + 1)) / 2 - R2
            #choose the smaller U value
            U_value <- min(U1, U2)
            if (summary_temp$p.value<0.05) {
              effsize_temp <- wilcox_effsize(DegS~genus, data = subset_temp, 
                                             exact = FALSE, 
                                             correct = FALSE, 
                                             conf.int = FALSE)
               pv <- paste0(toupper(effsize_temp$magnitude),", r = " , sub("^0+", "",sprintf("%.2f", effsize_temp$effsize)), ", U = " ,U_value,  ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
              temp_df <- rbind(temp_df,pv)
            } else {
              pv <-  paste0("INSIG., U = ",U_value, ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
             temp_df <- rbind(temp_df,pv)
            }
            
            
          }
        }else{
          pv <- "NA"
          temp_df <- rbind(temp_df,pv)
        }
      }
      result <- cbind(result,temp_df)
    }
    rownames(result) <- genera_list
    colnames(result) <- genera_list
    
  }
  
  assign(test_name,result)
  filename <- paste0(wd,"/save/Str/",test_name,".txt")
  write.table(get(test_name), file = filename, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
}
#LO females
{
  Str_temp <- (Str)
  Str_temp$DegS <-as.numeric(Str_temp$DegS)
  genera_list <- unique(Str$genus)
  
  #Name the test
  test_name <- "Str_LO_f"
  #don't forget to adjust subset filter settings! 
  #currently 13 lines down from this one!!!
  
  #start analysis:
  {
    result <- data.frame(matrix(ncol = 0, nrow = length(genera_list)))
    
    #i<- 1
    for (i in 1:length(genera_list)) {
      
      temp_df <- data.frame()
      #temp_df <- colnames(temp_df) <- genera_list[i]
      #j<- 2
      for (j in 1:length(genera_list)) {
        subset_temp <- filter(Str_temp, (genus == genera_list[i] | genus == genera_list[j]) & type == "LO" & sex == "f")
        
        if (length(unique(subset_temp$genus))==2){
          #Test for normal diStribution via qqplotr
          data <- as.data.frame(subset_temp$DegS)
          gg <- ggplot(subset_temp = data, mapping = aes(sample = subset_temp$DegS)) +
            stat_qq_band() +
            stat_qq_line() +
            stat_qq_point() +
            labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
          
          print(gg)
          
          #User Input to confirm normal diStribution in qqplot
          question_qqplot <- function(question) {
            answer <- readline(prompt = question)
            while (!(answer %in% c("y", "n"))) {
              cat("Is the data normally diStributed according to ggplot.\n answer: (y/n) ")
              answer <- readline(prompt = question)
            }
            return(answer)
          }
          
          answer <- question_qqplot("Is the data normally diStributed according to ggplot.\n answer: (y/n) ")
          
          if (answer == "y") {
            cat("Continue with population varinace test and corresponding t-test\n")
            
            var_result <- var.test(DegS ~ genus, subset_temp, 
                                   alternative = "two.sided")
            var_result <- var_result$p.value
            
            if (var_result>0.05) {
              cat("no significant difference in varianz -> >Student's t-test")
              
              a <-filter(subset_temp, genus == genera_list[i])
              b <-filter(subset_temp, genus == genera_list[j])
              
              summary_temp <- t.test(a$DegS,b$DegS, var.equal = TRUE)
              #pv <- paste0 (summary_temp$p.value,"s")
              #temp_df <- rbind(temp_df,pv)
              if (summary_temp$p.value<0.05) {
                
                effsize_temp <- cohens_d(subset_temp, DegS ~ genus, var.equal = TRUE)
                
                            pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
                temp_df <- rbind(temp_df,pv)
              } else {
                pv <- paste0 ("INSIG., " , " s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
                temp_df <- rbind(temp_df,pv)
              }
              
            }else{
              cat("significant difference in varianz -> >Welch's t-test")
              
              a <-filter(subset_temp, genus == genera_list[i])
              b <-filter(subset_temp, genus == genera_list[j])
              summary_temp <- t.test(a$DegS,b$DegS, var.equal = FALSE)
              #pv <- paste0(summary_temp$p.value,"w")
              #temp_df <- rbind(temp_df,pv)
              
              if (summary_temp$p.value<0.05) {
                
                effsize_temp <- cohens_d(subset_temp, DegS ~ genus, var.equal = FALSE)
                
               pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", w-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
                temp_df <- rbind(temp_df,pv)
              } else {
                pv <- paste0(round(summary_temp$p.value,4), ", INSIG, t: ",round(summary_temp$statistic,2), ", df: ",round(summary_temp$parameter,4), ",w")
                temp_df <- rbind(temp_df,pv)
              }
              
            }
            
            
          } else {
            cat("Continue with Mann-Whitney-U test\n")
            
            
            summary_temp <- wilcox.test(DegS~genus, data = subset_temp, 
                                        exact = FALSE, 
                                        correct = FALSE, 
                                        conf.int = FALSE)
            W_value <- summary_temp$statistic
            #U-value 
            #subsets
            a <-filter(subset_temp, genus == genera_list[i])
            b <-filter(subset_temp, genus == genera_list[j])
            #ranked value
            a <- a$DegS
            b <- b$DegS
            #ranks
            ranks <- rank(c(a,b))
            #nr of ranks
            a <- length(a)
            b <- length(b)
            #rank sum
            R1 <- sum(ranks[1:a])
            R2 <- sum(ranks[(a+1):(a+b)])
            # U-Values formula
            U1 <- a * b + (a * (a + 1)) / 2 - R1
            U2 <- a * b + (b * (b + 1)) / 2 - R2
            #choose the smaller U value
            U_value <- min(U1, U2)
            if (summary_temp$p.value<0.05) {
              effsize_temp <- wilcox_effsize(DegS~genus, data = subset_temp, 
                                             exact = FALSE, 
                                             correct = FALSE, 
                                             conf.int = FALSE)
               pv <- paste0(toupper(effsize_temp$magnitude),", r = " , sub("^0+", "",sprintf("%.2f", effsize_temp$effsize)), ", U = " ,U_value,  ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
              temp_df <- rbind(temp_df,pv)
            } else {
              pv <-  paste0("INSIG., U = ",U_value, ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
              temp_df <- rbind(temp_df,pv)
            }
            
            
          }
        }else{
          pv <- "NA"
          temp_df <- rbind(temp_df,pv)
        }
      }
      result <- cbind(result,temp_df)
    }
    rownames(result) <- genera_list
    colnames(result) <- genera_list
    
  }
  
  assign(test_name,result)
  filename <- paste0(wd,"/save/Str/",test_name,".txt")
  write.table(get(test_name), file = filename, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
}
####Summary by mod/unmod####

#female
subset_temp <-  filter(Str, sex == "f")
subset_n1 <- filter(subset_temp, mod == "y")
print(length(subset_n1$DegS))
print(mean(subset_n1$DegS))
print(sd(subset_n1$DegS))
subset_n2 <- filter(subset_temp, mod == "n")
print(length(subset_n2$DegS))
print(mean(subset_n2$DegS))
print(sd(subset_n2$DegS))

subset_n1 <- filter(subset_temp, mod == "n")
gg <- ggplot(data = subset_temp, mapping = aes(sample = subset_temp$DegS)) +
  stat_qq_band() +
  stat_qq_line() +
  stat_qq_point() +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")

print(gg)
# nicht normal verteilt

summary_temp <- wilcox.test(DegS~mod, data = subset_temp, 
                            exact = FALSE, 
                            correct = FALSE, 
                            conf.int = FALSE)
W_value <- summary_temp$statistic
#U-value 
#subsets
a <-filter(subset_temp, mod == "y")
b <-filter(subset_temp, mod == "n")
#ranked value
a <- a$DegS
b <- b$DegS
#ranks
ranks <- rank(c(a,b))
#nr of ranks
a <- length(a)
b <- length(b)
#rank sum
R1 <- sum(ranks[1:a])
R2 <- sum(ranks[(a+1):(a+b)])
# U-Values formula
U1 <- a * b + (a * (a + 1)) / 2 - R1
U2 <- a * b + (b * (b + 1)) / 2 - R2
#choose the smaller U value
U_value <- min(U1, U2)


effsize_temp <- wilcox_effsize(DegS~mod, data = subset_temp, 
                               exact = FALSE, 
                               correct = FALSE, 
                               conf.int = FALSE)
Subs_mean <- filter(subset_temp, mod =="y")
mean(Subs_mean$DegS)
sd(Subs_mean$DegS)
Subs_mean <- filter(subset_temp, mod =="n")
mean(Subs_mean$DegS)
sd(Subs_mean$DegS)

#male
subset_temp <-  filter(Str, sex == "m")
subset_n1 <- filter(subset_temp, mod == "y")
print(length(subset_n1$DegS))
print(mean(subset_n1$DegS))
print(sd(subset_n1$DegS))
subset_n2 <- filter(subset_temp, mod == "n")
print(length(subset_n2$DegS))
print(mean(subset_n2$DegS))
print(sd(subset_n2$DegS))

gg <- ggplot(data = subset_temp, mapping = aes(sample = subset_temp$DegS)) +
  stat_qq_band() +
  stat_qq_line() +
  stat_qq_point() +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")

print(gg)
# nicht normal verteilt

summary_temp <- wilcox.test(DegS~mod, data = subset_temp, 
                            exact = FALSE, 
                            correct = FALSE, 
                            conf.int = FALSE)
W_value <- summary_temp$statistic
#U-value 
#subsets
a <-filter(subset_temp, mod == "y")
b <-filter(subset_temp, mod == "n")
#ranked value
a <- a$DegS
b <- b$DegS
#ranks
ranks <- rank(c(a,b))
#nr of ranks
a <- length(a)
b <- length(b)
#rank sum
R1 <- sum(ranks[1:a])
R2 <- sum(ranks[(a+1):(a+b)])
# U-Values formula
U1 <- a * b + (a * (a + 1)) / 2 - R1
U2 <- a * b + (b * (b + 1)) / 2 - R2
#choose the smaller U value
U_value <- min(U1, U2)



effsize_temp <- wilcox_effsize(DegS~mod, data = subset_temp, 
                               exact = FALSE, 
                               correct = FALSE, 
                               conf.int = FALSE)
Subs_mean <- filter(subset_temp, mod =="y")
mean(Subs_mean$DegS)
sd(Subs_mean$DegS)
Subs_mean <- filter(subset_temp, mod =="n")
mean(Subs_mean$DegS)
sd(Subs_mean$DegS)
####------------------------------------------------------------------####
####Rhabdom alignment####

####Differences between sexes?####
####Stats####
{
  result <- data.frame(matrix(ncol = 0, nrow = length(genera_list)))
  temp_df <- data.frame() 
  test_name <- "Alg_fm"  
  genera_list <- unique(Alg$genus)
  i<- 5
  for (i in 1:length(genera_list)) {
    subset_temp <- filter (Alg, genus == genera_list[i] & type != "cMO") 
    if (length(unique(subset_temp$sex)) ==2 ){
      #Test for normal diAlgibution via qqplotr
      data <- as.data.frame(subset_temp$DegA)
      gg <- ggplot(subset_temp = data, mapping = aes(sample = subset_temp$DegA)) +
        stat_qq_band() +
        stat_qq_line() +
        stat_qq_point() +
        labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
      
      print(gg)
      
      #User Input to confirm normal diAlgibution in qqplot
      question_qqplot <- function(question) {
        answer <- readline(prompt = question)
        while (!(answer %in% c("y", "n"))) {
          cat("Is the data normally diAlgibuted according to ggplot.\n answer: (y/n) ")
          answer <- readline(prompt = question)
        }
        return(answer)
      }
      
      answer <- question_qqplot("Is the data normally diAlgibuted according to ggplot.\n answer: (y/n) ")
      
      if (answer == "y") {
        cat("Continue with population varinace test and corresponding t-test\n")
        
        var_result <- var.test(DegA ~ sex, subset_temp, 
                               alternative = "two.sided")
        var_result <- var_result$p.value
        
        if (var_result>0.05) {
          cat("no significant difference in varianz -> >Student's t-test")
          
          a <-filter(subset_temp, type == "MO")
          b <-filter(subset_temp, type == "LO")
          
          summary_temp <- t.test(a$DegS,b$DegS, var.equal = TRUE)
          #pv <- paste0 (summary_temp$p.value,"s")
          #temp_df <- rbind(temp_df,pv)
          if (summary_temp$p.value<0.05) {
            
            effsize_temp <- cohens_d(subset_temp, DegS ~ genus, var.equal = TRUE)
            
            pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          } else {
            pv <- paste0 ("INSIG., " , " s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          }
          
        }else{
          cat("significant difference in varianz -> >Welch's t-test")
          
          a <-filter(subset_temp, type == "MO")
          b <-filter(subset_temp, type == "LO")
          summary_temp <- t.test(a$DegS,b$DegS, var.equal = FALSE)
          #pv <- paste0(summary_temp$p.value,"w")
          #temp_df <- rbind(temp_df,pv)
          
          if (summary_temp$p.value<0.05) {
            
            effsize_temp <- cohens_d(subset_temp, DegS ~ genus, var.equal = FALSE)
            
            pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", w-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          } else {
            pv <- paste0 ("INSIG., " , " w-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          }
          
        }
        
        
        
      } else {
        cat("Continue with Mann-Whitney-U test\n")
        
        
        summary_temp <- wilcox.test(DegA~sex, data = subset_temp, 
                                    exact = FALSE, 
                                    correct = FALSE, 
                                    conf.int = FALSE)
        W_value <- summary_temp$statistic
        #U-value 
        #subsets
        a <- filter(subset_temp, sex == "f")
        b <- filter(subset_temp, sex == "m")
        #ranked value
        a <- a$DegA
        b <- b$DegA
        #ranks
        ranks <- rank(c(a,b))
        #nr of ranks
        a <- length(a)
        b <- length(b)
        #rank sum
        R1 <- sum(ranks[1:a])
        R2 <- sum(ranks[(a+1):(a+b)])
        # U-Values formula
        U1 <- a * b + (a * (a + 1)) / 2 - R1
        U2 <- a * b + (b * (b + 1)) / 2 - R2
        #choose the smaller U value
        U_value <- min(U1, U2)
        
        
        if (summary_temp$p.value<0.05) {
          effsize_temp <- wilcox_effsize(DegA~sex, data = subset_temp, 
                                         exact = FALSE, 
                                         correct = FALSE, 
                                         conf.int = FALSE)
           pv <- paste0(toupper(effsize_temp$magnitude),", r = " , sub("^0+", "",sprintf("%.2f", effsize_temp$effsize)), ", U = " ,U_value,  ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
          
          
          temp_df <- rbind(temp_df,pv)
        } else {
          pv <-  paste0("INSIG, U = ",U_value, ", p = ", round(summary_temp$p.value,3)) 
          temp_df <- rbind(temp_df,pv)
        }
        
        
      }
    }else{
      pv <- "NA"
      temp_df <- rbind(temp_df,pv)
      
    }
    
  }
  rownames(temp_df) <- genera_list
  colnames(temp_df) <- test_name
  assign(test_name,temp_df)
  
  filename <- paste0(wd,"/save/Alg/",test_name,".txt")
  write.table(get(test_name), file = filename, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
}
#MO only
{
  result <- data.frame(matrix(ncol = 0, nrow = length(genera_list)))
  temp_df <- data.frame() 
  test_name <- "Alg_fm_MO"  
  genera_list <- unique(Alg$genus)
  
  for (i in 1:length(genera_list)) {
    subset_temp <- filter (Alg, genus == genera_list[i] & (type != "cMO" | type == "MO")) 
    if (length(unique(subset_temp$sex)) ==2 ){
      #Test for normal diAlgibution via qqplotr
      data <- as.data.frame(subset_temp$DegA)
      gg <- ggplot(subset_temp = data, mapping = aes(sample = subset_temp$DegA)) +
        stat_qq_band() +
        stat_qq_line() +
        stat_qq_point() +
        labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
      
      print(gg)
      
      #User Input to confirm normal diAlgibution in qqplot
      question_qqplot <- function(question) {
        answer <- readline(prompt = question)
        while (!(answer %in% c("y", "n"))) {
          cat("Is the data normally diAlgibuted according to ggplot.\n answer: (y/n) ")
          answer <- readline(prompt = question)
        }
        return(answer)
      }
      
      answer <- question_qqplot("Is the data normally diAlgibuted according to ggplot.\n answer: (y/n) ")
      
      if (answer == "y") {
        cat("Continue with population varinace test and corresponding t-test\n")
        
        var_result <- var.test(DegA ~ sex, subset_temp, 
                               alternative = "two.sided")
        var_result <- var_result$p.value
        
        if (var_result>0.05) {
          cat("no significant difference in varianz -> >Student's t-test")
          
          a <- filter(subset_temp, sex == "f")
          b <- filter(subset_temp, sex == "m")
          
          summary_temp <- t.test(a$DegA,b$DegA, var.equal = TRUE)
          #pv <- paste0 (summary_temp$p.value,"s")
          #temp_df <- rbind(temp_df,pv)
          if (summary_temp$p.value<0.05) {
            
            effsize_temp <- cohens_d(subset_temp, DegA ~ genus, var.equal = TRUE)
            
            #pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
            pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          } else {
            #pv <- paste0 ("INSIG., " , " s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
            pv <- paste0 ("INSIG., " , " s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
            
          }
          
        }else{
          cat("significant difference in varianz -> >Welch's t-test")
          
          a <- filter(subset_temp, sex == "f")
          b <- filter(subset_temp, sex == "m")
          summary_temp <- t.test(a$DegA,b$DegA, var.equal = FALSE)
          #pv <- paste0(summary_temp$p.value,"w")
          #temp_df <- rbind(temp_df,pv)
          
          if (summary_temp$p.value<0.05) {
            
            effsize_temp <- cohens_d(subset_temp, DegA ~ genus, var.equal = FALSE)
            
            #pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", w-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
            pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", w-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          } else {
            #pv <- paste0 ("INSIG., " , " w-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
            pv <- paste0 ("INSIG., " , " w-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          }
          
        }
        
        
      } else {
        cat("Continue with Mann-Whitney-U test\n")
        
        
        summary_temp <- wilcox.test(DegA~sex, data = subset_temp, 
                                    exact = FALSE, 
                                    correct = FALSE, 
                                    conf.int = FALSE)
        W_value <- summary_temp$statistic
        #U-value 
        #subsets
        a <- filter(subset_temp, sex == "f")
        b <- filter(subset_temp, sex == "m")
        #ranked value
        a <- a$DegA
        b <- b$DegA
        #ranks
        ranks <- rank(c(a,b))
        #nr of ranks
        a <- length(a)
        b <- length(b)
        #rank sum
        R1 <- sum(ranks[1:a])
        R2 <- sum(ranks[(a+1):(a+b)])
        # U-Values formula
        U1 <- a * b + (a * (a + 1)) / 2 - R1
        U2 <- a * b + (b * (b + 1)) / 2 - R2
        #choose the smaller U value
        U_value <- min(U1, U2)
        
        
        if (summary_temp$p.value<0.05) {
          effsize_temp <- wilcox_effsize(DegA~sex, data = subset_temp, 
                                         exact = FALSE, 
                                         correct = FALSE, 
                                         conf.int = FALSE)
          # pv <- paste0(toupper(effsize_temp$magnitude),", r = " , sub("^0+", "",sprintf("%.2f", effsize_temp$effsize)), ", U = " ,U_value,  ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
          pv <- paste0(toupper(effsize_temp$magnitude),", r = " , sub("^0+", "",sprintf("%.2f", effsize_temp$effsize)), ", U = " ,U_value,  ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
          temp_df <- rbind(temp_df,pv)
        } else {
          #pv <-  paste0("INSIG., U = ",U_value, ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
          pv <-  paste0("INSIG., U = ",U_value, ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
          temp_df <- rbind(temp_df,pv)
        }
        
        
      }
    }else{
      pv <- "NA"
      temp_df <- rbind(temp_df,pv)
      
    }
    
  }
  rownames(temp_df) <- genera_list
  colnames(temp_df) <- test_name
  assign(test_name,temp_df)
  
  
 
  filename <- paste0(wd,"/save/Alg/",test_name,".txt")
  write.table(get(test_name), file = filename, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
}
#LO only
{
  result <- data.frame(matrix(ncol = 0, nrow = length(genera_list)))
  temp_df <- data.frame() 
  test_name <- "Alg_fm_LO"  
  genera_list <- unique(Alg$genus)
  
  for (i in 1:length(genera_list)) {
    subset_temp <- filter (Alg, genus == genera_list[i] & ( type == "LO")) 
    if (length(unique(subset_temp$sex)) ==2 ){
      #Test for normal diAlgibution via qqplotr
      data <- as.data.frame(subset_temp$DegA)
      gg <- ggplot(subset_temp = data, mapping = aes(sample = subset_temp$DegA)) +
        stat_qq_band() +
        stat_qq_line() +
        stat_qq_point() +
        labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
      
      print(gg)
      
      #User Input to confirm normal diAlgibution in qqplot
      question_qqplot <- function(question) {
        answer <- readline(prompt = question)
        while (!(answer %in% c("y", "n"))) {
          cat("Is the data normally diAlgibuted according to ggplot.\n answer: (y/n) ")
          answer <- readline(prompt = question)
        }
        return(answer)
      }
      
      answer <- question_qqplot("Is the data normally diAlgibuted according to ggplot.\n answer: (y/n) ")
      
      if (answer == "y") {
        cat("Continue with population varinace test and corresponding t-test\n")
        
        var_result <- var.test(DegA ~ sex, subset_temp, 
                               alternative = "two.sided")
        var_result <- var_result$p.value
        
        if (var_result>0.05) {
          cat("no significant difference in varianz -> >Student's t-test")
          
          a <- filter(subset_temp, sex == "f")
          b <- filter(subset_temp, sex == "m")
          
          summary_temp <- t.test(a$DegA,b$DegA, var.equal = TRUE)
          #pv <- paste0 (summary_temp$p.value,"s")
          #temp_df <- rbind(temp_df,pv)
          if (summary_temp$p.value<0.05) {
            
            effsize_temp <- cohens_d(subset_temp, DegA ~ genus, var.equal = TRUE)
            
            pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          } else {
            pv <- paste0 ("INSIG., " , " s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          }
          
        }else{
          cat("significant difference in varianz -> >Welch's t-test")
          
          a <- filter(subset_temp, sex == "f")
          b <- filter(subset_temp, sex == "m")
          summary_temp <- t.test(a$DegA,b$DegA, var.equal = FALSE)
          #pv <- paste0(summary_temp$p.value,"w")
          #temp_df <- rbind(temp_df,pv)
          
          if (summary_temp$p.value<0.05) {
            
            effsize_temp <- cohens_d(subset_temp, DegA ~ genus, var.equal = FALSE)
            
            pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", w-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          } else {
            pv <- paste0 ("INSIG., " , " w-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          }
          
        }
        
        
      } else {
        cat("Continue with Mann-Whitney-U test\n")
        
        
        summary_temp <- wilcox.test(DegA~sex, data = subset_temp, 
                                    exact = FALSE, 
                                    correct = FALSE, 
                                    conf.int = FALSE)
        W_value <- summary_temp$statistic
        #U-value 
        #subsets
        a <- filter(subset_temp, sex == "f")
        b <- filter(subset_temp, sex == "m")
        #ranked value
        a <- a$DegA
        b <- b$DegA
        #ranks
        ranks <- rank(c(a,b))
        #nr of ranks
        a <- length(a)
        b <- length(b)
        #rank sum
        R1 <- sum(ranks[1:a])
        R2 <- sum(ranks[(a+1):(a+b)])
        # U-Values formula
        U1 <- a * b + (a * (a + 1)) / 2 - R1
        U2 <- a * b + (b * (b + 1)) / 2 - R2
        #choose the smaller U value
        U_value <- min(U1, U2)
        
        
        if (summary_temp$p.value<0.05) {
          effsize_temp <- wilcox_effsize(DegA~sex, data = subset_temp, 
                                         exact = FALSE, 
                                         correct = FALSE, 
                                         conf.int = FALSE)
           pv <- paste0(toupper(effsize_temp$magnitude),", r = " , sub("^0+", "",sprintf("%.2f", effsize_temp$effsize)), ", U = " ,U_value,  ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
          temp_df <- rbind(temp_df,pv)
        } else {
          pv <-  paste0("INSIG., U = ",U_value, ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
          temp_df <- rbind(temp_df,pv)
        }
        
        
      }
    }else{
      pv <- "NA"
      temp_df <- rbind(temp_df,pv)
      
    }
    
  }
  rownames(temp_df) <- genera_list
  colnames(temp_df) <- test_name
  assign(test_name,temp_df)
  
  filename <- paste0(wd,"/save/Alg/",test_name,".txt")
  write.table(get(test_name), file = filename, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
}
####Differences between LO/MO?####
####Stats####
#Alg_MO/LO_f
{
  result <- data.frame(matrix(ncol = 0, nrow = length(genera_list)))
  temp_df <- data.frame() 
  test_name <- "Alg_MOLO_f"  
  genera_list <- unique(Alg$genus)
  
  for (i in 1:length(genera_list)) {
    subset_temp <- filter (Alg, genus == genera_list[i] & sex == "f" & type != "cMO") 
    if (length(unique(subset_temp$type)) ==2 ){
      #Test for normal diAlgibution via qqplotr
      data <- as.data.frame(subset_temp$DegA)
      gg <- ggplot(subset_temp = data, mapping = aes(sample = subset_temp$DegA)) +
        stat_qq_band() +
        stat_qq_line() +
        stat_qq_point() +
        labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
      
      print(gg)
      
      #User Input to confirm normal diAlgibution in qqplot
      question_qqplot <- function(question) {
        answer <- readline(prompt = question)
        while (!(answer %in% c("y", "n"))) {
          cat("Is the data normally diAlgibuted according to ggplot.\n answer: (y/n) ")
          answer <- readline(prompt = question)
        }
        return(answer)
      }
      
      answer <- question_qqplot("Is the data normally diAlgibuted according to ggplot.\n answer: (y/n) ")
      
      if (answer == "y") {
        cat("Continue with population varinace test and corresponding t-test\n")
        
        var_result <- var.test(DegA ~ type, subset_temp, 
                               alternative = "two.sided")
        var_result <- var_result$p.value
        
        if (var_result>0.05) {
          cat("no significant difference in varianz -> >Student's t-test")
          
          a <-filter(subset_temp, type == "MO")
          b <-filter(subset_temp, type == "LO")
          
          summary_temp <- t.test(a$DegA,b$DegA, var.equal = TRUE)
          #pv <- paste0 (summary_temp$p.value,"s")
          #temp_df <- rbind(temp_df,pv)
          if (summary_temp$p.value<0.05) {
            
            effsize_temp <- cohens_d(subset_temp, DegA ~ genus, var.equal = TRUE)
            
            pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          } else {
            pv <- paste0 ("INSIG., " , " s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          }
          
        }else{
          cat("significant difference in varianz -> >Welch's t-test")
          
          a <-filter(subset_temp, type == "MO")
          b <-filter(subset_temp, type == "LO")
          summary_temp <- t.test(a$DegA,b$DegA, var.equal = FALSE)
          #pv <- paste0(summary_temp$p.value,"w")
          #temp_df <- rbind(temp_df,pv)
          
          if (summary_temp$p.value<0.05) {
            
            effsize_temp <- cohens_d(subset_temp, DegA ~ genus, var.equal = FALSE)
            
            pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", w-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          } else {
            pv <- paste0 ("INSIG., " , " w-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          }
          
        }
        
        
      } else {
        cat("Continue with Mann-Whitney-U test\n")
        
        
        summary_temp <- wilcox.test(DegA~type, data = subset_temp, 
                                    exact = FALSE, 
                                    correct = FALSE, 
                                    conf.int = FALSE)
        W_value <- summary_temp$statistic
        #U-value 
        #subsets
        a <-filter(subset_temp, type == "MO")
        b <-filter(subset_temp, type == "LO")
        #ranked value
        a <- a$DegA
        b <- b$DegA
        #ranks
        ranks <- rank(c(a,b))
        #nr of ranks
        a <- length(a)
        b <- length(b)
        #rank sum
        R1 <- sum(ranks[1:a])
        R2 <- sum(ranks[(a+1):(a+b)])
        # U-Values formula
        U1 <- a * b + (a * (a + 1)) / 2 - R1
        U2 <- a * b + (b * (b + 1)) / 2 - R2
        #choose the smaller U value
        U_value <- min(U1, U2)
        
        
        if (summary_temp$p.value<0.05) {
          effsize_temp <- wilcox_effsize(DegA~type, data = subset_temp, 
                                         exact = FALSE, 
                                         correct = FALSE, 
                                         conf.int = FALSE)
           pv <- paste0(toupper(effsize_temp$magnitude),", r = " , sub("^0+", "",sprintf("%.2f", effsize_temp$effsize)), ", U = " ,U_value,  ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
          temp_df <- rbind(temp_df,pv)
        } else {
          pv <-  paste0("INSIG., U = ",U_value, ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
          temp_df <- rbind(temp_df,pv)
        }
        
        
      }
    }else{
      pv <- "NA"
      temp_df <- rbind(temp_df,pv)
      
    }
    
  }
  rownames(temp_df) <- genera_list
  colnames(temp_df) <- test_name
  assign(test_name,temp_df)
  
  filename <- paste0(wd,"/save/Alg/",test_name,".txt")
  write.table(get(test_name), file = filename, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
}
#Alg_MO/LO_m
{
  result <- data.frame(matrix(ncol = 0, nrow = length(genera_list)))
  temp_df <- data.frame() 
  test_name <- "Alg_MOLO_m"  
  genera_list <- unique(Alg$genus)
  
  for (i in 1:length(genera_list)) {
    subset_temp <- filter (Alg, genus == genera_list[i] & sex == "m" & type != "cMO") 
    if (length(unique(subset_temp$type)) ==2 ){
      #Test for normal diAlgibution via qqplotr
      data <- as.data.frame(subset_temp$DegA)
      gg <- ggplot(subset_temp = data, mapping = aes(sample = subset_temp$DegA)) +
        stat_qq_band() +
        stat_qq_line() +
        stat_qq_point() +
        labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
      
      print(gg)
      
      #User Input to confirm normal diAlgibution in qqplot
      question_qqplot <- function(question) {
        answer <- readline(prompt = question)
        while (!(answer %in% c("y", "n"))) {
          cat("Is the data normally diAlgibuted according to ggplot.\n answer: (y/n) ")
          answer <- readline(prompt = question)
        }
        return(answer)
      }
      
      answer <- question_qqplot("Is the data normally diAlgibuted according to ggplot.\n answer: (y/n) ")
      
      if (answer == "y") {
        cat("Continue with population varinace test and corresponding t-test\n")
        
        var_result <- var.test(DegA ~ type, subset_temp, 
                               alternative = "two.sided")
        var_result <- var_result$p.value
        
        if (var_result>0.05) {
          cat("no significant difference in varianz -> >Student's t-test")
          
          a <-filter(subset_temp, type == "MO")
          b <-filter(subset_temp, type == "LO")
          
          summary_temp <- t.test(a$DegA,b$DegA, var.equal = TRUE)
          #pv <- paste0 (summary_temp$p.value,"s")
          #temp_df <- rbind(temp_df,pv)
          if (summary_temp$p.value<0.05) {
            
            effsize_temp <- cohens_d(subset_temp, DegA ~ genus, var.equal = TRUE)
            
            pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          } else {
            pv <- paste0 ("INSIG., " , " s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          }
          
        }else{
          cat("significant difference in varianz -> >Welch's t-test")
          
          a <-filter(subset_temp, type == "MO")
          b <-filter(subset_temp, type == "LO")
          summary_temp <- t.test(a$DegA,b$DegA, var.equal = FALSE)
          #pv <- paste0(summary_temp$p.value,"w")
          #temp_df <- rbind(temp_df,pv)
          
          if (summary_temp$p.value<0.05) {
            
            effsize_temp <- cohens_d(subset_temp, DegA ~ genus, var.equal = FALSE)
            
            pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", w-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          } else {
            pv <- paste0 ("INSIG., " , " w-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
            temp_df <- rbind(temp_df,pv)
          }
          
        }
        
        
      } else {
        cat("Continue with Mann-Whitney-U test\n")
        
        
        summary_temp <- wilcox.test(DegA~type, data = subset_temp, 
                                    exact = FALSE, 
                                    correct = FALSE, 
                                    conf.int = FALSE)
        W_value <- summary_temp$statistic
        #U-value 
        #subsets
        a <-filter(subset_temp, type == "MO")
        b <-filter(subset_temp, type == "LO")
        #ranked value
        a <- a$DegA
        b <- b$DegA
        #ranks
        ranks <- rank(c(a,b))
        #nr of ranks
        a <- length(a)
        b <- length(b)
        #rank sum
        R1 <- sum(ranks[1:a])
        R2 <- sum(ranks[(a+1):(a+b)])
        # U-Values formula
        U1 <- a * b + (a * (a + 1)) / 2 - R1
        U2 <- a * b + (b * (b + 1)) / 2 - R2
        #choose the smaller U value
        U_value <- min(U1, U2)
        
        
        if (summary_temp$p.value<0.05) {
          effsize_temp <- wilcox_effsize(DegA~type, data = subset_temp, 
                                         exact = FALSE, 
                                         correct = FALSE, 
                                         conf.int = FALSE)
           pv <- paste0(toupper(effsize_temp$magnitude),", r = " , sub("^0+", "",sprintf("%.2f", effsize_temp$effsize)), ", U = " ,U_value,  ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
          temp_df <- rbind(temp_df,pv)
        } else {
          pv <-  paste0("INSIG., U = ",U_value, ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
          temp_df <- rbind(temp_df,pv)
        }
        
        
      }
    }else{
      pv <- "NA"
      temp_df <- rbind(temp_df,pv)
      
    }
    
  }
  rownames(temp_df) <- genera_list
  colnames(temp_df) <- test_name
  assign(test_name,temp_df)
  
  filename <- paste0(wd,"/save/Alg/",test_name,".txt")
  write.table(get(test_name), file = filename, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
}
####Differences betweem Mod/Unmod####
####Stats####

#MO males
{
  Alg_temp <- (Alg)
  Alg_temp$DegA <-as.numeric(Alg_temp$DegA)
  genera_list <- unique(Alg$genus)
  
  #Name the test
  test_name <- "Alg_MO_m"
  #don't forget to adjust subset filter settings! 
  #currently 13 lines down from this one!!!
  
  #start analysis:
  {
    result <- data.frame(matrix(ncol = 0, nrow = length(genera_list)))
    
    i<- 1
    for (i in 1:length(genera_list)) {
      
      temp_df <- data.frame()
      #temp_df <- colnames(temp_df) <- genera_list[i]
      j<- 10
      for (j in 1:length(genera_list)) {
        subset_temp <- filter(Alg_temp, (genus == genera_list[i] | genus == genera_list[j]) & type == "MO" & sex == "m")
        
        if (length(unique(subset_temp$genus))==2){
          #Test for normal diAlgibution via qqplotr
          data <- as.data.frame(subset_temp$DegA)
          gg <- ggplot(subset_temp = data, mapping = aes(sample = subset_temp$DegA)) +
            stat_qq_band() +
            stat_qq_line() +
            stat_qq_point() +
            labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
          
          print(gg)
          
          #User Input to confirm normal diAlgibution in qqplot
          question_qqplot <- function(question) {
            answer <- readline(prompt = question)
            while (!(answer %in% c("y", "n"))) {
              cat("Is the data normally diAlgibuted according to ggplot.\n answer: (y/n) ")
              answer <- readline(prompt = question)
            }
            return(answer)
          }
          
          answer <- question_qqplot("Is the data normally diAlgibuted according to ggplot.\n answer: (y/n) ")
          
          if (answer == "y") {
            cat("Continue with population varinace test and corresponding t-test\n")
            
            var_result <- var.test(DegA ~ genus, subset_temp, 
                                   alternative = "two.sided")
            var_result <- var_result$p.value
            
            if (var_result>0.05) {
              cat("no significant difference in varianz -> >Student's t-test")
              
              a <-filter(subset_temp, genus == genera_list[i])
              b <-filter(subset_temp, genus == genera_list[j])
              
              summary_temp <- t.test(a$DegA,b$DegA, var.equal = TRUE)
              #pv <- paste0 (summary_temp$p.value,"s")
              #temp_df <- rbind(temp_df,pv)
              if (summary_temp$p.value<0.05) {
                
                effsize_temp <- cohens_d(subset_temp, DegA ~ genus, var.equal = TRUE)
                
                            pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
                temp_df <- rbind(temp_df,pv)
              } else {
                pv <- paste0 ("INSIG., " , " s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
                temp_df <- rbind(temp_df,pv)
              }
              
            }else{
              cat("significant difference in varianz -> >Welch's t-test")
              
              a <-filter(subset_temp, genus == genera_list[i])
              b <-filter(subset_temp, genus == genera_list[j])
              summary_temp <- t.test(a$DegA,b$DegA, var.equal = FALSE)
              #pv <- paste0(summary_temp$p.value,"w")
              #temp_df <- rbind(temp_df,pv)
              
              if (summary_temp$p.value<0.05) {
                
                
                effsize_temp <- cohens_d(data = subset_temp, DegA ~ genus, var.equal = FALSE)
                
               pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", w-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
                temp_df <- rbind(temp_df,pv)
              } else {
                pv <- paste0 ("INSIG., " , " w-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value))) 
                temp_df <- rbind(temp_df,pv)
              }
              
            }
            
            
          } else {
            cat("Continue with Mann-Whitney-U test\n")
            
            
            summary_temp <- wilcox.test(DegA~genus, data = subset_temp, 
                                        exact = FALSE, 
                                        correct = FALSE, 
                                        conf.int = FALSE)
            W_value <- summary_temp$statistic
            #U-value 
            #subsets
            a <-filter(subset_temp, genus == genera_list[i])
            b <-filter(subset_temp, genus == genera_list[j])
            #ranked value
            a <- a$DegA
            b <- b$DegA
            #ranks
            ranks <- rank(c(a,b))
            #nr of ranks
            a <- length(a)
            b <- length(b)
            #rank sum
            R1 <- sum(ranks[1:a])
            R2 <- sum(ranks[(a+1):(a+b)])
            # U-Values formula
            U1 <- a * b + (a * (a + 1)) / 2 - R1
            U2 <- a * b + (b * (b + 1)) / 2 - R2
            #choose the smaller U value
            U_value <- min(U1, U2)
            if (summary_temp$p.value<0.05) {
              effsize_temp <- wilcox_effsize(DegA~genus, data = subset_temp, 
                                             exact = FALSE, 
                                             correct = FALSE, 
                                             conf.int = FALSE)
               pv <- paste0(toupper(effsize_temp$magnitude),", r = " , sub("^0+", "",sprintf("%.2f", effsize_temp$effsize)), ", U = " ,U_value,  ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
              temp_df <- rbind(temp_df,pv)
            } else {
              pv <-  paste0("INSIG., U = ",U_value, ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
              temp_df <- rbind(temp_df,pv)
            }
            
            
          }
        }else{
          pv <- "NA"
          temp_df <- rbind(temp_df,pv)
        }
      }
      result <- cbind(result,temp_df)
    }
    rownames(result) <- genera_list
    colnames(result) <- genera_list
    
  }
  
  assign(test_name,result)
  filename <- paste0(wd,"/save/Alg/",test_name,".txt")
  write.table(get(test_name), file = filename, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
}
#MO females 
{
  Alg_temp <- (Alg)
  Alg_temp$DegA <-as.numeric(Alg_temp$DegA)
  genera_list <- unique(Alg$genus)
  
  #Name the test
  test_name <- "Alg_MO_f"
  #don't forget to adjust subset filter settings! 
  #currently 13 lines down from this one!!!
  
  #start analysis:
  {
    result <- data.frame(matrix(ncol = 0, nrow = length(genera_list)))
    
    #i<- 1
    for (i in 1:length(genera_list)) {
      
      temp_df <- data.frame()
      #temp_df <- colnames(temp_df) <- genera_list[i]
      #j<- 2
      for (j in 1:length(genera_list)) {
        subset_temp <- filter(Alg_temp, (genus == genera_list[i] | genus == genera_list[j]) & (type == "MO") & sex == "f")
        
        if (length(unique(subset_temp$genus))==2){
          #Test for normal diAlgibution via qqplotr
          data <- as.data.frame(subset_temp$DegA)
          gg <- ggplot(subset_temp = data, mapping = aes(sample = subset_temp$DegA)) +
            stat_qq_band() +
            stat_qq_line() +
            stat_qq_point() +
            labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
          
          print(gg)
          
          #User Input to confirm normal diAlgibution in qqplot
          question_qqplot <- function(question) {
            answer <- readline(prompt = question)
            while (!(answer %in% c("y", "n"))) {
              cat("Is the data normally diAlgibuted according to ggplot.\n answer: (y/n) ")
              answer <- readline(prompt = question)
            }
            return(answer)
          }
          
          answer <- question_qqplot("Is the data normally diAlgibuted according to ggplot.\n answer: (y/n) ")
          
          if (answer == "y") {
            cat("Continue with population varinace test and corresponding t-test\n")
            
            var_result <- var.test(DegA ~ genus, subset_temp, 
                                   alternative = "two.sided")
            var_result <- var_result$p.value
            
            if (var_result>0.05) {
              cat("no significant difference in varianz -> >Student's t-test")
              
              a <-filter(subset_temp, genus == genera_list[i])
              b <-filter(subset_temp, genus == genera_list[j])
              
              summary_temp <- t.test(a$DegA,b$DegA, var.equal = TRUE)
              #pv <- paste0 (summary_temp$p.value,"s")
              #temp_df <- rbind(temp_df,pv)
              if (summary_temp$p.value<0.05) {
                
                effsize_temp <- cohens_d(subset_temp, DegA ~ genus, var.equal = TRUE)
                
                            pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
                temp_df <- rbind(temp_df,pv)
              } else {
                pv <- paste0 ("INSIG., " , " s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
                temp_df <- rbind(temp_df,pv)
              }
              
            }else{
              cat("significant difference in varianz -> >Welch's t-test")
              
              a <-filter(subset_temp, genus == genera_list[i])
              b <-filter(subset_temp, genus == genera_list[j])
              summary_temp <- t.test(a$DegA,b$DegA, var.equal = FALSE)
              #pv <- paste0(summary_temp$p.value,"w")
              #temp_df <- rbind(temp_df,pv)
              
              if (summary_temp$p.value<0.05) {
                
                effsize_temp <- cohens_d(subset_temp, DegA ~ genus, var.equal = FALSE)
                
               pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", w-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
                temp_df <- rbind(temp_df,pv)
              } else {
                pv <- paste0 ("INSIG., " , " w-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value))) 
                temp_df <- rbind(temp_df,pv)
              }
              
            }
            
            
          } else {
            cat("Continue with Mann-Whitney-U test\n")
            
            
            summary_temp <- wilcox.test(DegA~genus, data = subset_temp, 
                                        exact = FALSE, 
                                        correct = FALSE, 
                                        conf.int = FALSE)
            W_value <- summary_temp$statistic
            #U-value 
            #subsets
            a <-filter(subset_temp, genus == genera_list[i])
            b <-filter(subset_temp, genus == genera_list[j])
            #ranked value
            a <- a$DegA
            b <- b$DegA
            #ranks
            ranks <- rank(c(a,b))
            #nr of ranks
            a <- length(a)
            b <- length(b)
            #rank sum
            R1 <- sum(ranks[1:a])
            R2 <- sum(ranks[(a+1):(a+b)])
            # U-Values formula
            U1 <- a * b + (a * (a + 1)) / 2 - R1
            U2 <- a * b + (b * (b + 1)) / 2 - R2
            #choose the smaller U value
            U_value <- min(U1, U2)
            
            if (summary_temp$p.value<0.05) {
              effsize_temp <- wilcox_effsize(DegA~genus, data = subset_temp, 
                                             exact = FALSE, 
                                             correct = FALSE, 
                                             conf.int = FALSE)
               pv <- paste0(toupper(effsize_temp$magnitude),", r = " , sub("^0+", "",sprintf("%.2f", effsize_temp$effsize)), ", U = " ,U_value,  ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
              temp_df <- rbind(temp_df,pv)
            } else {
              pv <-  paste0("INSIG., U = ",U_value, ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
              temp_df <- rbind(temp_df,pv)
            }
            
            
          }
        }else{
          pv <- "NA"
          temp_df <- rbind(temp_df,pv)
        }
      }
      result <- cbind(result,temp_df)
    }
    rownames(result) <- genera_list
    colnames(result) <- genera_list
    
  }
  
  assign(test_name,result)
  filename <- paste0(wd,"/save/Alg/",test_name,".txt")
  write.table(get(test_name), file = filename, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
}
#LO males
{
  Alg_temp <- (Alg)
  Alg_temp$DegA <-as.numeric(Alg_temp$DegA)
  genera_list <- unique(Alg$genus)
  
  #Name the test
  test_name <- "Alg_LO_m"
  #don't forget to adjust subset filter settings! 
  #currently 13 lines down from this one!!!
  
  #start analysis:
  {
    result <- data.frame(matrix(ncol = 0, nrow = length(genera_list)))
    
    #i<- 1
    for (i in 1:length(genera_list)) {
      
      temp_df <- data.frame()
      #temp_df <- colnames(temp_df) <- genera_list[i]
      #j<- 2
      for (j in 1:length(genera_list)) {
        subset_temp <- filter(Alg_temp, (genus == genera_list[i] | genus == genera_list[j]) & type == "LO" & sex == "m")
        
        if (length(unique(subset_temp$genus))==2){
          #Test for normal diAlgibution via qqplotr
          data <- as.data.frame(subset_temp$DegA)
          gg <- ggplot(subset_temp = data, mapping = aes(sample = subset_temp$DegA)) +
            stat_qq_band() +
            stat_qq_line() +
            stat_qq_point() +
            labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
          
          print(gg)
          
          #User Input to confirm normal diAlgibution in qqplot
          question_qqplot <- function(question) {
            answer <- readline(prompt = question)
            while (!(answer %in% c("y", "n"))) {
              cat("Is the data normally diAlgibuted according to ggplot.\n answer: (y/n) ")
              answer <- readline(prompt = question)
            }
            return(answer)
          }
          
          answer <- question_qqplot("Is the data normally diAlgibuted according to ggplot.\n answer: (y/n) ")
          
          if (answer == "y") {
            cat("Continue with population varinace test and corresponding t-test\n")
            
            var_result <- var.test(DegA ~ genus, subset_temp, 
                                   alternative = "two.sided")
            var_result <- var_result$p.value
            
            if (var_result>0.05) {
              cat("no significant difference in varianz -> >Student's t-test")
              
              a <-filter(subset_temp, genus == genera_list[i])
              b <-filter(subset_temp, genus == genera_list[j])
              
              summary_temp <- t.test(a$DegA,b$DegA, var.equal = TRUE)
              #pv <- paste0 (summary_temp$p.value,"s")
              #temp_df <- rbind(temp_df,pv)
              if (summary_temp$p.value<0.05) {
                
                effsize_temp <- cohens_d(subset_temp, DegA ~ genus, var.equal = TRUE)
                
                            pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
                temp_df <- rbind(temp_df,pv)
              } else {
                pv <- paste0 ("INSIG., " , " s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
                temp_df <- rbind(temp_df,pv)
              }
              
            }else{
              cat("significant difference in varianz -> >Welch's t-test")
              
              a <-filter(subset_temp, genus == genera_list[i])
              b <-filter(subset_temp, genus == genera_list[j])
              summary_temp <- t.test(a$DegA,b$DegA, var.equal = FALSE)
              #pv <- paste0(summary_temp$p.value,"w")
              #temp_df <- rbind(temp_df,pv)
              
              if (summary_temp$p.value<0.05) {
                
                effsize_temp <- cohens_d(subset_temp, DegA ~ genus, var.equal = FALSE)
                
               pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", w-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
                temp_df <- rbind(temp_df,pv)
              } else {
                pv <- paste0 ("INSIG., " , " w-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value))) 
                temp_df <- rbind(temp_df,pv)
              }
              
            }
            
            
          } else {
            cat("Continue with Mann-Whitney-U test\n")
            
            
            summary_temp <- wilcox.test(DegA~genus, data = subset_temp, 
                                        exact = FALSE, 
                                        correct = FALSE, 
                                        conf.int = FALSE)
            W_value <- summary_temp$statistic
            #U-value 
            #subsets
            a <-filter(subset_temp, genus == genera_list[i])
            b <-filter(subset_temp, genus == genera_list[j])
            #ranked value
            a <- a$DegA
            b <- b$DegA
            #ranks
            ranks <- rank(c(a,b))
            #nr of ranks
            a <- length(a)
            b <- length(b)
            #rank sum
            R1 <- sum(ranks[1:a])
            R2 <- sum(ranks[(a+1):(a+b)])
            # U-Values formula
            U1 <- a * b + (a * (a + 1)) / 2 - R1
            U2 <- a * b + (b * (b + 1)) / 2 - R2
            #choose the smaller U value
            U_value <- min(U1, U2)
            if (summary_temp$p.value<0.05) {
              effsize_temp <- wilcox_effsize(DegA~genus, data = subset_temp, 
                                             exact = FALSE, 
                                             correct = FALSE, 
                                             conf.int = FALSE)
               pv <- paste0(toupper(effsize_temp$magnitude),", r = " , sub("^0+", "",sprintf("%.2f", effsize_temp$effsize)), ", U = " ,U_value,  ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
              temp_df <- rbind(temp_df,pv)
            } else {
              pv <-  paste0("INSIG., U = ",U_value, ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
              temp_df <- rbind(temp_df,pv)
            }
            
            
          }
        }else{
          pv <- "NA"
          temp_df <- rbind(temp_df,pv)
        }
      }
      result <- cbind(result,temp_df)
    }
    rownames(result) <- genera_list
    colnames(result) <- genera_list
    
  }
  
  assign(test_name,result)
  filename <- paste0(wd,"/save/Alg/",test_name,".txt")
  write.table(get(test_name), file = filename, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
}
#LO females
{
  Alg_temp <- (Alg)
  Alg_temp$DegA <-as.numeric(Alg_temp$DegA)
  genera_list <- unique(Alg$genus)
  
  #Name the test
  test_name <- "Alg_LO_f"
  #don't forget to adjust subset filter settings! 
  #currently 13 lines down from this one!!!
  
  #start analysis:
  {
    result <- data.frame(matrix(ncol = 0, nrow = length(genera_list)))
    
    #i<- 1
    for (i in 1:length(genera_list)) {
      
      temp_df <- data.frame()
      #temp_df <- colnames(temp_df) <- genera_list[i]
      #j<- 2
      for (j in 1:length(genera_list)) {
        subset_temp <- filter(Alg_temp, (genus == genera_list[i] | genus == genera_list[j]) & type == "LO" & sex == "f")
        
        if (length(unique(subset_temp$genus))==2){
          #Test for normal diAlgibution via qqplotr
          data <- as.data.frame(subset_temp$DegA)
          gg <- ggplot(subset_temp = data, mapping = aes(sample = subset_temp$DegA)) +
            stat_qq_band() +
            stat_qq_line() +
            stat_qq_point() +
            labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
          
          print(gg)
          
          #User Input to confirm normal diAlgibution in qqplot
          question_qqplot <- function(question) {
            answer <- readline(prompt = question)
            while (!(answer %in% c("y", "n"))) {
              cat("Is the data normally diAlgibuted according to ggplot.\n answer: (y/n) ")
              answer <- readline(prompt = question)
            }
            return(answer)
          }
          
          answer <- question_qqplot("Is the data normally diAlgibuted according to ggplot.\n answer: (y/n) ")
          
          if (answer == "y") {
            cat("Continue with population varinace test and corresponding t-test\n")
            
            var_result <- var.test(DegA ~ genus, subset_temp, 
                                   alternative = "two.sided")
            var_result <- var_result$p.value
            
            if (var_result>0.05) {
              cat("no significant difference in varianz -> >Student's t-test")
              
              a <-filter(subset_temp, genus == genera_list[i])
              b <-filter(subset_temp, genus == genera_list[j])
              
              summary_temp <- t.test(a$DegA,b$DegA, var.equal = TRUE)
              #pv <- paste0 (summary_temp$p.value,"s")
              #temp_df <- rbind(temp_df,pv)
              if (summary_temp$p.value<0.05) {
                
                effsize_temp <- cohens_d(subset_temp, DegA ~ genus, var.equal = TRUE)
                
                pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
                temp_df <- rbind(temp_df,pv)
              } else {
                pv <- paste0 ("INSIG., " , " s-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
                temp_df <- rbind(temp_df,pv)
              }
              
            }else{
              cat("significant difference in varianz -> >Welch's t-test")
              
              a <-filter(subset_temp, genus == genera_list[i])
              b <-filter(subset_temp, genus == genera_list[j])
              summary_temp <- t.test(a$DegA,b$DegA, var.equal = FALSE)
              #pv <- paste0(summary_temp$p.value,"w")
              #temp_df <- rbind(temp_df,pv)
              
              if (summary_temp$p.value<0.05) {
                
                effsize_temp <- cohens_d(subset_temp, DegA ~ genus, var.equal = FALSE)
                
                pv <- paste0(toupper(effsize_temp$magnitude),", d = ",round(effsize_temp$effsize,2), ", w-t(" ,round(summary_temp$parameter,1), "),  = " , round (summary_temp$statistic, 2), ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
                temp_df <- rbind(temp_df,pv)
              } else {
                pv <- paste0(round(summary_temp$p.value,4), ", INSIG, t: ",round(summary_temp$statistic,2), ", df: ",round(summary_temp$parameter,4), ",w")
                temp_df <- rbind(temp_df,pv)
              }
              
            }
            
            
          } else {
            cat("Continue with Mann-Whitney-U test\n")
            
            
            summary_temp <- wilcox.test(DegA~genus, data = subset_temp, 
                                        exact = FALSE, 
                                        correct = FALSE, 
                                        conf.int = FALSE)
            W_value <- summary_temp$statistic
            #U-value 
            #subsets
            a <-filter(subset_temp, genus == genera_list[i])
            b <-filter(subset_temp, genus == genera_list[j])
            #ranked value
            a <- a$DegA
            b <- b$DegA
            #ranks
            ranks <- rank(c(a,b))
            #nr of ranks
            a <- length(a)
            b <- length(b)
            #rank sum
            R1 <- sum(ranks[1:a])
            R2 <- sum(ranks[(a+1):(a+b)])
            # U-Values formula
            U1 <- a * b + (a * (a + 1)) / 2 - R1
            U2 <- a * b + (b * (b + 1)) / 2 - R2
            #choose the smaller U value
            U_value <- min(U1, U2)
            if (summary_temp$p.value<0.05) {
              effsize_temp <- wilcox_effsize(DegA~genus, data = subset_temp, 
                                             exact = FALSE, 
                                             correct = FALSE, 
                                             conf.int = FALSE)
               pv <- paste0(toupper(effsize_temp$magnitude),", r = " , sub("^0+", "",sprintf("%.2f", effsize_temp$effsize)), ", U = " ,U_value,  ", p = ", sub("^0+", "",sprintf("%.3f", summary_temp$p.value)))
              temp_df <- rbind(temp_df,pv)
            } else {
              pv <-  paste0("INSIG., U = ",U_value, ", p = ", sub("^0+", "",sprintf("%.2f", summary_temp$p.value)))
              temp_df <- rbind(temp_df,pv)
            }
            
            
          }
        }else{
          pv <- "NA"
          temp_df <- rbind(temp_df,pv)
        }
      }
      result <- cbind(result,temp_df)
    }
    rownames(result) <- genera_list
    colnames(result) <- genera_list
    
  }
  
  assign(test_name,result)
  filename <- paste0(wd,"/save/Alg/",test_name,".txt")
  write.table(get(test_name), file = filename, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
}
####Summary by mod/unmod####

#female
subset_temp <-  filter(Alg, sex == "f")
subset_n1 <- filter(subset_temp, mod == "y")
print(length(subset_n1$DegA))
print(mean(subset_n1$DegA))
print(sd(subset_n1$DegA))
subset_n2 <- filter(subset_temp, mod == "n")
print(length(subset_n2$DegA))
print(mean(subset_n2$DegA))
print(sd(subset_n2$DegA))

subset_n1 <- filter(subset_temp, mod == "n")
gg <- ggplot(data = subset_temp, mapping = aes(sample = subset_temp$DegA)) +
  stat_qq_band() +
  stat_qq_line() +
  stat_qq_point() +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")

print(gg)
# nicht normal verteilt

summary_temp <- wilcox.test(DegA~mod, data = subset_temp, 
                            exact = FALSE, 
                            correct = FALSE, 
                            conf.int = FALSE)
W_value <- summary_temp$statistic
#U-value 
#subsets
a <-filter(subset_temp, mod == "y")
b <-filter(subset_temp, mod == "n")
#ranked value
a <- a$DegA
b <- b$DegA
#ranks
ranks <- rank(c(a,b))
#nr of ranks
a <- length(a)
b <- length(b)
#rank sum
R1 <- sum(ranks[1:a])
R2 <- sum(ranks[(a+1):(a+b)])
# U-Values formula
U1 <- a * b + (a * (a + 1)) / 2 - R1
U2 <- a * b + (b * (b + 1)) / 2 - R2
#choose the smaller U value
U_value <- min(U1, U2)


effsize_temp <- wilcox_effsize(DegA~mod, data = subset_temp, 
                               exact = FALSE, 
                               correct = FALSE, 
                               conf.int = FALSE)
Subs_mean <- filter(subset_temp, mod =="y")
mean(Subs_mean$DegA)
sd(Subs_mean$DegA)
Subs_mean <- filter(subset_temp, mod =="n")
mean(Subs_mean$DegA)
sd(Subs_mean$DegA)

#male
subset_temp <-  filter(Alg, sex == "m")
subset_n1 <- filter(subset_temp, mod == "y")
print(length(subset_n1$DegA))
print(mean(subset_n1$DegA))
print(sd(subset_n1$DegA))
subset_n2 <- filter(subset_temp, mod == "n")
print(length(subset_n2$DegA))
print(mean(subset_n2$DegA))
print(sd(subset_n2$DegA))

gg <- ggplot(data = subset_temp, mapping = aes(sample = subset_temp$DegA)) +
  stat_qq_band() +
  stat_qq_line() +
  stat_qq_point() +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")

print(gg)
# nicht normal verteilt

summary_temp <- wilcox.test(DegA~mod, data = subset_temp, 
                            exact = FALSE, 
                            correct = FALSE, 
                            conf.int = FALSE)
W_value <- summary_temp$statistic
#U-value 
#subsets
a <-filter(subset_temp, mod == "y")
b <-filter(subset_temp, mod == "n")
#ranked value
a <- a$DegA
b <- b$DegA
#ranks
ranks <- rank(c(a,b))
#nr of ranks
a <- length(a)
b <- length(b)
#rank sum
R1 <- sum(ranks[1:a])
R2 <- sum(ranks[(a+1):(a+b)])
# U-Values formula
U1 <- a * b + (a * (a + 1)) / 2 - R1
U2 <- a * b + (b * (b + 1)) / 2 - R2
#choose the smaller U value
U_value <- min(U1, U2)



effsize_temp <- wilcox_effsize(DegA~mod, data = subset_temp, 
                               exact = FALSE, 
                               correct = FALSE, 
                               conf.int = FALSE)
Subs_mean <- filter(subset_temp, mod =="y")
mean(Subs_mean$DegA)
sd(Subs_mean$DegA)
Subs_mean <- filter(subset_temp, mod =="n")
mean(Subs_mean$DegA)
sd(Subs_mean$DegA)

####------------------------------------------------------------------####
####Rhabdom length####
####Difference Mod/Unmod

#female
subset_temp <-  filter(Str, sex == "f")
subset_n1 <- filter(subset_temp, mod == "y")
print(length(subset_n1$rhlength))
print(mean(subset_n1$rhlength))
print(sd(subset_n1$rhlength))
subset_n2 <- filter(subset_temp, mod == "n")
print(length(subset_n2$rhlength))
print(mean(subset_n2$rhlength))
print(sd(subset_n2$rhlength))

subset_n1 <- filter(subset_temp, mod == "n")
gg <- ggplot(data = subset_temp, mapping = aes(sample = subset_temp$rhlength)) +
  stat_qq_band() +
  stat_qq_line() +
  stat_qq_point() +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")

print(gg)
# nicht normal verteilt

summary_temp <- wilcox.test(rhlength~mod, data = subset_temp, 
                            exact = FALSE, 
                            correct = FALSE, 
                            conf.int = FALSE)
W_value <- summary_temp$statistic
#U-value 
#subsets
a <-filter(subset_temp, mod == "y")
b <-filter(subset_temp, mod == "n")
#ranked value
a <- a$rhlength
b <- b$rhlength
#ranks
ranks <- rank(c(a,b))
#nr of ranks
a <- length(a)
b <- length(b)
#rank sum
R1 <- sum(ranks[1:a])
R2 <- sum(ranks[(a+1):(a+b)])
# U-Values formula
U1 <- a * b + (a * (a + 1)) / 2 - R1
U2 <- a * b + (b * (b + 1)) / 2 - R2
#choose the smaller U value
U_value <- min(U1, U2)


effsize_temp <- wilcox_effsize(rhlength~mod, data = subset_temp, 
                                 exact = FALSE, 
                                 correct = FALSE, 
                                 conf.int = FALSE)
Subs_mean <- filter(subset_temp, mod =="y")
mean(Subs_mean$rhlength)
sd(Subs_mean$rhlength)
Subs_mean <- filter(subset_temp, mod =="n")
mean(Subs_mean$rhlength)
sd(Subs_mean$rhlength)

#male
subset_temp <-  filter(Str, sex == "m")
subset_n1 <- filter(subset_temp, mod == "y")
print(length(subset_n1$rhlength))
print(mean(subset_n1$rhlength))
print(sd(subset_n1$rhlength))
subset_n2 <- filter(subset_temp, mod == "n")
print(length(subset_n2$rhlength))
print(mean(subset_n2$rhlength))
print(sd(subset_n2$rhlength))

gg <- ggplot(data = subset_temp, mapping = aes(sample = subset_temp$rhlength)) +
  stat_qq_band() +
  stat_qq_line() +
  stat_qq_point() +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")

print(gg)
# nicht normal verteilt

summary_temp <- wilcox.test(rhlength~mod, data = subset_temp, 
                            exact = FALSE, 
                            correct = FALSE, 
                            conf.int = FALSE)
W_value <- summary_temp$statistic
#U-value 
#subsets
a <-filter(subset_temp, mod == "y")
b <-filter(subset_temp, mod == "n")
#ranked value
a <- a$rhlength
b <- b$rhlength
#ranks
ranks <- rank(c(a,b))
#nr of ranks
a <- length(a)
b <- length(b)
#rank sum
R1 <- sum(ranks[1:a])
R2 <- sum(ranks[(a+1):(a+b)])
# U-Values formula
U1 <- a * b + (a * (a + 1)) / 2 - R1
U2 <- a * b + (b * (b + 1)) / 2 - R2
#choose the smaller U value
U_value <- min(U1, U2)



effsize_temp <- wilcox_effsize(rhlength~mod, data = subset_temp, 
                               exact = FALSE, 
                               correct = FALSE, 
                               conf.int = FALSE)
Subs_mean <- filter(subset_temp, mod =="y")
mean(Subs_mean$rhlength)
sd(Subs_mean$rhlength)
Subs_mean <- filter(subset_temp, mod =="n")
mean(Subs_mean$rhlength)
sd(Subs_mean$rhlength)


####------------------------------------------------------------------####
####Str Vis####
#MO/LO females
{subset_vis <- filter(Str_temp, sex == "f")
p_raincloud <- ggplot(data = subset_vis, aes(x = interaction(type,genus), y = DegS, fill = mod))+
  
  geom_flat_violin( width = 1.25, position=position_nudge(x=0.15, y=0), alpha = 0.5) +
  
  geom_boxplot(color = "white", outlier.shape=NA, alpha=1, width=0.075, position=position_nudge(x=0.15, y=0)) +
  
  geom_point(shape = 16, color = "white", position=position_jitter(width=.075, height=0), size=0.4) +
  
  scale_fill_manual(values = c("n" = "blue", "y" = "orange"))+
  
  theme_classic() +
  
  theme(axis.title.x=element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
        axis.title.y=element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
        panel.background = element_rect(fill = "black"),
        plot.background = element_rect(fill = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(color = "white"),
        axis.text = element_text(color = "white", size = 30),
        axis.ticks = element_line(color = "white"),
        legend.background = element_rect(fill = "black"),
        legend.text = element_text(color = "white", size = 20),
        legend.title = element_text(color = "white", size = 20),
        plot.title = element_text(color = "white", hjust = 0.5),
        plot.subtitle = element_text(color = "white", hjust = 0.5),
        plot.caption = element_text(color = "white"),
        axis.line = element_line(color = "white")
  ) +
  
  guides(fill = guide_legend(title = "Modified ocelli")) + 
  
  scale_x_discrete(expand = expansion(mult = c(0.04, 0.0))) +
  
  ylab("mean deviation in °" ) +
  
  xlab("") +
  
  coord_flip() 
print(p_raincloud)
ggsave(paste0(wd,"/save/Str/","Str_MOuLO_females",".png"), width = 12, height = 12, plot = p_raincloud, dpi = 300)
}
#MO/LO males
{subset_vis <- filter(Str_temp, sex == "m")
  p_raincloud <- ggplot(data = subset_vis, aes(x = interaction(type,genus), y = DegS, fill = mod))+
    
    geom_flat_violin( width = 1.25, position=position_nudge(x=0.15, y=0), alpha = 0.5) +
    
    geom_boxplot(color = "white", outlier.shape=NA, alpha=1, width=0.075, position=position_nudge(x=0.15, y=0)) +
    
    geom_point(shape = 16, color = "white", position=position_jitter(width=.075, height=0), size=0.4) +
    
    scale_fill_manual(values = c("n" = "blue", "y" = "orange"))+
    
    theme_classic() +
    
    theme(axis.title.x=element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
          axis.title.y=element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
          panel.background = element_rect(fill = "black"),
          plot.background = element_rect(fill = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(color = "white"),
          axis.text = element_text(color = "white", size = 30),
          axis.ticks = element_line(color = "white"),
          legend.background = element_rect(fill = "black"),
          legend.text = element_text(color = "white", size = 20),
          legend.title = element_text(color = "white", size = 20),
          plot.title = element_text(color = "white", hjust = 0.5),
          plot.subtitle = element_text(color = "white", hjust = 0.5),
          plot.caption = element_text(color = "white"),
          axis.line = element_line(color = "white")
    ) +
    
    guides(fill = guide_legend(title = "Modified ocelli")) + 
    
    scale_x_discrete(expand = expansion(mult = c(0.04, 0.0))) +
    
    ylab("mean deviation in °" ) +
    
    xlab("") +
    
    coord_flip() 
  print(p_raincloud)
  ggsave(paste0(wd,"/save/Str/","Str_MOuLO_females",".png"), width = 12, height = 12, plot = p_raincloud, dpi = 300)
}
#MO females
{subset_vis <- filter(Str_temp, (type == "MO"| type == "hMO" | type == "cMO") & sex == "f")
  p_raincloud <- ggplot(data = subset_vis, aes(x = genus, y = DegS, fill = mod))+
    
    geom_flat_violin( width = 1.25, position=position_nudge(x=0.15, y=0), alpha = 0.5) +
    
    geom_boxplot(color = "white", outlier.shape=NA, alpha=1, width=0.075, position=position_nudge(x=0.15, y=0)) +
    
    geom_point(shape = 16, color = "white", position=position_jitter(width=.075, height=0), size=0.4) +
    
    scale_fill_manual(values = c("n" = "blue", "y" = "orange"))+
    
    theme_classic() +
    
    theme(axis.title.x=element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
          axis.title.y=element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
          panel.background = element_rect(fill = "black"),
          plot.background = element_rect(fill = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(color = "white"),
          axis.text = element_text(color = "white", size = 30),
          axis.ticks = element_line(color = "white"),
          legend.background = element_rect(fill = "black"),
          legend.text = element_text(color = "white", size = 20),
          legend.title = element_text(color = "white", size = 20),
          plot.title = element_text(color = "white", hjust = 0.5),
          plot.subtitle = element_text(color = "white", hjust = 0.5),
          plot.caption = element_text(color = "white"),
          axis.line = element_line(color = "white")
    ) +
    
    guides(fill = guide_legend(title = "Modified ocelli")) + 
    
    scale_x_discrete(expand = expansion(mult = c(0.04, 0.0))) +
    
    ylab("mean deviation in °" ) +
    
    xlab("") +
    
    coord_flip() 
  print(p_raincloud)
  ggsave(paste0(wd,"/save/Str/","Str_MO_females",".png"), width = 12, height = 12, plot = p_raincloud, dpi = 300)
}
#MO males
{subset_vis <- filter(Str_temp, (type == "MO"| type == "hMO" | type == "cMO") & sex == "m")
  p_raincloud <- ggplot(data = subset_vis, aes(x = genus, y = DegS, fill = mod))+
    
    geom_flat_violin( width = 1.25, position=position_nudge(x=0.15, y=0), alpha = 0.5) +
    
    geom_boxplot(color = "white", outlier.shape=NA, alpha=1, width=0.075, position=position_nudge(x=0.15, y=0)) +
    
    geom_point(shape = 16, color = "white", position=position_jitter(width=.075, height=0), size=0.4) +
    
    scale_fill_manual(values = c("n" = "blue", "y" = "orange"))+
    
    theme_classic() +
    
    theme(axis.title.x=element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
          axis.title.y=element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
          panel.background = element_rect(fill = "black"),
          plot.background = element_rect(fill = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(color = "white"),
          axis.text = element_text(color = "white", size = 30),
          axis.ticks = element_line(color = "white"),
          legend.background = element_rect(fill = "black"),
          legend.text = element_text(color = "white", size = 20),
          legend.title = element_text(color = "white", size = 20),
          plot.title = element_text(color = "white", hjust = 0.5),
          plot.subtitle = element_text(color = "white", hjust = 0.5),
          plot.caption = element_text(color = "white"),
          axis.line = element_line(color = "white")
    ) +
    
    guides(fill = guide_legend(title = "Modified ocelli")) + 
    
    scale_x_discrete(expand = expansion(mult = c(0.04, 0.0))) +
    
    ylab("mean deviation in °") +
    
    xlab("") +
    
    coord_flip() 
  print(p_raincloud)
  ggsave(paste0(wd,"/save/Str/","Str_MO_males",".png"), width = 12, height = 12, plot = p_raincloud, dpi = 300)
}
#LO females
{subset_vis <- filter(Str_temp, (type == "LO") & sex == "f")
  p_raincloud <- ggplot(data = subset_vis, aes(x = genus, y = DegS, fill = mod))+
    
    geom_flat_violin( width = 1.25, position=position_nudge(x=0.15, y=0), alpha = 0.5) +
    
    geom_boxplot(color = "white", outlier.shape=NA, alpha=1, width=0.075, position=position_nudge(x=0.15, y=0)) +
    
    geom_point(shape = 16, color = "white", position=position_jitter(width=.075, height=0), size=0.4) +
    
    scale_fill_manual(values = c("n" = "blue", "y" = "orange"))+
    
    theme_classic() +
    
    theme(axis.title.x=element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
          axis.title.y=element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
          panel.background = element_rect(fill = "black"),
          plot.background = element_rect(fill = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(color = "white"),
          axis.text = element_text(color = "white", size = 30),
          axis.ticks = element_line(color = "white"),
          legend.background = element_rect(fill = "black"),
          legend.text = element_text(color = "white", size = 20),
          legend.title = element_text(color = "white", size = 20),
          plot.title = element_text(color = "white", hjust = 0.5),
          plot.subtitle = element_text(color = "white", hjust = 0.5),
          plot.caption = element_text(color = "white"),
          axis.line = element_line(color = "white")
    ) +
    
    guides(fill = guide_legend(title = "Modified ocelli")) + 
    
    scale_x_discrete(expand = expansion(mult = c(0.04, 0.0))) +
    
    ylab("mean deviation in °" ) +
    
    xlab("") +
    
    coord_flip() 
  print(p_raincloud)
  ggsave(paste0(wd,"/save/Str/","Str_LO_females",".png"), width = 12, height = 12, plot = p_raincloud, dpi = 300)
}
#LO males 
{subset_vis <- filter(Str_temp, (type == "LO") & sex == "m")
  p_raincloud <- ggplot(data = subset_vis, aes(x = genus, y = DegS, fill = mod))+
    
    geom_flat_violin( width = 1.25, position=position_nudge(x=0.15, y=0), alpha = 0.5) +
    
    geom_boxplot(color = "white", outlier.shape=NA, alpha=1, width=0.075, position=position_nudge(x=0.15, y=0)) +
    
    geom_point(shape = 16, color = "white", position=position_jitter(width=.075, height=0), size=0.4) +
    
    scale_fill_manual(values = c("n" = "blue", "y" = "orange"))+
    
    theme_classic() +
    
    theme(axis.title.x=element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
          axis.title.y=element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
          panel.background = element_rect(fill = "black"),
          plot.background = element_rect(fill = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(color = "white"),
          axis.text = element_text(color = "white", size = 30),
          axis.ticks = element_line(color = "white"),
          legend.background = element_rect(fill = "black"),
          legend.text = element_text(color = "white", size = 20),
          legend.title = element_text(color = "white", size = 20),
          plot.title = element_text(color = "white", hjust = 0.5),
          plot.subtitle = element_text(color = "white", hjust = 0.5),
          plot.caption = element_text(color = "white"),
          axis.line = element_line(color = "white")
    ) +
    
    guides(fill = guide_legend(title = "Modified ocelli")) + 
    
    scale_x_discrete(expand = expansion(mult = c(0.04, 0.0))) +
    
    ylab("mean deviation in °") +
    
    xlab("") +
    
    coord_flip() 
  print(p_raincloud)
  ggsave(paste0(wd,"/save/Str/","Str_LO_males",".png"), width = 12, height = 12, plot = p_raincloud, dpi = 300)
}
####mod####
#MO/LO females
{subset_vis <- filter(Str_temp, sex == "f")


p_raincloud <- ggplot(data = subset_vis, aes(x = mod, y = DegS, fill = mod))+
  
  geom_flat_violin( width = 1.6, position=position_nudge(x=0, y=0), alpha = 0.5) +
  
  geom_boxplot(color = "black", outlier.shape=NA, alpha=1, width=0.2, position=position_nudge(x=0, y=0)) +
  
  #geom_point(shape = 16, color = "black", position=position_jitter(width=.075, height=0), size=0.4) +
  
  scale_fill_manual(values = c("n" = "white", "y" = "#4F4F4F"))+
  
  theme_classic() +
  
  theme(axis.title.x=element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
        axis.title.y=element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        panel.grid.major.x =  element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black", size = 30),
        axis.ticks = element_line(color = "black"),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(color = "black", size = 20),
        legend.title = element_text(color = "black", size = 20),
        plot.title = element_text(color = "black", hjust = 0.5),
        plot.subtitle = element_text(color = "black", hjust = 0.5),
        plot.caption = element_text(color = "black"),
        axis.line = element_line(color = "black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()
  ) +
  
  guides(fill = guide_legend(title = "Modified ocelli")) + 
  
  stat_summary(fun = mean, geom = "crossbar", width = 0.6, color = "black", 
               linetype = "dashed", size = 0.2, position = position_nudge(x = 0, y = 0)) +
  
  ylab("rhabdom length in µm" ) +
  
  guides(fill = "none")+
  
  xlab("") +
  
  
  #scale_x_discrete(labels = labels_for_interaction, 
  #expand = expansion(mult = c(0.05, 0))) +
  scale_y_continuous(
    limits = c(0, NA),               # Behalte dein Limit von -5 und NA (automatische Obergrenze)
    breaks = seq(0, 15, by = 2.5),     # Behalte die Breaks von -5 bis 80 in 5er-Schritten
    #labels = function(x) ifelse(x %% 5 == 0 & x >= 0 & x <= 80, x, "")  # Zeige nur 0, 20, 40, 60, 80
  )+
  
  coord_flip() 
print(p_raincloud)
ggsave(paste0(wd,"/save/Str/","Str_mod",".png"), width = 12, height = 12, plot = p_raincloud, dpi = 300)
}
#MO/LO males
{subset_vis <- filter(Str_temp, sex == "m")
  
  p_raincloud <- ggplot(data = subset_vis, aes(x = mod, y = DegS, fill = mod))+
    
    geom_flat_violin( width = 1.6, position=position_nudge(x=0, y=0), alpha = 0.5) +
    
    geom_boxplot(color = "black", outlier.shape=NA, alpha=1, width=0.2, position=position_nudge(x=0, y=0)) +
    
    #geom_point(shape = 16, color = "black", position=position_jitter(width=.075, height=0), size=0.4) +
    
    scale_fill_manual(values = c("n" = "white", "y" = "#4F4F4F"))+
    
    theme_classic() +
    
    theme(axis.title.x=element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
          axis.title.y=element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"),
          panel.grid.major.x =  element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.title = element_text(color = "black"),
          axis.text = element_text(color = "black", size = 30),
          axis.ticks = element_line(color = "black"),
          legend.background = element_rect(fill = "white"),
          legend.text = element_text(color = "black", size = 20),
          legend.title = element_text(color = "black", size = 20),
          plot.title = element_text(color = "black", hjust = 0.5),
          plot.subtitle = element_text(color = "black", hjust = 0.5),
          plot.caption = element_text(color = "black"),
          axis.line = element_line(color = "black"),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank()
    ) +
    
    guides(fill = guide_legend(title = "Modified ocelli")) + 
    
    stat_summary(fun = mean, geom = "crossbar", width = 0.6, color = "black", 
                 linetype = "dashed", size = 0.2, position = position_nudge(x = 0, y = 0)) +
    
    ylab("rhabdom length in µm" ) +
    
    guides(fill = "none")+
    
    xlab("") +
    
    
    #scale_x_discrete(labels = labels_for_interaction, 
    #expand = expansion(mult = c(0.05, 0))) +
    scale_y_continuous(
      limits = c(0, NA),               # Behalte dein Limit von -5 und NA (automatische Obergrenze)
      breaks = seq(0, 15, by = 2.5),     # Behalte die Breaks von -5 bis 80 in 5er-Schritten
      #labels = function(x) ifelse(x %% 5 == 0 & x >= 0 & x <= 80, x, "")  # Zeige nur 0, 20, 40, 60, 80
    )+
    coord_flip() 
  print(p_raincloud)
  ggsave(paste0(wd,"/save/Str/","Str_mod",".png"), width = 12, height = 12, plot = p_raincloud, dpi = 300)
}
####------------------------------------------------------------------####
####Alg Vis####
#MO/LO females
{subset_vis <- filter(Alg_temp, sex == "f")
p_raincloud <- ggplot(data = subset_vis, aes(x = interaction(type,genus), y = DegA, fill = mod))+
  
  geom_flat_violin( width = 1.25, position=position_nudge(x=0.15, y=0), alpha = 0.5) +
  
  geom_boxplot(color = "white", outlier.shape=NA, alpha=1, width=0.075, position=position_nudge(x=0.15, y=0)) +
  
  geom_point(shape = 16, color = "white", position=position_jitter(width=.075, height=0), size=0.4) +
  
  scale_fill_manual(values = c("n" = "blue", "y" = "orange"))+
  
  theme_classic() +
  
  theme(axis.title.x=element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
        axis.title.y=element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
        panel.background = element_rect(fill = "black"),
        plot.background = element_rect(fill = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(color = "white"),
        axis.text = element_text(color = "white", size = 30),
        axis.ticks = element_line(color = "white"),
        legend.background = element_rect(fill = "black"),
        legend.text = element_text(color = "white", size = 20),
        legend.title = element_text(color = "white", size = 20),
        plot.title = element_text(color = "white", hjust = 0.5),
        plot.subtitle = element_text(color = "white", hjust = 0.5),
        plot.caption = element_text(color = "white"),
        axis.line = element_line(color = "white")
  ) +
  
  guides(fill = guide_legend(title = "Modified ocelli")) + 
  
  scale_x_discrete(expand = expansion(mult = c(0.04, 0.0))) +
  
  ylab("mean deviation in °" ) +
  
  xlab("") +
  
  coord_flip() 
print(p_raincloud)
ggsave(paste0(wd,"/save/Alg/","Alg_MOuLO_females",".png"), width = 12, height = 12, plot = p_raincloud, dpi = 300)
}
#MO/LO males
{subset_vis <- filter(Alg_temp, sex == "m")
  p_raincloud <- ggplot(data = subset_vis, aes(x = interaction(type,genus), y = DegA, fill = mod))+
    
    geom_flat_violin( width = 1.25, position=position_nudge(x=0.15, y=0), alpha = 0.5) +
    
    geom_boxplot(color = "white", outlier.shape=NA, alpha=1, width=0.075, position=position_nudge(x=0.15, y=0)) +
    
    geom_point(shape = 16, color = "white", position=position_jitter(width=.075, height=0), size=0.4) +
    
    scale_fill_manual(values = c("n" = "blue", "y" = "orange"))+
    
    theme_classic() +
    
    theme(axis.title.x=element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
          axis.title.y=element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
          panel.background = element_rect(fill = "black"),
          plot.background = element_rect(fill = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(color = "white"),
          axis.text = element_text(color = "white", size = 30),
          axis.ticks = element_line(color = "white"),
          legend.background = element_rect(fill = "black"),
          legend.text = element_text(color = "white", size = 20),
          legend.title = element_text(color = "white", size = 20),
          plot.title = element_text(color = "white", hjust = 0.5),
          plot.subtitle = element_text(color = "white", hjust = 0.5),
          plot.caption = element_text(color = "white"),
          axis.line = element_line(color = "white")
    ) +
    
    guides(fill = guide_legend(title = "Modified ocelli")) + 
    
    scale_x_discrete(expand = expansion(mult = c(0.04, 0.0))) +
    
    ylab("mean deviation in °" ) +
    
    xlab("") +
    
    coord_flip() 
  print(p_raincloud)
  ggsave(paste0(wd,"/save/Alg/","Alg_MOuLO_females",".png"), width = 12, height = 12, plot = p_raincloud, dpi = 300)
}
#MO females
{subset_vis <- filter(Alg_temp, (type == "MO"| type == "hMO" | type == "cMO") & sex == "f")
  p_raincloud <- ggplot(data = subset_vis, aes(x = genus, y = DegA, fill = mod))+
    
    geom_flat_violin( width = 1.25, position=position_nudge(x=0.15, y=0), alpha = 0.5) +
    
    geom_boxplot(color = "white", outlier.shape=NA, alpha=1, width=0.075, position=position_nudge(x=0.15, y=0)) +
    
    geom_point(shape = 16, color = "white", position=position_jitter(width=.075, height=0), size=0.4) +
    
    scale_fill_manual(values = c("n" = "blue", "y" = "orange"))+
    
    theme_classic() +
    
    theme(axis.title.x=element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
          axis.title.y=element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
          panel.background = element_rect(fill = "black"),
          plot.background = element_rect(fill = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(color = "white"),
          axis.text = element_text(color = "white", size = 30),
          axis.ticks = element_line(color = "white"),
          legend.background = element_rect(fill = "black"),
          legend.text = element_text(color = "white", size = 20),
          legend.title = element_text(color = "white", size = 20),
          plot.title = element_text(color = "white", hjust = 0.5),
          plot.subtitle = element_text(color = "white", hjust = 0.5),
          plot.caption = element_text(color = "white"),
          axis.line = element_line(color = "white")
    ) +
    
    guides(fill = guide_legend(title = "Modified ocelli")) + 
    
    scale_x_discrete(expand = expansion(mult = c(0.04, 0.0))) +
    
    ylab("mean deviation in °" ) +
    
    xlab("") +
    
    coord_flip() 
  print(p_raincloud)
  ggsave(paste0(wd,"/save/Alg/","Alg_MO_females",".png"), width = 12, height = 12, plot = p_raincloud, dpi = 300)
}
#MO males
{subset_vis <- filter(Alg_temp, (type == "MO"| type == "hMO" | type == "cMO") & sex == "m")
  p_raincloud <- ggplot(data = subset_vis, aes(x = genus, y = DegA, fill = mod))+
    
    geom_flat_violin( width = 1.25, position=position_nudge(x=0.15, y=0), alpha = 0.5) +
    
    geom_boxplot(color = "white", outlier.shape=NA, alpha=1, width=0.075, position=position_nudge(x=0.15, y=0)) +
    
    geom_point(shape = 16, color = "white", position=position_jitter(width=.075, height=0), size=0.4) +
    
    scale_fill_manual(values = c("n" = "blue", "y" = "orange"))+
    
    theme_classic() +
    
    theme(axis.title.x=element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
          axis.title.y=element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
          panel.background = element_rect(fill = "black"),
          plot.background = element_rect(fill = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(color = "white"),
          axis.text = element_text(color = "white", size = 30),
          axis.ticks = element_line(color = "white"),
          legend.background = element_rect(fill = "black"),
          legend.text = element_text(color = "white", size = 20),
          legend.title = element_text(color = "white", size = 20),
          plot.title = element_text(color = "white", hjust = 0.5),
          plot.subtitle = element_text(color = "white", hjust = 0.5),
          plot.caption = element_text(color = "white"),
          axis.line = element_line(color = "white")
    ) +
    
    guides(fill = guide_legend(title = "Modified ocelli")) + 
    
    scale_x_discrete(expand = expansion(mult = c(0.04, 0.0))) +
    
    ylab("mean deviation in °") +
    
    xlab("") +
    
    coord_flip() 
  print(p_raincloud)
  ggsave(paste0(wd,"/save/Alg/","Alg_MO_males",".png"), width = 12, height = 12, plot = p_raincloud, dpi = 300)
}
#LO females
{subset_vis <- filter(Alg_temp, (type == "LO") & sex == "f")
  p_raincloud <- ggplot(data = subset_vis, aes(x = genus, y = DegA, fill = mod))+
    
    geom_flat_violin( width = 1.25, position=position_nudge(x=0.15, y=0), alpha = 0.5) +
    
    geom_boxplot(color = "white", outlier.shape=NA, alpha=1, width=0.075, position=position_nudge(x=0.15, y=0)) +
    
    geom_point(shape = 16, color = "white", position=position_jitter(width=.075, height=0), size=0.4) +
    
    scale_fill_manual(values = c("n" = "blue", "y" = "orange"))+
    
    theme_classic() +
    
    theme(axis.title.x=element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
          axis.title.y=element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
          panel.background = element_rect(fill = "black"),
          plot.background = element_rect(fill = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(color = "white"),
          axis.text = element_text(color = "white", size = 30),
          axis.ticks = element_line(color = "white"),
          legend.background = element_rect(fill = "black"),
          legend.text = element_text(color = "white", size = 20),
          legend.title = element_text(color = "white", size = 20),
          plot.title = element_text(color = "white", hjust = 0.5),
          plot.subtitle = element_text(color = "white", hjust = 0.5),
          plot.caption = element_text(color = "white"),
          axis.line = element_line(color = "white")
    ) +
    
    guides(fill = guide_legend(title = "Modified ocelli")) + 
    
    scale_x_discrete(expand = expansion(mult = c(0.04, 0.0))) +
    
    ylab("mean deviation in °" ) +
    
    xlab("") +
    
    coord_flip() 
  print(p_raincloud)
  ggsave(paste0(wd,"/save/Alg/","Alg_LO_females",".png"), width = 12, height = 12, plot = p_raincloud, dpi = 300)
}
#LO males 
{subset_vis <- filter(Alg_temp, (type == "LO") & sex == "m")
  
  p_raincloud <- ggplot(data = subset_vis, aes(x = genus, y = DegA, fill = mod))+
    
    geom_flat_violin( width = 1.25, position=position_nudge(x=0.15, y=0), alpha = 0.5) +
    
    geom_boxplot(color = "white", outlier.shape=NA, alpha=1, width=0.075, position=position_nudge(x=0.15, y=0)) +
    
    geom_point(shape = 16, color = "white", position=position_jitter(width=.075, height=0), size=0.4) +
    
    scale_fill_manual(values = c("n" = "blue", "y" = "orange"))+
    
    theme_classic() +
    
    theme(axis.title.x=element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
          axis.title.y=element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
          panel.background = element_rect(fill = "black"),
          plot.background = element_rect(fill = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(color = "white"),
          axis.text = element_text(color = "white", size = 30),
          axis.ticks = element_line(color = "white"),
          legend.background = element_rect(fill = "black"),
          legend.text = element_text(color = "white", size = 20),
          legend.title = element_text(color = "white", size = 20),
          plot.title = element_text(color = "white", hjust = 0.5),
          plot.subtitle = element_text(color = "white", hjust = 0.5),
          plot.caption = element_text(color = "white"),
          axis.line = element_line(color = "white")
    ) +
    
    guides(fill = guide_legend(title = "Modified ocelli")) + 
    
    scale_x_discrete(expand = expansion(mult = c(0.04, 0.0))) +
    
    ylab("mean deviation in °") +
    
    xlab("") +
    
    coord_flip() 
  print(p_raincloud)
  ggsave(paste0(wd,"/save/Alg/","Alg_LO_males",".png"), width = 12, height = 12, plot = p_raincloud, dpi = 300)
}
####mod####
#MO/LO females
{subset_vis <- filter(Alg_temp, sex == "f")


p_raincloud <- ggplot(data = subset_vis, aes(x = mod, y = DegA, fill = mod))+
  
  geom_flat_violin( width = 1.6, position=position_nudge(x=0, y=0), alpha = 0.5) +
  
  geom_boxplot(color = "black", outlier.shape=NA, alpha=1, width=0.2, position=position_nudge(x=0, y=0)) +
  
  #geom_point(shape = 16, color = "black", position=position_jitter(width=.075, height=0), size=0.4) +
  
  scale_fill_manual(values = c("n" = "white", "y" = "#4F4F4F"))+
  
  theme_classic() +
  
  theme(axis.title.x=element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
        axis.title.y=element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        panel.grid.major.x =  element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black", size = 30),
        axis.ticks = element_line(color = "black"),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(color = "black", size = 20),
        legend.title = element_text(color = "black", size = 20),
        plot.title = element_text(color = "black", hjust = 0.5),
        plot.subtitle = element_text(color = "black", hjust = 0.5),
        plot.caption = element_text(color = "black"),
        axis.line = element_line(color = "black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()
  ) +
  
  guides(fill = guide_legend(title = "Modified ocelli")) + 
  
  stat_summary(fun = mean, geom = "crossbar", width = 0.6, color = "black", 
               linetype = "dashed", size = 0.2, position = position_nudge(x = 0, y = 0)) +
  
  ylab("rhabdom length in µm" ) +
  
  guides(fill = "none")+
  
  xlab("") +
  
  
  #scale_x_discrete(labels = labels_for_interaction, 
  #expand = expansion(mult = c(0.05, 0))) +
  scale_y_continuous(
    limits = c(0, NA),               # Behalte dein Limit von -5 und NA (automatische Obergrenze)
    breaks = seq(0, 15, by = 2.5),     # Behalte die Breaks von -5 bis 80 in 5er-Schritten
    #labels = function(x) ifelse(x %% 5 == 0 & x >= 0 & x <= 80, x, "")  # Zeige nur 0, 20, 40, 60, 80
  )+
  
  coord_flip() 
print(p_raincloud)
ggsave(paste0(wd,"/save/Alg/","Alg_mod",".png"), width = 12, height = 12, plot = p_raincloud, dpi = 300)
}
#MO/LO males
{subset_vis <- filter(Alg_temp, sex == "m")
  
  p_raincloud <- ggplot(data = subset_vis, aes(x = mod, y = DegA, fill = mod))+
    
    geom_flat_violin( width = 1.6, position=position_nudge(x=0, y=0), alpha = 0.5) +
    
    geom_boxplot(color = "black", outlier.shape=NA, alpha=1, width=0.2, position=position_nudge(x=0, y=0)) +
    
    #geom_point(shape = 16, color = "black", position=position_jitter(width=.075, height=0), size=0.4) +
    
    scale_fill_manual(values = c("n" = "white", "y" = "#4F4F4F"))+
    
    theme_classic() +
    
    theme(axis.title.x=element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
          axis.title.y=element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"),
          panel.grid.major.x =  element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.title = element_text(color = "black"),
          axis.text = element_text(color = "black", size = 30),
          axis.ticks = element_line(color = "black"),
          legend.background = element_rect(fill = "white"),
          legend.text = element_text(color = "black", size = 20),
          legend.title = element_text(color = "black", size = 20),
          plot.title = element_text(color = "black", hjust = 0.5),
          plot.subtitle = element_text(color = "black", hjust = 0.5),
          plot.caption = element_text(color = "black"),
          axis.line = element_line(color = "black"),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank()
    ) +
    
    guides(fill = guide_legend(title = "Modified ocelli")) + 
    
    stat_summary(fun = mean, geom = "crossbar", width = 0.6, color = "black", 
                 linetype = "dashed", size = 0.2, position = position_nudge(x = 0, y = 0)) +
    
    ylab("rhabdom length in µm" ) +
    
    guides(fill = "none")+
    
    xlab("") +
    
    
    #scale_x_discrete(labels = labels_for_interaction, 
    #expand = expansion(mult = c(0.05, 0))) +
    scale_y_continuous(
      limits = c(0, NA),               # Behalte dein Limit von -5 und NA (automatische Obergrenze)
      breaks = seq(0, 15, by = 2.5),     # Behalte die Breaks von -5 bis 80 in 5er-Schritten
      #labels = function(x) ifelse(x %% 5 == 0 & x >= 0 & x <= 80, x, "")  # Zeige nur 0, 20, 40, 60, 80
    )+
    coord_flip() 
  print(p_raincloud)
  ggsave(paste0(wd,"/save/Alg/","Alg_mod",".png"), width = 12, height = 12, plot = p_raincloud, dpi = 300)
}
####------------------------------------------------------------------####
####Black and white for Thesis####
####STR####
#MO/LO females
{subset_vis <- filter(Str_temp, sex == "f")

# Gewünschte Reihenfolge der Genera und des Typs setzen

genus_order <- c("Ampulex", "Sceliphron", "Ammophila", "Bembix", "Microbembex", "Liris", "Tachysphex", "Tachytes")
genus_order <- rev(genus_order)
subset_vis$genus <- factor(subset_vis$genus, levels = genus_order)
subset_vis$type <- factor(subset_vis$type, levels = c("MO", "LO"))

subset_vis$genus_type <- interaction( subset_vis$type,subset_vis$genus)

#genus_labels <- unique(subset_vis$genus)
#labels_for_interaction <- rep(genus_labels, each = 2)
#labels_for_interaction <- c("", labels_for_interaction[-length(labels_for_interaction)])
p_raincloud <- ggplot(data = subset_vis, aes(x = genus_type, y = DegS, fill = mod))+
  
  geom_flat_violin( width = 1.6, position=position_nudge(x=0, y=0), alpha = 0.5) +
  
  geom_boxplot(color = "black", outlier.shape=NA, alpha=1, width=0.2, position=position_nudge(x=0, y=0)) +
  
  #geom_point(shape = 16, color = "black", position=position_jitter(width=.075, height=0), size=0.4) +
  
  scale_fill_manual(values = c("n" = "white", "y" = "#4F4F4F"))+
  
  theme_classic() +
  
  theme(axis.title.x=element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
        axis.title.y=element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        panel.grid.major.x =  element_line(),
        panel.grid.minor.x = element_line(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black", size = 30),
        axis.ticks = element_line(color = "black"),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(color = "black", size = 20),
        legend.title = element_text(color = "black", size = 20),
        plot.title = element_text(color = "black", hjust = 0.5),
        plot.subtitle = element_text(color = "black", hjust = 0.5),
        plot.caption = element_text(color = "black"),
        axis.line = element_line(color = "black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()
  ) +
  
  guides(fill = guide_legend(title = "Modified ocelli")) + 
  
  stat_summary(fun = mean, geom = "crossbar", width = 0.6, color = "black", 
               linetype = "dashed", size = 0.2, position = position_nudge(x = 0, y = 0)) +
  
  ylab("mean straightness deviation in °" ) +
  
  guides(fill = "none")+
  
  xlab("") +
 
  
  #scale_x_discrete(labels = labels_for_interaction, 
                   #expand = expansion(mult = c(0.05, 0))) +
  scale_y_continuous(
    limits = c(-5, NA),               # Behalte dein Limit von -5 und NA (automatische Obergrenze)
    breaks = seq(0, 80, by = 5),     # Behalte die Breaks von -5 bis 80 in 5er-Schritten
    labels = function(x) ifelse(x %% 10 == 0 & x >= 0 & x <= 80, x, "")  # Zeige nur 0, 20, 40, 60, 80
  )+
  
  coord_flip() 
print(p_raincloud)
ggsave(paste0(wd,"/save/Str/","Str_MOuLO_femalesBW",".png"), width = 12, height = 12, plot = p_raincloud, dpi = 300)
}


#MO/LO males
{subset_vis <- filter(Str_temp, sex == "m")
  
  # Gewünschte Reihenfolge der Genera und des Typs setzen
  
  genus_order <- c("Ampulex", "Sceliphron", "Ammophila", "Bembix", "Bicyrtes", "Microbembex", "Stictiella", "Liris", "Tachysphex")
  genus_order <- rev(genus_order)
  subset_vis$genus <- factor(subset_vis$genus, levels = genus_order)
  subset_vis$type <- factor(subset_vis$type, levels = c("MO", "LO"))
  
  subset_vis$genus_type <- interaction( subset_vis$type,subset_vis$genus)
  
  #genus_labels <- unique(subset_vis$genus)
  #labels_for_interaction <- rep(genus_labels, each = 2)
  #labels_for_interaction <- c("", labels_for_interaction[-length(labels_for_interaction)])
  p_raincloud <- ggplot(data = subset_vis, aes(x = genus_type, y = DegS, fill = mod))+
    
    geom_flat_violin( width = 1.6, position=position_nudge(x=0, y=0), alpha = 0.5) +
    
    geom_boxplot(color = "black", outlier.shape=NA, alpha=1, width=0.2, position=position_nudge(x=0, y=0)) +
    
    #geom_point(shape = 16, color = "black", position=position_jitter(width=.075, height=0), size=0.4) +
    
    scale_fill_manual(values = c("n" = "white", "y" = "#4F4F4F"))+
    
    theme_classic() +
    
    theme(axis.title.x=element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
          axis.title.y=element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"),
          panel.grid.major.x =  element_line(),
          panel.grid.minor.x = element_line(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.title = element_text(color = "black"),
          axis.text = element_text(color = "black", size = 30),
          axis.ticks = element_line(color = "black"),
          legend.background = element_rect(fill = "white"),
          legend.text = element_text(color = "black", size = 20),
          legend.title = element_text(color = "black", size = 20),
          plot.title = element_text(color = "black", hjust = 0.5),
          plot.subtitle = element_text(color = "black", hjust = 0.5),
          plot.caption = element_text(color = "black"),
          axis.line = element_line(color = "black"),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank()
    ) +
    
    guides(fill = guide_legend(title = "Modified ocelli")) + 
    
    stat_summary(fun = mean, geom = "crossbar", width = 0.6, color = "black", 
                 linetype = "dashed", size = 0.2, position = position_nudge(x = 0, y = 0)) +
    
    ylab("mean straightness deviation in °" ) +
    
    guides(fill = "none")+
    
    xlab("") +
    
    
    #scale_x_discrete(labels = labels_for_interaction, 
    #expand = expansion(mult = c(0.05, 0))) +
    scale_y_continuous(
      limits = c(-5, NA),               # Behalte dein Limit von -5 und NA (automatische Obergrenze)
      breaks = seq(0, 80, by = 5),     # Behalte die Breaks von -5 bis 80 in 5er-Schritten
      labels = function(x) ifelse(x %% 10 == 0 & x >= 0 & x <= 80, x, "")  # Zeige nur 0, 20, 40, 60, 80
    )+
    
    coord_flip() 
  print(p_raincloud)
  ggsave(paste0(wd,"/save/Str/","Str_MOuLO_malesBW",".png"), width = 12, height = 13.5, plot = p_raincloud, dpi = 300)
}
####ALG####
#MO/LO females
{subset_vis <- filter(Alg_temp, sex == "f")

# Gewünschte Reihenfolge der Genera und des Typs setzen

genus_order <- c("Ampulex", "Sceliphron", "Ammophila", "Bembix", "Microbembex", "Liris", "Tachysphex", "Tachytes")
genus_order <- rev(genus_order)
subset_vis$genus <- factor(subset_vis$genus, levels = genus_order)
subset_vis$type <- factor(subset_vis$type, levels = c("MO", "LO"))

subset_vis$genus_type <- interaction( subset_vis$type,subset_vis$genus)

#genus_labels <- unique(subset_vis$genus)
#labels_for_interaction <- rep(genus_labels, each = 2)
#labels_for_interaction <- c("", labels_for_interaction[-length(labels_for_interaction)])
p_raincloud <- ggplot(data = subset_vis, aes(x = genus_type, y = DegA, fill = mod))+
  
  geom_flat_violin( width = 1.6, position=position_nudge(x=0, y=0), alpha = 0.5) +
  
  geom_boxplot(color = "black", outlier.shape=NA, alpha=1, width=0.2, position=position_nudge(x=0, y=0)) +
  
  #geom_point(shape = 16, color = "black", position=position_jitter(width=.075, height=0), size=0.4) +
  
  scale_fill_manual(values = c("n" = "white", "y" = "#4F4F4F"))+
  
  theme_classic() +
  
  theme(axis.title.x=element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
        axis.title.y=element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        panel.grid.major.x =  element_line(),
        panel.grid.minor.x = element_line(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black", size = 30),
        axis.ticks = element_line(color = "black"),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(color = "black", size = 20),
        legend.title = element_text(color = "black", size = 20),
        plot.title = element_text(color = "black", hjust = 0.5),
        plot.subtitle = element_text(color = "black", hjust = 0.5),
        plot.caption = element_text(color = "black"),
        axis.line = element_line(color = "black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()
  ) +
  
  guides(fill = guide_legend(title = "Modified ocelli")) + 
  
  stat_summary(fun = mean, geom = "crossbar", width = 0.6, color = "black", 
               linetype = "dashed", size = 0.2, position = position_nudge(x = 0, y = 0)) +
  
  ylab("mean alignment deviation in °" ) +
  
  guides(fill = "none")+
  
  xlab("") +
  
  
  #scale_x_discrete(labels = labels_for_interaction, 
  #expand = expansion(mult = c(0.05, 0))) +
  scale_y_continuous(
    limits = c(-10, NA),               # Behalte dein Limit von -5 und NA (automatische Obergrenze)
    breaks = seq(0, 80, by = 10),     # Behalte die Breaks von -5 bis 80 in 5er-Schritten
    labels = function(x) ifelse(x %% 20 == 0 & x >= 0 & x <= 80, x, "")  # Zeige nur 0, 20, 40, 60, 80
  )+
  
  coord_flip() 
print(p_raincloud)
ggsave(paste0(wd,"/save/Alg/","Alg_MOuLO_femalesBW",".png"), width = 12, height = 12, plot = p_raincloud, dpi = 300)
}
#MO/LO males
{subset_vis <- filter(Alg_temp, sex == "m")
  
  # Gewünschte Reihenfolge der Genera und des Typs setzen
  
  genus_order <- c("Ampulex", "Sceliphron", "Ammophila", "Bembix", "Bicyrtes", "Microbembex", "Stictiella", "Liris", "Tachysphex")
  genus_order <- rev(genus_order)
  subset_vis$genus <- factor(subset_vis$genus, levels = genus_order)
  subset_vis$type <- factor(subset_vis$type, levels = c("MO", "LO"))
  
  subset_vis$genus_type <- interaction( subset_vis$type,subset_vis$genus)
  
  #genus_labels <- unique(subset_vis$genus)
  #labels_for_interaction <- rep(genus_labels, each = 2)
  #labels_for_interaction <- c("", labels_for_interaction[-length(labels_for_interaction)])
  p_raincloud <- ggplot(data = subset_vis, aes(x = genus_type, y = DegA, fill = mod))+
    
    geom_flat_violin( width = 1.6, position=position_nudge(x=0, y=0), alpha = 0.5) +
    
    geom_boxplot(color = "black", outlier.shape=NA, alpha=1, width=0.2, position=position_nudge(x=0, y=0)) +
    
    #geom_point(shape = 16, color = "black", position=position_jitter(width=.075, height=0), size=0.4) +
    
    scale_fill_manual(values = c("n" = "white", "y" = "#4F4F4F"))+
    
    theme_classic() +
    
    theme(axis.title.x=element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
          axis.title.y=element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"),
          panel.grid.major.x =  element_line(),
          panel.grid.minor.x = element_line(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.title = element_text(color = "black"),
          axis.text = element_text(color = "black", size = 30),
          axis.ticks = element_line(color = "black"),
          legend.background = element_rect(fill = "white"),
          legend.text = element_text(color = "black", size = 20),
          legend.title = element_text(color = "black", size = 20),
          plot.title = element_text(color = "black", hjust = 0.5),
          plot.subtitle = element_text(color = "black", hjust = 0.5),
          plot.caption = element_text(color = "black"),
          axis.line = element_line(color = "black"),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank()
    ) +
    
    guides(fill = guide_legend(title = "Modified ocelli")) + 
    
    stat_summary(fun = mean, geom = "crossbar", width = 0.6, color = "black", 
                 linetype = "dashed", size = 0.2, position = position_nudge(x = 0, y = 0)) +
    
    ylab("mean alignment deviation in °" ) +
    
    guides(fill = "none")+
    
    xlab("") +
    
    
    #scale_x_discrete(labels = labels_for_interaction, 
    #expand = expansion(mult = c(0.05, 0))) +
    scale_y_continuous(
      limits = c(-10, NA),               # Behalte dein Limit von -5 und NA (automatische Obergrenze)
      breaks = seq(0, 80, by = 10),     # Behalte die Breaks von -5 bis 80 in 5er-Schritten
      labels = function(x) ifelse(x %% 20 == 0 & x >= 0 & x <= 80, x, "")  # Zeige nur 0, 20, 40, 60, 80
    )+
    coord_flip() 
  print(p_raincloud)
  ggsave(paste0(wd,"/save/Alg/","Alg_MOuLO_malesBW",".png"), width = 12, height = 13.5, plot = p_raincloud, dpi = 300)
}
####LEGTH####
#MO/LO females
{subset_vis <- filter(Str_temp, sex == "f")

# Gewünschte Reihenfolge der Genera und des Typs setzen

genus_order <- c("Ampulex", "Sceliphron", "Ammophila", "Bembix", "Microbembex", "Liris", "Tachysphex", "Tachytes")
genus_order <- rev(genus_order)
subset_vis$genus <- factor(subset_vis$genus, levels = genus_order)
subset_vis$type <- factor(subset_vis$type, levels = c("MO", "LO"))

subset_vis$genus_type <- interaction( subset_vis$type,subset_vis$genus)

#genus_labels <- unique(subset_vis$genus)
#labels_for_interaction <- rep(genus_labels, each = 2)
#labels_for_interaction <- c("", labels_for_interaction[-length(labels_for_interaction)])
p_raincloud <- ggplot(data = subset_vis, aes(x = genus_type, y = rhlength, fill = mod))+
  
  geom_flat_violin( width = 1.6, position=position_nudge(x=0, y=0), alpha = 0.5) +
  
  geom_boxplot(color = "black", outlier.shape=NA, alpha=1, width=0.2, position=position_nudge(x=0, y=0)) +
  
  #geom_point(shape = 16, color = "black", position=position_jitter(width=.075, height=0), size=0.4) +
  
  scale_fill_manual(values = c("n" = "white", "y" = "#4F4F4F"))+
  
  theme_classic() +
  
  theme(axis.title.x=element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
        axis.title.y=element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        panel.grid.major.x =  element_line(),
        panel.grid.minor.x = element_line(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black", size = 30),
        axis.ticks = element_line(color = "black"),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(color = "black", size = 20),
        legend.title = element_text(color = "black", size = 20),
        plot.title = element_text(color = "black", hjust = 0.5),
        plot.subtitle = element_text(color = "black", hjust = 0.5),
        plot.caption = element_text(color = "black"),
        axis.line = element_line(color = "black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()
  ) +
  
  guides(fill = guide_legend(title = "Modified ocelli")) + 
  
  stat_summary(fun = mean, geom = "crossbar", width = 0.6, color = "black", 
               linetype = "dashed", size = 0.2, position = position_nudge(x = 0, y = 0)) +
  
  ylab("rhabdom length in µm" ) +
  
  guides(fill = "none")+
  
  xlab("") +
  
  
  #scale_x_discrete(labels = labels_for_interaction, 
  #expand = expansion(mult = c(0.05, 0))) +
  scale_y_continuous(
    limits = c(0, NA),               # Behalte dein Limit von -5 und NA (automatische Obergrenze)
    breaks = seq(0, 15, by = 2.5),     # Behalte die Breaks von -5 bis 80 in 5er-Schritten
    #labels = function(x) ifelse(x %% 5 == 0 & x >= 0 & x <= 80, x, "")  # Zeige nur 0, 20, 40, 60, 80
  )+
  
  coord_flip() 
print(p_raincloud)
ggsave(paste0(wd,"/save/Str/","Length_MOuLO_femalesBW",".png"), width = 12, height = 12, plot = p_raincloud, dpi = 300)
}
#MO/LO males
{subset_vis <- filter(Str_temp, sex == "m")
  
  # Gewünschte Reihenfolge der Genera und des Typs setzen
  
  genus_order <- c("Ampulex", "Sceliphron", "Ammophila", "Bembix", "Bicyrtes", "Microbembex", "Stictiella", "Liris", "Tachysphex")
  genus_order <- rev(genus_order)
  subset_vis$genus <- factor(subset_vis$genus, levels = genus_order)
  subset_vis$type <- factor(subset_vis$type, levels = c("MO", "LO"))
  
  subset_vis$genus_type <- interaction( subset_vis$type,subset_vis$genus)
  
  #genus_labels <- unique(subset_vis$genus)
  #labels_for_interaction <- rep(genus_labels, each = 2)
  #labels_for_interaction <- c("", labels_for_interaction[-length(labels_for_interaction)])
  p_raincloud <- ggplot(data = subset_vis, aes(x = genus_type, y = rhlength, fill = mod))+
    
    geom_flat_violin( width = 1.6, position=position_nudge(x=0, y=0), alpha = 0.5) +
    
    geom_boxplot(color = "black", outlier.shape=NA, alpha=1, width=0.2, position=position_nudge(x=0, y=0)) +
    
    #geom_point(shape = 16, color = "black", position=position_jitter(width=.075, height=0), size=0.4) +
    
    scale_fill_manual(values = c("n" = "white", "y" = "#4F4F4F"))+
    
    theme_classic() +
    
    theme(axis.title.x=element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
          axis.title.y=element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"),
          panel.grid.major.x =  element_line(),
          panel.grid.minor.x = element_line(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.title = element_text(color = "black"),
          axis.text = element_text(color = "black", size = 30),
          axis.ticks = element_line(color = "black"),
          legend.background = element_rect(fill = "white"),
          legend.text = element_text(color = "black", size = 20),
          legend.title = element_text(color = "black", size = 20),
          plot.title = element_text(color = "black", hjust = 0.5),
          plot.subtitle = element_text(color = "black", hjust = 0.5),
          plot.caption = element_text(color = "black"),
          axis.line = element_line(color = "black"),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank()
    ) +
    
    guides(fill = guide_legend(title = "Modified ocelli")) + 
    
    stat_summary(fun = mean, geom = "crossbar", width = 0.6, color = "black", 
                 linetype = "dashed", size = 0.2, position = position_nudge(x = 0, y = 0)) +
    
    ylab("rhabdom length in µm" ) +
    
    guides(fill = "none")+
    
    xlab("") +
    
    
    #scale_x_discrete(labels = labels_for_interaction, 
    #expand = expansion(mult = c(0.05, 0))) +
    scale_y_continuous(
      limits = c(0, NA),               # Behalte dein Limit von -5 und NA (automatische Obergrenze)
      breaks = seq(0, 15, by = 2.5),     # Behalte die Breaks von -5 bis 80 in 5er-Schritten
      #labels = function(x) ifelse(x %% 5 == 0 & x >= 0 & x <= 80, x, "")  # Zeige nur 0, 20, 40, 60, 80
    )+
    coord_flip() 
  print(p_raincloud)
  ggsave(paste0(wd,"/save/Str/","Length_MOuLO_malesBW",".png"), width = 12, height = 12, plot = p_raincloud, dpi = 300)
}
####Length - mod####
#MO/LO females
{subset_vis <- filter(Str_temp, sex == "f")


p_raincloud <- ggplot(data = subset_vis, aes(x = mod, y = rhlength, fill = mod))+
  
  geom_flat_violin( width = 1.6, position=position_nudge(x=0, y=0), alpha = 0.5) +
  
  geom_boxplot(color = "black", outlier.shape=NA, alpha=1, width=0.2, position=position_nudge(x=0, y=0)) +
  
  #geom_point(shape = 16, color = "black", position=position_jitter(width=.075, height=0), size=0.4) +
  
  scale_fill_manual(values = c("n" = "white", "y" = "#4F4F4F"))+
  
  theme_classic() +
  
  theme(axis.title.x=element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
        axis.title.y=element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        panel.grid.major.x =  element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black", size = 30),
        axis.ticks = element_line(color = "black"),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(color = "black", size = 20),
        legend.title = element_text(color = "black", size = 20),
        plot.title = element_text(color = "black", hjust = 0.5),
        plot.subtitle = element_text(color = "black", hjust = 0.5),
        plot.caption = element_text(color = "black"),
        axis.line = element_line(color = "black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()
  ) +
  
  guides(fill = guide_legend(title = "Modified ocelli")) + 
  
  stat_summary(fun = mean, geom = "crossbar", width = 0.6, color = "black", 
               linetype = "dashed", size = 0.2, position = position_nudge(x = 0, y = 0)) +
  
  ylab("rhabdom length in µm - females" ) +
  
  guides(fill = "none")+
  
  xlab("") +
  
  
  #scale_x_discrete(labels = labels_for_interaction, 
  #expand = expansion(mult = c(0.05, 0))) +
  scale_y_continuous(
    limits = c(0, NA),               # Behalte dein Limit von -5 und NA (automatische Obergrenze)
    breaks = seq(0, 15, by = 2.5),     # Behalte die Breaks von -5 bis 80 in 5er-Schritten
    #labels = function(x) ifelse(x %% 5 == 0 & x >= 0 & x <= 80, x, "")  # Zeige nur 0, 20, 40, 60, 80
  )+
  
  coord_flip() 
print(p_raincloud)
ggsave(paste0(wd,"/save/Str/","Length_mod_femalesBW",".png"), width = 12, height = 12, plot = p_raincloud, dpi = 300)
}
#MO/LO males
{subset_vis <- filter(Str_temp, sex == "m")
  
  p_raincloud <- ggplot(data = subset_vis, aes(x = mod, y = rhlength, fill = mod))+
    
    geom_flat_violin( width = 1.6, position=position_nudge(x=0, y=0), alpha = 0.5) +
    
    geom_boxplot(color = "black", outlier.shape=NA, alpha=1, width=0.2, position=position_nudge(x=0, y=0)) +
    
    #geom_point(shape = 16, color = "black", position=position_jitter(width=.075, height=0), size=0.4) +
    
    scale_fill_manual(values = c("n" = "white", "y" = "#4F4F4F"))+
    
    theme_classic() +
    
    theme(axis.title.x=element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
          axis.title.y=element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"),
          panel.grid.major.x =  element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.title = element_text(color = "black"),
          axis.text = element_text(color = "black", size = 30),
          axis.ticks = element_line(color = "black"),
          legend.background = element_rect(fill = "white"),
          legend.text = element_text(color = "black", size = 20),
          legend.title = element_text(color = "black", size = 20),
          plot.title = element_text(color = "black", hjust = 0.5),
          plot.subtitle = element_text(color = "black", hjust = 0.5),
          plot.caption = element_text(color = "black"),
          axis.line = element_line(color = "black"),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank()
    ) +
    
    guides(fill = guide_legend(title = "Modified ocelli")) + 
    
    stat_summary(fun = mean, geom = "crossbar", width = 0.6, color = "black", 
                 linetype = "dashed", size = 0.2, position = position_nudge(x = 0, y = 0)) +
    
    ylab("rhabdom length in µm - males" ) +
    
    guides(fill = "none")+
    
    xlab("") +
    
    
    #scale_x_discrete(labels = labels_for_interaction, 
    #expand = expansion(mult = c(0.05, 0))) +
    scale_y_continuous(
      limits = c(0, NA),               # Behalte dein Limit von -5 und NA (automatische Obergrenze)
      breaks = seq(0, 15, by = 2.5),     # Behalte die Breaks von -5 bis 80 in 5er-Schritten
      #labels = function(x) ifelse(x %% 5 == 0 & x >= 0 & x <= 80, x, "")  # Zeige nur 0, 20, 40, 60, 80
    )+
    coord_flip() 
  print(p_raincloud)
  ggsave(paste0(wd,"/save/Str/","Length_mod_malesBW",".png"), width = 12, height = 12, plot = p_raincloud, dpi = 300)
}
####Str_mod####
#MO/LO females
{subset_vis <- filter(Str_temp, sex == "f")


p_raincloud <- ggplot(data = subset_vis, aes(x = mod, y = DegS, fill = mod))+
  
  geom_flat_violin( width = 1.6, position=position_nudge(x=0, y=0), alpha = 0.5) +
  
  geom_boxplot(color = "black", outlier.shape=NA, alpha=1, width=0.2, position=position_nudge(x=0, y=0)) +
  
  #geom_point(shape = 16, color = "black", position=position_jitter(width=.075, height=0), size=0.4) +
  
  scale_fill_manual(values = c("n" = "white", "y" = "#4F4F4F"))+
  
  theme_classic() +
  
  theme(axis.title.x=element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
        axis.title.y=element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        panel.grid.major.x =  element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black", size = 30),
        axis.ticks = element_line(color = "black"),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(color = "black", size = 20),
        legend.title = element_text(color = "black", size = 20),
        plot.title = element_text(color = "black", hjust = 0.5),
        plot.subtitle = element_text(color = "black", hjust = 0.5),
        plot.caption = element_text(color = "black"),
        axis.line = element_line(color = "black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()
  ) +
  
  guides(fill = guide_legend(title = "Modified ocelli")) + 
  
  stat_summary(fun = mean, geom = "crossbar", width = 0.6, color = "black", 
               linetype = "dashed", size = 0.2, position = position_nudge(x = 0, y = 0)) +
  
  ylab("rhabdom straightness deviation - females" ) +
  
  guides(fill = "none")+
  
  xlab("") +
  
  
  #scale_x_discrete(labels = labels_for_interaction, 
  #expand = expansion(mult = c(0.05, 0))) +
  scale_y_continuous(
    limits = c(-5, NA),               # Behalte dein Limit von -5 und NA (automatische Obergrenze)
    breaks = seq(0, 80, by = 5),     # Behalte die Breaks von -5 bis 80 in 5er-Schritten
    labels = function(x) ifelse(x %% 10 == 0 & x >= 0 & x <= 80, x, "")  # Zeige nur 0, 20, 40, 60, 80
  )+
  
  coord_flip() 
print(p_raincloud)
ggsave(paste0(wd,"/save/Str/","Str_mod_f",".png"), width = 12, height = 12, plot = p_raincloud, dpi = 300)
}
#MO/LO males
{subset_vis <- filter(Str_temp, sex == "m")
  
  p_raincloud <- ggplot(data = subset_vis, aes(x = mod, y = DegS, fill = mod))+
    
    geom_flat_violin( width = 1.6, position=position_nudge(x=0, y=0), alpha = 0.5) +
    
    geom_boxplot(color = "black", outlier.shape=NA, alpha=1, width=0.2, position=position_nudge(x=0, y=0)) +
    
    #geom_point(shape = 16, color = "black", position=position_jitter(width=.075, height=0), size=0.4) +
    
    scale_fill_manual(values = c("n" = "white", "y" = "#4F4F4F"))+
    
    theme_classic() +
    
    theme(axis.title.x=element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
          axis.title.y=element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"),
          panel.grid.major.x =  element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.title = element_text(color = "black"),
          axis.text = element_text(color = "black", size = 30),
          axis.ticks = element_line(color = "black"),
          legend.background = element_rect(fill = "white"),
          legend.text = element_text(color = "black", size = 20),
          legend.title = element_text(color = "black", size = 20),
          plot.title = element_text(color = "black", hjust = 0.5),
          plot.subtitle = element_text(color = "black", hjust = 0.5),
          plot.caption = element_text(color = "black"),
          axis.line = element_line(color = "black"),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank()
    ) +
    
    guides(fill = guide_legend(title = "Modified ocelli")) + 
    
    stat_summary(fun = mean, geom = "crossbar", width = 0.6, color = "black", 
                 linetype = "dashed", size = 0.2, position = position_nudge(x = 0, y = 0)) +
    
    ylab("rhabdom straightness deviation - males" ) +
    
    guides(fill = "none")+
    
    xlab("") +
    
    
    #scale_x_discrete(labels = labels_for_interaction, 
    #expand = expansion(mult = c(0.05, 0))) +
    scale_y_continuous(
      limits = c(-5, NA),               # Behalte dein Limit von -5 und NA (automatische Obergrenze)
      breaks = seq(0, 80, by = 5),     # Behalte die Breaks von -5 bis 80 in 5er-Schritten
      labels = function(x) ifelse(x %% 10 == 0 & x >= 0 & x <= 80, x, "")  # Zeige nur 0, 20, 40, 60, 80
    )+
    coord_flip() 
  print(p_raincloud)
  ggsave(paste0(wd,"/save/Str/","Str_mod_m",".png"), width = 12, height = 12, plot = p_raincloud, dpi = 300)
}
####Alg mod####
#MO/LO females
{subset_vis <- filter(Alg_temp, sex == "f")


p_raincloud <- ggplot(data = subset_vis, aes(x = mod, y = DegA, fill = mod))+
  
  geom_flat_violin( width = 1.6, position=position_nudge(x=0, y=0), alpha = 0.5) +
  
  geom_boxplot(color = "black", outlier.shape=NA, alpha=1, width=0.2, position=position_nudge(x=0, y=0)) +
  
  #geom_point(shape = 16, color = "black", position=position_jitter(width=.075, height=0), size=0.4) +
  
  scale_fill_manual(values = c("n" = "white", "y" = "#4F4F4F"))+
  
  theme_classic() +
  
  theme(axis.title.x=element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
        axis.title.y=element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        panel.grid.major.x =  element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black", size = 30),
        axis.ticks = element_line(color = "black"),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(color = "black", size = 20),
        legend.title = element_text(color = "black", size = 20),
        plot.title = element_text(color = "black", hjust = 0.5),
        plot.subtitle = element_text(color = "black", hjust = 0.5),
        plot.caption = element_text(color = "black"),
        axis.line = element_line(color = "black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()
  ) +
  
  guides(fill = guide_legend(title = "Modified ocelli")) + 
  
  stat_summary(fun = mean, geom = "crossbar", width = 0.6, color = "black", 
               linetype = "dashed", size = 0.2, position = position_nudge(x = 0, y = 0)) +
  
  ylab("rhabdom alignment deviation - females" ) +
  
  guides(fill = "none")+
  
  xlab("") +
  
  
  #scale_x_discrete(labels = labels_for_interaction, 
  #expand = expansion(mult = c(0.05, 0))) +
  scale_y_continuous(
    limits = c(-10, NA),               # Behalte dein Limit von -5 und NA (automatische Obergrenze)
    breaks = seq(0, 80, by = 10),     # Behalte die Breaks von -5 bis 80 in 5er-Schritten
    labels = function(x) ifelse(x %% 20 == 0 & x >= 0 & x <= 80, x, "")  # Zeige nur 0, 20, 40, 60, 80
  )+
  
  coord_flip() 
print(p_raincloud)
ggsave(paste0(wd,"/save/Alg/","Alg_mod_f",".png"), width = 12, height = 12, plot = p_raincloud, dpi = 300)
}
#MO/LO males
{subset_vis <- filter(Alg_temp, sex == "m")
  
  p_raincloud <- ggplot(data = subset_vis, aes(x = mod, y = DegA, fill = mod))+
    
    geom_flat_violin( width = 1.6, position=position_nudge(x=0, y=0), alpha = 0.5) +
    
    geom_boxplot(color = "black", outlier.shape=NA, alpha=1, width=0.2, position=position_nudge(x=0, y=0)) +
    
    #geom_point(shape = 16, color = "black", position=position_jitter(width=.075, height=0), size=0.4) +
    
    scale_fill_manual(values = c("n" = "white", "y" = "#4F4F4F"))+
    
    theme_classic() +
    
    theme(axis.title.x=element_text(margin = margin(t = 25, r = 0, b = 0, l = 0), size = 30),
          axis.title.y=element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), size = 30),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"),
          panel.grid.major.x =  element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.title = element_text(color = "black"),
          axis.text = element_text(color = "black", size = 30),
          axis.ticks = element_line(color = "black"),
          legend.background = element_rect(fill = "white"),
          legend.text = element_text(color = "black", size = 20),
          legend.title = element_text(color = "black", size = 20),
          plot.title = element_text(color = "black", hjust = 0.5),
          plot.subtitle = element_text(color = "black", hjust = 0.5),
          plot.caption = element_text(color = "black"),
          axis.line = element_line(color = "black"),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank()
    ) +
    
    guides(fill = guide_legend(title = "Modified ocelli")) + 
    
    stat_summary(fun = mean, geom = "crossbar", width = 0.6, color = "black", 
                 linetype = "dashed", size = 0.2, position = position_nudge(x = 0, y = 0)) +
    
    ylab("rhabdom alignment deviation - males" ) +
    
    guides(fill = "none")+
    
    xlab("") +
    
    
    #scale_x_discrete(labels = labels_for_interaction, 
    #expand = expansion(mult = c(0.05, 0))) +
    scale_y_continuous(
      limits = c(-10, NA),               # Behalte dein Limit von -5 und NA (automatische Obergrenze)
      breaks = seq(0, 80, by = 10),     # Behalte die Breaks von -5 bis 80 in 5er-Schritten
      labels = function(x) ifelse(x %% 20 == 0 & x >= 0 & x <= 80, x, "")  # Zeige nur 0, 20, 40, 60, 80
    )+
    coord_flip() 
  print(p_raincloud)
  ggsave(paste0(wd,"/save/Alg/","Alg_mod_m",".png"), width = 12, height = 12, plot = p_raincloud, dpi = 300)
}

####------------------------------------------------------------------####
ggplot(data = )
#### Plot Str~Alg ####
comb <- cbind(Alg,Str)

comb <- comb[ , !duplicated(names(comb))] # remove double columns
comb <- comb %>% select(-c(nr,xmid,ymid,x1,x2,y1,y2)) # drop nr column with smal"n"
comb <- comb[, c("Nr", "ID", "genus","species","sex","type","mod", "DegA","DegS","rhlength", "complete")] #reorder columns

ggplot(data = comb, mapping = aes(x = DegA, y = DegS))+
  geom_point()

comb <- filter(comb, mod == "y")

ggplot(data = comb, mapping = aes(x = DegA, y = DegS))+
  geom_point()
