**####WORK IN PROGRESS####** <br>
**Rhabom straightness and alignment deviation measurements in frontal sections of Hymenoptera.** <br>
All code is written in R studio (version 2022.07.2, build 576) and R (version 4.4.1) <br> 
<br>

**Introduction** <br>
<br>
    These scripts were written as part of the Dissertation “Inner morphology, phylogenetic significance and evolution as well as possible functions of modified ocelli in apoid wasps”. It provides a data pipeline to analyse the straightness and alignment of rhabdoms from frontal histological sections of the retina. The aim of this analysis is to determine whether polarisation detection is likely in the specimen and allow for a statistical comparison. In short, polarisation detection capabilities depend on the rhabdoms being straight and their alignment organised. In most Hymenoptera this means that all rhabdoms are, with individual deviations, horizontally aligned. The traditional method to evaluate this is an artificial horizon in relation to which the deviation of each rhabdom can be measured. As an alternative, this pipeline uses rhabdom alignment deviation, for which the deviation of a set number of neighbours from the rhabdoms orientation is measured. This is necessary, as it ensures comparability across the highly irregular modified retinae of some apoid wasps. For details on the use case please refer to the aforementioned thesis. 

**How to use the scripts** <br>
<br>
The scripts “Rhabdom Straightness Script Final” (RSS) and “Rhabdom Alignment by Neighbour” (RAN) can easily be used for different projects. “Stats and Vis Final” as well as “lateral sections R script” are provided for completions sake and as part of good scientific practice during submission of the Dissertation. Their application is very narrow and specific to the thesis. Thus, this read.me will focus on explaining the first two scripts and how to use them. Details on the other two scripts can be found in the Dissertation. 
<br>
<br>
     RSS should be run first. It establishes a folder structure and RAN will require its output. RSS requires a .txt file of a specific format and the source picture of the retina in .tif format in the work directory. The .tif files should have a dpi of 400, or you have to adjust the dpi at the very top of the script through the dpi_user variable. This ensure that the aspect ratio of your pictures is automatically used for some ouput graphs. The text files should have the following format: <br>
<br>
The first line contains, without blanks but with commas as separators: <br>
<br>
    “Genus,Species,ID,ocellus_type,scalebar_size,sex,modified_y_n,complete” <br>
<br>

 See “Reference .txt for input” for an example. <br>
<br>

   Genus: enter genus name <br>
   Species: enter species name <br>
   ID: your project ID, could be a s simple as numbering your samples <br>
   ocellus_type: MO for “median ocellus” or LO for “lateral ocellus” (mostly used for statistics and visualisation as well as folder names, could be NA or any type of category) <br>
   scalebar_size: you will add coordinates for your scalebar downstream, this defines the actual distance they represent, e.g. 10 µm <br>
   sex: sex of your specimen <br>
   modified_y_n: enter “y” if your specimen has modified ocelli or “n” if it has unmodified ocelli (this is also mostly used in the statistics and visualisation script, could be NA or any type of category) <br>
   complete: enter “c” if you used a picture of a complete retina, or “h” if you used a picture of half a retina (this is also mostly used in the statistics and visualisation script and specific to the usecase, could be NA or any type of category) <br>
<br>
    The next six line contain an x and y coordinate each. Line two to five should contain the coordinates of the four corners of your picture. Line six and seven contain the coordinates of two points on your scalebar, the real length of which you entered under “scalebar_size”. <br>
    The current version of this script is purpose written to deal with Image J files. Thus you only need to name the first column in the original header (now line 8) “time”. You can leave the rest of the data as is. For some files Image J creates an additional column called “Ch” which you can ignore as well. The name of your modified .txt files is irrelevant, the respective .tif should be called “IDpoints.tif” replacing “ID” with your specimen ID.  <br>
     Now you can run the script. <br>
**Outpout of RSS** <br>
<br>
**Principles of Analysis** <br>
<br>
    Image J was used to manually place five equidistant points along the middle of each rhabdom. Rhabdom straightness will be accessed by calculating the incline between point one and point five as the baseline using the arctan. Then the incline and, based on it, the  deviation of each segment (point one to two, point two to three and so on) in degree from the baseline will be calculated. The absolute values of these deviations will be used to calculate the mean deviation per rhabdom, ensuring that e.g. deviations of -45° and 45° do not cancel each other out. This principle is taken from Ribi, W., Warrant, E., & Zeil, J. (2011). The organization of honeybee ocelli: Regional specializations and rhabdom arrangements. Arthropod Structure and Development, 40(6), 509–520. https://doi.org/https://doi.org/10.1016/j.asd.2011.06.004. <br>
<br>
**Output of the scripts**
