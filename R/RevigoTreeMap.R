# A treemap R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","freqInDbPercent","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0017113","dihydropyrimidine dehydrogenase (NADP+) activity",0.006,2.789,0.871,0.000,"dihydropyrimidine dehydrogenase (NADP+) activity"),
c("GO:0052851","ferric-chelate reductase (NADPH) activity",0.017,2.198,0.780,0.328,"dihydropyrimidine dehydrogenase (NADP+) activity"),
c("GO:0038048","dynorphin receptor activity",0.006,2.492,0.998,0.000,"dynorphin receptor activity"),
c("GO:0046911","metal chelating activity",0.006,2.920,0.999,0.000,"metal chelating activity"),
c("GO:0002058","uracil binding",0.011,2.466,0.773,0.001,"uracil binding"),
c("GO:0050661","NADP binding",0.293,2.258,0.875,0.230,"uracil binding"),
c("GO:0002054","nucleobase binding",0.028,2.088,0.886,0.245,"uracil binding"),
c("GO:0050682","AF-2 domain binding",0.017,2.120,0.999,0.001,"AF-2 domain binding"),
c("GO:0060228","phosphatidylcholine-sterol O-acyltransferase activator activity",0.033,2.143,0.998,0.012,"phosphatidylcholine-sterol O-acyltransferase activator activity"),
c("GO:0004098","cerebroside-sulfatase activity",0.006,2.773,0.956,0.079,"cerebroside-sulfatase activity"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = "REVIGO TreeMap",
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  fontsize.labels = c(0, 12),  # don't draw main group names as they are redundant
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  bg.labels = 0,
  position.legend = "none"
)

dev.off()

