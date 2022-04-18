##install the libraries:
.libPaths()
.libPaths('C:/R/library')
library("ape")
library(gdata)
library(phangorn)
library(ggplot2)
library(scales)
library(otuSummary)


####
setwd("C:/") #Set the way to the sequence allignment file

dna <- read.dna("file_name.fas", format="fasta") 
#load the sequence allignment file "file_name.fas" as a "dna" variable

dna2=dist.dna(dna, as.matrix = TRUE, model = "raw")
#create a distance matrix (rows - names of all viruses, collumns - names of all viruses)
#cell - pairwise distance for the corresponding viruses pair
#save the matrix as a "dna2" variable

dna3=lowerTriangle(dna2) 
#take matrix lower traingle (upper traingle is the same)
#display as three columns: id1,id2, distance
#save the matrix as a "dna3" variable


#the same is for the protein allignment:
protein=read.aa("file_name.fas", format = "fasta")
protein2=as.matrix(dist.ml(protein, model="Blosum62"))
protein3=lowerTriangle(protein2)

#create a plot:
dna_protein=
  ggplot(data.frame(dna3*100,protein3*100),aes(dna3*100,protein3*100))+
  geom_bin2d(bins=85)+
  scale_fill_gradientn(colours=c("blue","red"),trans = "log10")+ 
  theme(legend.justification=c(1,0))+
  labs(fill = "log10 count")+ 
  xlab('nucleotide distance, %')+ylab('protein distance, %')+
  scale_x_continuous(expand = c(0, 0),limits = c(0, max(dna3)*100+3)) + scale_y_continuous(expand = c(0, 0),limits = c(0, max(protein3)*100+3), breaks= pretty_breaks())+
  theme(
    panel.background = element_rect(fill = "grey90", colour = "grey100"),
    axis.line = element_line(colour = 'black', size = 0.5),
    axis.title.x = element_text( size = 12, angle = 0, hjust = .5, vjust = 1.5, face = "plain"),
    axis.title.y = element_text( size = 12, angle = 90, hjust = .5, vjust = 1.5, face = "plain"),
    axis.text.x = element_text(color = "black", size = 10,  hjust = .5, vjust = .5, face = "plain"),
    axis.text.y = element_text(color = "black", size = 10,  hjust = .5, vjust = .5, face = "plain"),
    axis.ticks.length=unit(.25, "cm"),
    plot.background=element_rect(fill="grey100")
  )+ 
  geom_vline(xintercept = 10, linetype="dashed",
             color = "black", size=2)
dna_protein

