
if(dev.cur() <= 1) get(getOption("device"))()

opar <-
    par(ask = interactive() &&
        (.Device %in% c("X11", "GTK", "windows","quartz"))
        )

#datadir <- system.file("examples", package = "DNAcopy")
#Read in two example by Snijders et al.

data(coriell)

#Combine into one CNA object to prepare for analysis on Chromosomes 1-23

CNA.object <- CNA(cbind(coriell$Coriell.05296,coriell$Coriell.13330),
                  coriell$Chromosome,coriell$Position,
                  data.type="logratio",sampleid=c("c05296","c13330"))

#We generally recommend smoothing single point outliers before analysis
#Make sure to check that the smoothing is proper

smoothed.CNA.object <- smooth.CNA(CNA.object)

#Segmentation at default parameters

segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1)

#Plot whole studies

plot(segment.smoothed.CNA.object, plot.type="w")

#Plot each study by chromosome

plot(segment.smoothed.CNA.object, plot.type="s")

#Plot each chromosome across studies (6 per page)

plot(segment.smoothed.CNA.object, plot.type="c", cbys.layout=c(2,1), cbys.nchrom=6)

#Plot by plateaus

plot(segment.smoothed.CNA.object, plot.type="p")

#Segment again but making sure that splits correspond are 3SDs separated

segment.smoothed.CNA.object <- segment(smoothed.CNA.object, undo.splits="sdundo", verbose=1)

#All the non-obvious splits have been removed

plot(segment.smoothed.CNA.object,plot.type="s")





