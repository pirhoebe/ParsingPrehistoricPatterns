###########################################################################
#### Modeltest #### 
###########################################################################
if(exists("myModel")==FALSE){
  myModel <- "exponential"
  message("exponential modeltest is chosen as default
change model by defining myModel")
}

###########################################################################
# perform modeltest
message("preparing model...")
message(paste("Chosen parameters:",
              "\n simulations : ", simulations,
              "\n running mean: ", runningmean,
              "\n time range  : ", lowerlim, " - ", upperlim,
              "\n model       : ", myModel,
              "\n redundant normalisation: ", normalised,
              "\n processor cores: ", ncores
              ))

set.mod <- modelTest(set.cal, 
                        errors=set.cal$metadata$Error, 
                        bins=bins,
                        nsim=simulations, 
                        runm=runningmean, 
                        timeRange=c(max((lowerlim)), 
                                    min((upperlim))), 
                        model=myModel, 
                        ncores=ncores, 
                        datenormalised=normalised)
message("complete!")