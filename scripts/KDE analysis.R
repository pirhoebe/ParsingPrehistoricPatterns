###############################################################################
#### Composite Kernel Density Estimate ####
###############################################################################
message(paste("Chosen parameters:",
              "\nsampleDates:",
              "\n simulations    : ", simulations,"\n",
              "\nckde: ",
              "\n timerange      : ", lowerlim, " - ", upperlim,
              "\n bandwidth      : ", bandwidth,
              "\n processor cores: ", ncores
))

message("sampling kernels...")
set.randates = sampleDates(set.cal,
                           bins=bins,
                           nsim=simulations,
                           verbose=FALSE,
                           boot = TRUE)
message("Kernel Density Estimate...")
set.ckde = ckde(set.randates,timeRange=c(lowerlim,upperlim),bw=bandwidth)
message("complete!")
