###########################################################################
#### plot phase rectangles #### 
###########################################################################

# prepare cold climate phases as following Rasmussen et al. 2014 + 
# 10.3ka (Bond et al 2001)

event <- data.frame(
  name = c("GS-2", "GI\n1d", "GI\n1c2","GI\n1b","GS-1", 
           "11.4ka\nevent", "10.3ka\nevent", 
            "9.3ka\nevent", "8.2ka\nevent"),
  start = c(lowerlim, 14025, 13610, 13261, 12846, 
            11470, 10375,  9300, 8250),
  end = c(14642, 13904, 13550, 13049, 11653, 
          11350, 10225,  9190, 8090))
event$centre <- event$start-((event$start-event$end)/2)

###########################################################################
par(xpd=NA)
#placement options ####
if(exists("phaserect")==FALSE){message(
"define phaserect to choose placement\nphase rectangle placement options:
1: throughout plot area (selected)
2: part of plot area
3: below plot area
4: above plot area
5: top of plot area
6: bottom of plot area\n   ")
phaserect <- 1
  }

#phase rectangles throughout graph
if(phaserect==1){ message(
"define phaserect to choose placement\nphase rectangle placement options:
1: throughout plot area (selected)
2: part of plot area
3: below plot area
4: above plot area
5: top of plot area
6: bottom of plot area\n   ")
rect(x= event$start,
     xright = event$end,
     ybottom = 0,
     ytop = par("usr")[4],
     col= rgb(0,0,0,alpha=0.2),
     border=NA)
}
#phase rectangles as part of plot area
if(phaserect==2){message(
"define phaserect to choose placement\nphase rectangle placement options:
1: throughout plot area 
2: part of plot area (selected)
3: below plot area
4: above plot area
5: top of plot area
6: bottom of plot area\n   ")
rect(x= event$start,
     xright = event$end,
          ybottom = par("usr")[4]-(0.15*(par("usr")[4])),
          ytop = par("usr")[4]-(0.85*(par("usr")[4])),
     col= rgb(0,0,0,alpha=0.2),
     border=NA)
}

#phase rectangles below plot area
if(phaserect==3){message(
"define phaserect to choose placement\nphase rectangle placement options:
1: throughout plot area 
2: part of plot area 
3: below plot area (selected)
4: above plot area
5: top of plot area
6: bottom of plot area\n   ")
rect(x= event$start,
     xright = event$end,
     ybottom =  par("usr")[3]-(0.07*(par("usr")[4])),
     ytop = par("usr")[3]-(0.02*(par("usr")[4])),
     col= rgb(0,0,1,alpha=0.5),
     border=NA)
}

#phase rectangles above top of graph
if(phaserect==4){message(
"define phaserect to choose placement\nphase rectangle placement options:
1: throughout plot area 
2: part of plot area 
3: below plot area 
4: above plot area (selected)
5: top of plot area
6: bottom of plot area\n   ")
rect(x= event$start,
     xright = event$end,
     ybottom = par("usr")[4]+(0.05*(par("usr")[4])),
     ytop = par("usr")[4]+(0.10*(par("usr")[4])),
     col= rgb(0,0,1,alpha=0.5),
     border=NA)
}

if(exists("l")==FALSE){
  l <- 0.05}

#phase rectangles at top of graph
if(phaserect==5){message(
"define phaserect to choose placement\nphase rectangle placement options:
1: throughout plot area 
2: part of plot area 
3: below plot area 
4: above plot area
5: top of plot area (selected)
6: bottom of plot area\n   ")
rect(x= event$start,
     xright = event$end,
     ybottom = par("usr")[4]-(l*(par("usr")[4])),
     ytop = par("usr")[4],
     col= rgb(0,0,0,alpha=0.3),
     border=NA)
}

#phase rectangles at bottom of graph
if(phaserect==6){message(
"define phaserect to choose placement\nphase rectangle placement options:
1: throughout plot area 
2: part of plot area 
3: below plot area 
4: above plot area
5: top of plot area 
6: bottom of plot area (selected)\n   ")
rect(x= event$start,
     xright = event$end,
     ybottom = 0,
     ytop = par("usr")[3]+(l*(par("usr")[4])),
     col= rgb(0,0,0,alpha=0.3),
     border=NA)
}

#active placement option ####

if(exists("phaselabels")==FALSE ){
message("note: phase labels are not turned on
select value between 0 and 1 for position on y scale")
}else{
text(x=event$centre, 
     y=par("usr")[4]-(phaselabels # plot text at % of plot area
                      *(par("usr")[4])), 
     labels=event$name, col="black")
}

par(xpd=FALSE)
