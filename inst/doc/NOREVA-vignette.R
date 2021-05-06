## ----setup, message=FALSE-----------------------------------------------------
library(NOREVA)

## ----figurea, echo = FALSE, fig.align='center', fig.cap ='NOREVA-All-Criteria-Output-Figures', out.width ="650px"----
knitr::include_graphics ("./sampledata/NOREVA-Output-All-Criteria-Figures.png")

## ----kable2-------------------------------------------------------------------
allrankings <- read.csv(file = "./sampledata/OUTPUT-NOREVA-Overall.Ranking.Data.csv",header = T)
head(allrankings)

## ----pressure, echo = FALSE, fig.align='center', fig.cap ='NOREVA-Ranking-Top.100.workflows', out.width ="650px"----
knitr::include_graphics ("./sampledata/NOREVA-Ranking-Top.100.workflows.png")

