#############################################################################################################################################
###  plot.bfast function   ##################################################################################################################
# bfast uses the embeded function plot.bfast to plot results from timeseries decomposition. 
# To bring up code for plot.bfast: bfast:::plot.bfast
# I have modidied sections of this to allow plotting of NPP and precip data
#############################################################################################################################################

plot.bfast <- function (x, type = c("components", "all", "data", "seasonal", "trend", "noise"), sim = NULL, largest = FALSE, main, ANOVA = FALSE, npp, precip, ...) 
{
  type <- match.arg(type)
  realdata <- is.null(sim)
  Trend.bp <- !x$nobp$Vt
  if (type == "largest" & !Trend.bp) 
    stop("No trend breakpoints")
  title <- !missing(main)
  niter <- length(x$output)
  out <- x$output[[niter]]
  Tt <- out$Tt
  St <- out$St
  noise <- out$Nt
  npp <- npp
  precip <- precip
  if (type == "data") {
    if (!title) 
      main <- "Yt"
    plot(x$Yt, main = main, ...)
  }
  else if (type == "components") {
    ft <- cbind(seasonal = out$St, trend = out$Tt, remainder = out$Nt)
    tsp(ft) <- tsp(x$Yt)
    ft <- list(time.series = ft)
    if (!title) 
      main <- paste("no. iterations to estimate breakpoints:", 
        niter)
    if (ANOVA == TRUE) {
      seasonal(ft, out, npp, precip, sim = sim, main = main, fit = x)
    }
    else {
      seasonal(ft, out, npp, precip,  sim = sim, main = main)
    }
  }
  else if (type == "noise") {
    require(forecast)
    if (!title) 
      main <- "Noise component"
    tsdisplay(noise, main = main)
  }
  else {
    if (type == "all") {
      idx <- 1:niter
      opar <- par(mfrow = c(2, niter))
    }
    else idx <- niter
    for (i in idx) {
      out <- x$output[[i]]
      if (type != "seasonal") {
        if (type == "trend" & !title) 
          main <- "Trend component"
        else if (!title) 
          main <- paste("Iteration ", i, ": Trend", sep = "")
        plot(out$Vt, main = main, ylab = "Vt", ...)
        lines(out$Tt, col = 4)
        if (Trend.bp) {
          lines(out$bp.Vt)
          lines(out$ci.Vt)
          legend("topright", paste("Time of BP(s)", paste(out$Vt.bp, 
            collapse = ",")), col = 2)
        }
        if (!realdata) {
          lines(sim$time.series[, "abrupt"], col = 1, 
            lty = 2)
          legend("bottomleft", c("estimated", "simulated"), 
            lty = c(1, 2), col = 1)
        }
        if (largest) {
          legend("bottomright", c("Magnitude of most sign change"), 
            lty = c(1), col = 6)
          lines(x$jump, col = 6)
          points(x$jump, pch = 14, cex = 1, col = 6)
        }
      }
      if (type != "trend") {
        if (type == "seasonal" & !title) 
          main <- "Seasonal component"
        else if (!title) 
          main <- paste("Iteration ", i, ": Seasonal", 
            sep = "")
        plot(out$Wt, main = main, ylab = "Wt", ...)
        lines(out$St, col = 2)
        Seas.bp <- !x$nobp$Wt
        if (Seas.bp) {
          lines(out$bp.Wt)
          lines(out$ci.Wt)
          legend("topright", paste("Time of BP(s)", paste(out$Wt.bp, 
            collapse = ",")), col = 2)
        }
        if (!realdata) {
          lines(sim$time.series[, "seasonal"], col = 1, 
            lty = 2)
          legend("bottomleft", c("first run seasonality", 
            "first run estimated", "simulated"), lty = c(1, 
              1, 2), col = c(1, 2, 1))
        }
      }
    }
    if (type == "all") 
      par(opar)
  }
}





##################################################################################
#Plotting all components of the bfast output uses the seasonal function. I have modified this function to allow plotting of NPP and precip.
##Seasonal function


seasonal <- function (x, out, npp, precip, sim = NULL, labels = c(colnames(X), "Precipitation", "NPP"), set.pars = list(tck = -0.01, 
  mar = c(0, 6, 0, 6), oma = c(6, 0, 4, 0)), main = NULL, range.bars = FALSE, 
  col.range = "light gray", fit = NULL) 
{
  layout(matrix(c(1, 2, 3, 4,5,6), 6, 1, byrow = TRUE), heights = c(1, 
    0.75, 1, 0.75, 1, 1), TRUE) #modify trend plot height, change from 1.5 to 1
  sers <- x$time.series
  npp<-npp
  precip<-precip
  ncomp <- ncol(sers)
  data <- drop(sers %*% rep(1, ncomp))
  X <- cbind(data, sers)
  X<-ts(as.zoo(X)[1:291], start=2000.130, frequency=22.82812)
  colnames(X) <- c("Data", "Seasonal", "Trend", "Resid") #change plot x axis names from c("Yt", "St", "Tt", "et")
  nplot <- ncomp + 3
  if (range.bars) 
    mx <- min(apply(rx <- apply(X, 2, range), 2, diff))
  if (length(set.pars)) {
    oldpar <- do.call("par", as.list(names(set.pars)))
    on.exit(par(oldpar))
    do.call("par", set.pars)
  }
  if (is.null(fit) == F) {
    niter <- length(fit$output)
    out <- fit$output[[niter]]
    out_ANOVA <- array()
    out_breakdates <- array()
    if (out$Vt.bp[1] > 0) {
      breaks <- length(out$Vt.bp)
    }
    else {
      breaks <- 0
    }
    if (breaks > 0) {
      breakdates <- out$Vt.bp
      coefs <- coef(out$bp.Vt)
      sl <- coefs[, 2]
    }
    TS_anova <- fit$Yt - out$St
    dataframe <- data.frame(TIME = c(1:length(fit$Yt)), DATA = TS_anova)
    for (m in 1:(breaks + 1)) {
      startpoint <- if (m == 1) 
        1
      else breakdates[[m - 1]]
      endpoint <- if (m == (breaks + 1)) 
        length(fit$Yt)
      else breakdates[m] - 1
      df2 <- dataframe[startpoint:endpoint, ]
      model <- lm(DATA ~ TIME, data = df2)
      modelAnova <- anova(model)
      out_ANOVA[m] <- modelAnova$Pr[1]
      if (breaks == 0) {
        sl <- model$coefficients[2]
      }
    }
  }
  
  for (i in 1:nplot) {
    if (i == 6) {
      par(mar = c(0, 6, 0, 6))
    }
    if (i<5) {
      plot(X[, i], col = if (i == 1) 
        "black"
        else "red", ylim = if (i == 1 | i == 3) 
          range(X[, 1])
        else range(X[, i], sim$time.series[, i - 1], na.rm = TRUE), 
        type = if (i != 4) 
          "l"
        else "h", xlab = "", ylab = "", axes = FALSE) }
    
    if (i==5) {
      plot(precip, ann=FALSE, xaxt="n", yaxt="n", )
      
    }
    if (i==6) {
      plot(npp, type="b", ann=FALSE, yaxt="n", ylim=c(0,round(max(npp)+15))) 
      
    }    
    
    if (range.bars) {
      dx <- 1/64 * diff(ux <- par("usr")[1:2])
      y <- mean(rx[, i])
      rect(ux[2] - dx, y + mx/2, ux[2] - 0.4 * dx, y - 
          mx/2, col = col.range, xpd = TRUE)
    }
    if (i == 1 && !is.null(main)) {
      mtext(main, side = 3, font = 1.5, line = 1, cex = 1.1)
      if (!is.null(sim)) {
        lines(X[, i + 1] + X[, i + 2], col = "red", type = "l")
        legend("topleft", c("input", "estimated seasonal + trend "), 
          col = c("black", "red"), lty = 1, cex=0.7)
      }
    }
    if (i == 2) {
      lines(sim$time.series[, "seasonal"], col = "black")
      lines(out$bp.Wt)
      lines(out$ci.Wt)
    }
    if (i == 3) {
      lines(sim$time.series[, "abrupt"], col = "black")
      lines(out$bp.Vt)
      lines(out$ci.Vt)
      if (is.null(fit) == FALSE) {
        for (m in 1:(breaks + 1)) {
          x_coor <- out$bp.Wt$datatsp[[1]]
          if (m > 1) {
            x_coor <- x_coor + breakdates[[m - 1]]/frequency(fit$Yt)
          }
          y_range <- range(X[, 1])
          y_sl <- y_range[2] - (y_range[2] - y_range[1])/10
          y_Pr <- y_range[2] - (y_range[2] - y_range[1])/5
          beta <- formatC(sl[m], format = "f", digits = 3)
          text(x_coor, y_sl, bquote(beta == .(beta)), 
            pos = 4, cex=0.8)
          Pr <- formatC(out_ANOVA[m], format = "f", digits = 3)
          text(x_coor, y_Pr, bquote(p == .(Pr)), pos = 4, cex=0.8)
        }
      }
    }
    if (i == 4) {
      abline(h = 0)
      lines(sim$time.series[, "remainder"], col = "black")
    }
    if (i == 5) {
      lines(out$bp.Wt, col="blue", lty=2)    #this doesn't work
      lines(out$bp.Vt, col="green", lty=2)   #this doesn't work
    }
    if (i == 6) {
      lines(out$bp.Wt, col="blue", lty=2)   #this doesn't work
      lines(out$bp.Vt, col="green", lty=2)  #this doesn't work
    }
    
    
    box()
    right <- i%%2 == 0
    par(las=1)
    axis(2, labels = !right)
    axis(4, labels = right)
    #axis(1, labels = i == nplot)
    par(las=0)
    mtext(labels[i], side = 2, 3, cex=0.7)
  }
  mtext("Time", side = 1, line = 3, cex=0.8)
  invisible()
  if (is.null(fit) == FALSE) {
    return(data.frame(slope = sl, prob = out_ANOVA))
  }
  layout(matrix(1))
}

############################################################################################################################################
### Plotting NPP data by plant functional tyep ###############
############################################################################################################################################

jrn_func <- read.csv("G:/ARS_Static_Data/Jornada/JORNADA_NPP/biomass_funcgrp_site_2000_2012.csv", header=T, sep=',')

BASN.jrn_func <- jrn_func[jrn_func$SITE %in% c("BASN"), ]
CALI.jrn_func  <- jrn_func[jrn_func$SITE %in% c("CALI"), ]
COLL.jrn_func  <- jrn_func[jrn_func$SITE %in% c("COLL"), ]
EAST.jrn_func  <- jrn_func[jrn_func$SITE %in% c("EAST"), ]
GRAV.jrn_func  <- jrn_func[jrn_func$SITE %in% c("GRAV"), ]
IBPE.jrn_func  <- jrn_func[jrn_func$SITE %in% c("IBPE"), ]
NORT.jrn_func  <- jrn_func[jrn_func$SITE %in% c("NORT"), ]
RABB.jrn_func  <- jrn_func[jrn_func$SITE %in% c("RABB"), ]
SAND.jrn_func  <- jrn_func[jrn_func$SITE %in% c("SAND"), ]
SMAL.jrn_func  <- jrn_func[jrn_func$SITE %in% c("SMAL"), ]
SUMM.jrn_func  <- jrn_func[jrn_func$SITE %in% c("SUMM"), ]
TAYL.jrn_func  <- jrn_func[jrn_func$SITE %in% c("TAYL"), ]
TOBO.jrn_func  <- jrn_func[jrn_func$SITE %in% c("TOBO"), ]
WELL.jrn_func  <- jrn_func[jrn_func$SITE %in% c("WELL"), ]
WEST.jrn_func  <- jrn_func[jrn_func$SITE %in% c("WEST"), ]


######################################################################################################################################
#######################################################################################################################################
#####################################################################
############  Example plotting for BASN   ###############################################
#This code was used to create the NPP functional plots for the other 14 sites

#BASN Absolute Contribution
BASN.jrn_func$median_date <- as.Date(as.character(BASN.jrn_func$median_date))

ylim <- c(0, max(BASN.jrn_func$bio_annual + BASN.jrn_func$bio_pforb +  BASN.jrn_func$bio_pgrass +  BASN.jrn_func$bio_shrub))
xx <- c(BASN.jrn_func$median_date, rev(BASN.jrn_func$median_date))

basn.ab.ann <-BASN.jrn_func$bio_annual
basn.ab.pforb <- BASN.jrn_func$bio_pforb
basn.ab.grass <- BASN.jrn_func$bio_pgrass
basn.ab.shrub <- BASN.jrn_func$bio_shrub


######################################################order grass, shrub,  ann, pforb
CairoPDF("BASN_AbsCont_FunGrps.pdf",
  7, 5, bg = "transparent", pointsize=14 )
par(mar=c(6.5,4,3.5,4)+.1, las=1)
#order herb,  shrub, grass, ann, pforb
yy_grass <- c(rep(0, nrow(BASN.jrn_func)), rev(basn.ab.grass))
plot(x=BASN.jrn_func$median_date, y=basn.ab.grass, ylim=ylim, col='red', type='l', xaxt='n',
  ylab='', xlab='Date', main='BASN')
mtext(expression(paste('Mean Biomass (g ',m^2,')')), side=2, line=2.5, las=0)

polygon(xx, yy_grass, col="#00A600FF")

yy_shrub <- c( basn.ab.grass , rev(basn.ab.shrub) + rev(basn.ab.grass))
polygon(xx, yy_shrub, col="#E6E600FF")

yy_ann <- c(basn.ab.grass +basn.ab.shrub,  rev(basn.ab.grass) + rev(basn.ab.shrub) + rev(basn.ab.ann))
polygon(xx, yy_ann, col="#EBB25EFF")

yy_pforb <-  c( basn.ab.grass +basn.ab.shrub + basn.ab.ann,  rev(basn.ab.grass) + rev(basn.ab.shrub) + rev(basn.ab.ann) + rev(basn.ab.pforb))
polygon(xx, yy_pforb, col="#F0C9C0FF")

# x axis labels
labdates <- as.Date('1970-01-01') + min(BASN.jrn_func$median_date):max(BASN.jrn_func$median_date)
labdates <- labdates[ format(labdates, '%m') == '01']
labdates <- labdates[ format(labdates, '%d') == '01']
labdates <- unique(c(labdates, max(BASN.jrn_func$median_date)))

labnames <- format(labdates, '%Y')
axis(1, at=labdates, labels=labnames)

# # black lines first day of month
# for(a in labdates[ format(labdates, '%y') == '01']){
#   abline(v=a)
# }

# make the first 3 lower to avoid our legend "#00A600FF" "#63C600FF"  "#E6E600FF" "#EBB25EFF" "#F2F2F2FF"
cnt <- 0

# set alpha = 80 so it's relatively transparent
color <- rgb(190, 190, 190, alpha=80, maxColorValue=255)

# check if the first data point is a Sunday
BASN.jrn_func$mon <- format( BASN.jrn_func$median_date, '%m' )

# plot strip along sept biomass max
for( i in BASN.jrn_func[ BASN.jrn_func$mon == '09', ]$median_date ){
  rect(xleft=i-30, xright=i+30, ybottom=-1000, ytop=1.1*ylim[2], density=100, col=color)
  cnt <- cnt + 1
}
par(xpd=TRUE)
#order herb,  shrub, grass, ann, pforb
legend(x=BASN.jrn_func[1,]$median_date -100 , y=ylim[2]+5, c('Grass', 'Shrubs','Annuals', 'Forbs'), fill=c("#00A600FF",  "#E6E600FF", "#EBB25EFF" ,"#F0C9C0FF"), cex=0.67, bg = "white")
par(xpd=F)
dev.off()

######################################################################################################################################
#######################################################################################################################################

