#Process data.frame
df.trans <- function(dat){
  require(tidyverse)
  dat <- dat[,c("Well","Treatment_Drug","Concentration_1","Unit_1","Concentration_2",
                "Unit_2","Remark","Count","Patient","Plate")]
  dat$Concentration_1[dat$Remark == "Neg_Control"] <- 0
  dat$Concentration_2[dat$Remark == "Neg_Control"] <- 0
  dat <- dat %>%
    mutate(Max = mean(.[Concentration_2 == 0, "Count"])) %>%
    mutate(Percentage = Count*100/Max)
  return(dat)
}

drm.trans <- function(dat){
  drm.dat <- split(dat, dat$Treatment_Drug)
  dmso <- drm.dat$DMSO
  drm.dat <- lapply(drm.dat,function(x) rbind(x,dmso))
  drm.dat <- drm.dat[-which(names(drm.dat) == "DMSO")]
  return(drm.dat)
}

#plots creation in list
plot.list <- function(drc, dat){
  res <- lapply(names(drc), function(x) drpplot(drc[[x]],dat %>% filter(Treatment_Drug == x)))
  names(res) <- names(drc)
  return(res)
}

#turn two layer (deck) list of df into one big dataframe
tdlc <- function(list){
  for(i in seq_along(list)){
    for(j in seq_along(list[[i]])){
      if (j == 1){
        res1 <- list[[i]][[j]]
      } else {
        res1 <- rbind(res1,list[[i]][[j]])
      }
    }
    if (i == 1){
      res2 <- res1
    } else{
      res2 <- rbind(res2,res1)
    }
  }
  return(res2)
}

#turn one layer(deck) list of df into one dataframe
odlc <- function(list){
  for(i in seq_along(list)){
    if (i == 1){
      res1 <- list[[i]]
    } else {
      res1 <- rbind(res1,list[[i]])
    }
  }
  return(res1)
}

#merge list
merge.list <- function(list,by,all = T){
  for(i in seq_along(list)){
    if (i == 1){
      res1 <- list[[i]]
    } else {
      res1 <- merge(res1,list[[i]], by = by, all = all)
    }
  }
  return(res1)
}

# Using predict number - for printing
drpplot <- function(drc.model,dat){
  require(drc)
  pt <- backfit(drc.model)[,"dose"]
  
  fits <- expand.grid(conc=exp(seq(log(pt[2]/10), log(pt[length(pt)]), length=1000)))
  
  pm <- predict(drc.model, newdata=fits, interval="confidence")
  fits$p <- pm[,1]
  fits$pmin <- pm[,2]
  fits$pmax <- pm[,3]
  
  return(fits)
}

plot.drug <- function(dat.mod,fits.mod){
  require(ggplot2)
  g1 <- ggplot(dat.mod,aes(x = Concentration_2,y = Percentage)) +
    geom_point(aes(colour = factor(Patient))) +
    geom_line(data = fits.mod, aes(x=conc, y=p, colour = factor(Patient))) +
    scale_y_continuous(breaks = seq(0, 125, by = 25), limits=c(-10,125)) +
    scale_x_log10(limits = c(min(fits.mod$conc),max(fits.mod$conc)),
                  breaks = 10^(seq(log10(min(fits.mod$conc)),log10(max(fits.mod$conc)),length = log10(max(fits.mod$conc)/min(fits.mod$conc))+1)))+
    geom_hline(yintercept=0)+
    xlab(paste0("Concentration (",unique(dat.mod$Unit_2),")")) +
    labs(colour = "Patient")+
    theme_bw()
  
  return(g1)
}

# Saving ggplot
plotsave <- function(plot, filename, device = "png", scale = 1,
                     width = 1800,
                     height = 1000,
                     units = "px",
                     dpi = 300){
  filename.dev <- paste0(filename,".",device)
  
  pngPath <- normalizePath(file.path(tempdir(), filename.dev))
  
  ggsave(pngPath, plot = plot, device = device, scale = scale,
         width = width,
         height = height,
         units = units,
         dpi = dpi)
  
  return(pngPath)
}

