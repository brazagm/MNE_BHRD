dispeRsal2 <- function(data.predict, model, CI = F, random = T, tax = "family",
         write.result=F){
  require(nlme)
  require(AICcmodavg)
  thisVersion <- 0.2
#  currentVersion <- scan("http://www.botany.ut.ee/dispersal/version.txt")
  currentVersion <- 0.2 # modified by Alan Braz
  TPL_nomatch <- NULL
  
  if (is.numeric(model)){
    traits <- switch(model, c("DS", "GF","TV"), c("DS", "GF", "SM","RH"), c("DS", "GF", "RH"),
                     c("DS", "GF", "SM"), c("DS","GF"))
  } else{
    traits <- model
  }
  
  if(any(is.na(match(traits, names(data.predict))))) {
    stop("variable names do not match required names")
  }
  
  if(any(is.na(match(levels(data.predict$GF), levels(model.data$GF))))) {
    stop("wrong category name in column GF")
  }
  
  if(any(is.na(match(levels(data.predict$DS), levels(model.data$DS))))) {
    stop("wrong category name in column DS")
  }
  
  model.data.traits <- model.data[,2:6]
  model.data.rest <- model.data[,c(1,7:9)]
  ind <- names(model.data) %in% traits
  md <- na.omit(data.frame(model.data.rest, model.data[ind]))
  pd <- data.predict[,na.omit(match(names(md),
                                    names(data.predict)))]
  if(random == TRUE){
    pd.species.tpl <- TPLMod(as.character(pd$Species))
    pd$FamilyTPL <- pd.species.tpl$Family
    TPL_nomatch <- pd[which(model.data$FamilyTPL=="" |
                              is.na(model.data$FamilyTPL)),]$Species
    if (any(pd$FamilyTPL=="" | is.na(pd$FamilyTPL)) ==F) {
      pd <- pd
    } else{
      pd <- pd[-which(pd$FamilyTPL=="" | is.na(pd$FamilyTPL)),]
    }
    ln <- length(pd$FamilyTPL)
    Order <- vector("character", ln)
    for (i in 1:ln) {
      Order[i] <- as.character(OrderFamilies$Order[match(pd$FamilyTPL[i],
                                                         OrderFamilies$Family)])
    }
    pd$OrderTPL <- Order
  }
  
  pd <- na.omit(pd)
  
  if(dim(pd)[1] == 0) {
    stop(paste("not enough data to run specified model with traits", paste(traits, collapse = " + ")))
  }
  
  formu <- as.formula(paste("logMax ~", paste(traits, collapse="+")))
  m <- lm(formu,  data=md)
  
  if(CI == T) {
    log10MDD <- predict(m, na.omit(pd), interval = "confidence")
    colnames(log10MDD) <- c("log10MDD", "log10MDD_lwrCL", "log10MDD_uppCL")
  } else{
    log10MDD <- predict(m, na.omit(pd))
  }
  
  Species <- as.character(na.omit(pd)$Species)
  assign("formu", formu,  envir = .GlobalEnv)
  
  if(random == T) {
    
    if (tax == "family") {
      rs <- reStruct(object = ~ 1 | OrderTPL/FamilyTPL, pdClass="pdDiag")
      level.max <- 2
      m <- lme(formu, data = md, random = rs)
      m$call$fixed <- as.call(formu)
      
      if (dim(na.omit(pd))[1] == 1) {
        log10MDD_Family <- data.frame(predict(m, na.omit(pd),
                                              level = 0:level.max))[3, 3]
        log10MDD_Order <- data.frame(predict(m, na.omit(pd),
                                             level = 0:level.max))[2, 3]
      } else{
        log10MDD_Family <- data.frame(predict(m, na.omit(pd),
                                              level = 0:level.max))[,-c(1,2)][,3]
        log10MDD_Order <- data.frame(predict(m, na.omit(pd),
                                             level = 0:level.max))[,-c(1,2)][,2]
      }
      
      if (CI == TRUE){
        se <- predictSE.lme(m, na.omit(pd))$se.fit
        log10MDD_Family <- as.data.frame(log10MDD_Family)
        log10MDD_Family$log10MDD_Family_lwrCL <-
          log10MDD_Family$log10MDD_Family - 2 * se
        log10MDD_Family$log10MDD_Family_uppCL <-
          log10MDD_Family$log10MDD_Family + 2 * se
        log10MDD_Order <- as.data.frame(log10MDD_Order)
        log10MDD_Order$log10MDD_Order_lwrCL <-
          log10MDD_Order$log10MDD_Order - 2 * se
        log10MDD_Order$log10MDD_Order_uppCL <-
          log10MDD_Order$log10MDD_Order + 2 * se
      }
      
      if (any(is.na(log10MDD_Family))) {
        p = data.frame(Species, log10MDD_Family, log10MDD_Order)
      } else{
        p = data.frame(Species, log10MDD_Family)
      }
      
      if (any(is.na(log10MDD_Order))) {
        p = data.frame(Species, log10MDD_Family, log10MDD_Order, log10MDD)
      }
    }
    
    if (tax == "order") {
      rs <- reStruct(object = ~ 1 | OrderTPL, pdClass="pdDiag")
      level.max <- 1
      m <- lme(formu, data = md, random = rs)
      m$call$fixed <- as.call(formu)
      
      if (dim(na.omit(pd))[1] == 1) {
        log10MDD_Order <- data.frame(predict(m, na.omit(pd),
                                             level = 0:level.max))[2,2]
      } else{
        log10MDD_Order <- data.frame(predict(m, na.omit(pd),
                                             level = 0:level.max))[,-c(1,2)]
      }
      
      if (CI == TRUE){
        se <- predictSE.lme(m, na.omit(pd))$se.fit
        log10MDD_Order <- as.data.frame(log10MDD_Order)
        log10MDD_Order$log10MDD_Order_lwrCL <-
          log10MDD_Order$log10MDD_Order - 2 * se
        log10MDD_Order$log10MDD_Order_uppCL <-
          log10MDD_Order$log10MDD_Order + 2 * se
      }
      
      if (any(is.na(log10MDD_Order))) {
        p = data.frame(Species, log10MDD_Order, log10MDD)
      } else{
        p = data.frame(Species, log10MDD_Order)
      }
    }
    
  } else{
    p <-  data.frame(Species, log10MDD)
  }
  
  if(any(traits == "DS")) {
    DS <- pd[na.omit(pmatch(p$Species, pd$Species)),]$DS
  } else {
    DS <- NULL
  }
  SpecDS.p <- paste(p$Species, DS)
  SpecDS.md <- paste(md$Species, md$DS)
  md.s <- md[which(match(SpecDS.md, SpecDS.p)!="NA"),]
  p.s <- p[which(match(SpecDS.p, SpecDS.md)!="NA"),]
  r <-  rownames(p[which(match(SpecDS.p, SpecDS.md)!="NA"),])
  if(length(r) > 0){
    log10MDD_measured <- as.vector(rep(NA, dim(p)[1]), mode="numeric")
    p <- cbind(p, log10MDD_measured)
    p[as.character(r),]$log10MDD_measured <- md.s[match(p.s$Species, md.s$Species),]$logMax
  } 
  
  if(random == T) {
    Order <- pd[na.omit(match(p$Species, pd$Species)),]$OrderTPL
    Family <- pd[na.omit(match(p$Species, pd$Species)),]$FamilyTPL
    if(is.null(DS)) {
      p <- cbind(p, Order, Family)
      cs <- dim(p)[2]
      p <- p[,c(1, cs-1, cs, c(2:(cs-2)))]
      p <- p[order(p$Family, p$Species),]    
    } else {
      p <- cbind(p, DS, Order, Family) 
      cs <- dim(p)[2]
      p <- p[,c(1, cs-1, cs, cs-2, c(2:(cs-3)))]
      p <- p[order(p$Family, p$Species),]
    }
    
  } else{
    p <- cbind(p, DS)
    cs <- dim(p)[2]
    p <- p[,c(1, cs, c(2:(cs-1)))]
  }
  
  remove(formu,  envir = .GlobalEnv)
  
  if (thisVersion != currentVersion) {
    ver = cat( "\n\n\nATTENTION\n",
               "You are using dispeRsal version ", thisVersion , ".\n",
               "You can download the more up to date\n",
               "version ", currentVersion , " at www.botany.ut/dispersal\n\n\n", sep="")
  } else{
    ver = cat("You are using dispeRsal version ", thisVersion, ".\n\n", sep="")
  }
  out <- list(p, TPL_nomatch)
  names(out) <- c("predictions", "unmatched_species")
  if(write.result==T) {
    write.table(out[[1]] , "predictedDD.txt")
    write.table(out[[2]] , "unmatched.txt")
  }
  out
}
