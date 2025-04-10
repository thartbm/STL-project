
# generic plot functions -----

setupFigureFile <- function(target='inline',width=8,height=6,dpi=300,filename) {
  
  #filename=sprintf('doc/fig3_taskerror_parmin.%s',target)
  
  if (target == 'pdf') {
    pdf(file   = filename, 
        width  = width, 
        height = height)
  }
  if (target == 'svg') {
    svglite::svglite( filename = filename,
                      width = width,
                      height = height,
                      fix_text_size = FALSE) 
                      # fix_text_size messes up figures on my machine... 
                      # maybe it's better on yours?
  }
  if (target == 'png') {
    png( filename = filename,
         width = width*dpi,
         height = height*dpi,
         res = dpi
    )
  }
  if (target == 'tiff') {
    tiff( filename = filename,
          compression = 'lzw',
          width = width*dpi,
          height = height*dpi,
          res = dpi
    )
  }
}


# task error plots -----

plotDataParmin <- function(models=FALSE, target='inline') {

  if (models) {
    filename=sprintf('doc/fig4_model_fits_parmin.%s',target)
  } else {
    filename=sprintf('doc/fig2_alldata_parmin.%s',target)
  }
  width = 8
  if (target %in% c('pdf','svg','png','tiff')) {
    setupFigureFile( target   = target,
                     width    = width,
                     height   = 4,
                     dpi      = 300,
                     filename = filename)
  }
  
  
  left  <- 0.75
  right <- 0.05
  
  deg15 <- (width - ((left + right) * 2)) / 7
  
  layout(mat=matrix(c(1:2), ncol=2), widths=c(deg15*3,deg15*4))
  par(mai=c(0.75,left,0.4,right))
  
  
  for (maxrot in c(45,60)) {
    
    plot(x=-1000,y=-1000,
         main=sprintf('%d° max rotation',maxrot),ylab='',xlab='',
         xlim=c(0,maxrot+1),
         # xlim=c(0,91),
         ylim=c(-0.5,8.5),
         ax=F,bty='n')
    
    title(ylab='initial change [°]',line=2.25)
    title(xlab='rotation [°]',line=2.25)
    
    df <- loadSTLdata(average=median, maxrots=c(maxrot))
    
    if (models) {
      fits <- read.csv(sprintf('data/modelFits_%d.csv', maxrot), stringsAsFactors = FALSE)
      fits <- fits[which(fits$participant=='all'),]
    }
    
    CIdf <- aggregate(response ~ target + rotation, data=df, FUN=Reach::getConfidenceInterval, method='b', resamples=2000)
    avgdf <- aggregate(response ~ target + rotation, data=df, FUN=mean)
    for (tt in c('point','arc')) {
      
      if (tt == 'point') {
        colors <- c('#0066FFFF', '#0066FF33')
      }
      if (tt == 'arc') {
        colors <- c('#FF6600FF', '#FF660033')
      }
      
      tCI <- CIdf[which(CIdf$target == tt),]
      tavg <- avgdf[which(avgdf$target == tt),]
      
      X <- c(tCI$rotation, rev(tCI$rotation))
      Y <- c(tCI$response[,1], rev(tCI$response[,2]))
      polygon(X,Y,border=NA,col=colors[2])
      
      if (models) {
        fit <- fits[which(fits$target == tt),]
        # cross <- fit$c/fit$r
        
        lines(x=c(0,fit$c/fit$r,maxrot),
              y=c(0,fit$c,fit$c),
              lty=3,
              col=colors[1],
              lw=2)
        
        attpar = c('s'=fit$s, 'w'=fit$w)
        rots = seq(0,maxrot,0.5)
        attribution_predictions <- STLPredict(par=attpar,
                                              rotations=rots)
        
        lines(x=rots,
              y=attribution_predictions,
              lty=2,
              col=colors[1],
              lw=2)
      } else {
        lines(x=tavg$rotation, y=tavg$response, col=colors[1], lw=2)
      }
      
    }
    
    if (maxrot == 45) {
      axis(side=1,at=c(1,5,10,15,20,25,30,35,40,45),cex.axis=0.75)
    }
    if (maxrot == 60) {
      axis(side=1,at=c(1,5,10,15,20,25,30,40,50,60),cex.axis=0.75)
    }
    if (maxrot == 90) {
      axis(side=1,at=c(1,5,10,15,20,30,40,50,70,90),cex.axis=0.75)
    }
    axis(side=2,at=c(0,2,4,6,8),cex.axis=0.75)
    
    
    if (maxrot == 60) {
      legend(x=35,y=2,
             legend=c('dot target', 'arc target'),
             lwd=c(2,2),
             col=c('#0066FFFF','#FF6600FF'),
             bty='n')
    }
    if (maxrot == 60) {
      if (models) {
        legend(x=20,y=2,
               legend=c('attribution model', 'capped model'),
               lwd=c(2,2),
               lty=c(2,3),
               col=c('#999999','#999999'),
               bty='n')
      }
    }
    
  }
  
  if (target %in% c('pdf','svg','png','tiff')) {
    dev.off()
  }
  
}


plotTaskErrorEffects <- function(target='inline',posthocs=TRUE) {

  
  if (target %in% c('pdf','svg','png','tiff')) {
    setupFigureFile( target   = target,
                     width    = 4,
                     height   = 8,
                     dpi      = 300,
                     filename = sprintf('doc/fig3_taskerror.%s',target))
  }
  
  
  layout(mat=matrix(c(1:3), ncol=1))
  par(mar=c(3.25,3.25,2,0.1))
  
  pval.df <- taskErrorPostHocs()
  
  for (maxrot in c(45,60,90)) {
    
    plot(x=-1000,y=-1000,
         main=sprintf('%d° max rotation',maxrot),ylab='',xlab='',
         # xlim=c(0,maxrot+1),
         xlim=c(0,91),
         ylim=c(-2.5,4.5),
         ax=F,bty='n')
    
    title(ylab='reach aftereffect [°]',line=2.25)
    title(xlab='rotation [°]',line=2.25)
    
    df <- loadSTLdata(average=median, maxrots=c(maxrot))
    
    diffdf <- aggregate(response ~ participant + rotation, data=df, FUN=diff)
    diffdf$response <- diffdf$response
    
    CIdf <- aggregate(response ~ rotation, data=diffdf, FUN=Reach::getConfidenceInterval, method='b', resamples=2000)
    avgdf <- aggregate(response ~ rotation, data=diffdf, FUN=mean)
    
    colors <- c('#9900FFFF', '#9900FF33')
    
    lines(x=c(1,maxrot),y=c(0,0),col='#999999',lty=2)
    
    # tCI <- CIdf[which(CIdf$target == target),]
    # tavg <- avgdf[which(avgdf$target == target),]
    tCI <- CIdf
    tavg <- avgdf
    
    X <- c(tCI$rotation, rev(tCI$rotation))
    Y <- c(tCI$response[,1], rev(tCI$response[,2]))
    polygon(X,Y,border=NA,col=colors[2])
    
    lines(x=tavg$rotation, y=tavg$response, col=colors[1])
    
    # add FDR corrected p-values:
    mr.pval.df <- pval.df[which(pval.df$group == maxrot),]
    
    if (posthocs) {
      for (rotation in mr.pval.df$rotation) {
        
        pval <- mr.pval.df$pval[which(mr.pval.df$rotation == rotation)]
        pval.adj <- mr.pval.df$pval.adj[which(mr.pval.df$rotation == rotation)]
        
        pch <- 1
        col <- '#9900FF33'
        if (pval < 0.05) {
          col <- '#9900FFFF'
        }
        if (pval.adj < 0.05) {
          pch <- 16
        }
        
        points( x = rotation,
                y = -2,
                pch=pch,
                col=col)
        
      }
    }
    
    if (maxrot == 45) {
      axis(side=1,at=c(1,5,10,15,20,25,30,35,40,45),cex.axis=0.75)
    }
    if (maxrot == 60) {
      axis(side=1,at=c(1,5,10,15,20,25,30,40,50,60),cex.axis=0.75)
    }
    if (maxrot == 90) {
      axis(side=1,at=c(1,5,10,15,20,30,40,50,70,90),cex.axis=0.75)
    }
    axis(side=2,at=c(-2,0,2,4),cex.axis=0.75)
    
    
  }
  
  if (target %in% c('pdf','svg','png','tiff')) {
    dev.off()
  }
  
  
}

plotTaskErrorEffectsParmin <- function(target='inline',posthocs=TRUE) {
  
  width = 8
  if (target %in% c('pdf','svg','png','tiff')) {
    setupFigureFile( target   = target,
                     width    = width,
                     height   = 4,
                     dpi      = 300,
                     filename = sprintf('doc/fig3_taskerror_parmin.%s',target))
  }
  
  left  <- 0.75
  right <- 0.05
  
  deg15 <- (width - ((left + right) * 2)) / 7
  
  layout(mat=matrix(c(1:2), ncol=2), widths=c(deg15*3,deg15*4))
  par(mai=c(0.75,left,0.4,right))
  
  
  pval.df <- taskErrorPostHocs()
  
  for (maxrot in c(45,60)) {
    
    plot(x=-1000,y=-1000,
         main=sprintf('%d° max rotation',maxrot),ylab='',xlab='',
         xlim=c(0,maxrot+1),
         # xlim=c(0,91),
         ylim=c(-2.5,4.5),
         ax=F,bty='n')
    
    title(ylab='initial change difference [°]',line=2.25)
    title(xlab='rotation [°]',line=2.25)
    
    df <- loadSTLdata(average=median, maxrots=c(maxrot))
    
    diffdf <- aggregate(response ~ participant + rotation, data=df, FUN=diff)
    diffdf$response <- diffdf$response
    
    CIdf <- aggregate(response ~ rotation, data=diffdf, FUN=Reach::getConfidenceInterval, method='b', resamples=2000)
    avgdf <- aggregate(response ~ rotation, data=diffdf, FUN=mean)
    
    colors <- c('#9900FFFF', '#9900FF33')
    
    lines(x=c(1,maxrot),y=c(0,0),col='#999999',lty=2)
    
    # tCI <- CIdf[which(CIdf$target == target),]
    # tavg <- avgdf[which(avgdf$target == target),]
    tCI <- CIdf
    tavg <- avgdf
    
    X <- c(tCI$rotation, rev(tCI$rotation))
    Y <- c(tCI$response[,1], rev(tCI$response[,2]))
    polygon(X,Y,border=NA,col=colors[2])
    
    lines(x=tavg$rotation, y=tavg$response, col=colors[1])
    
    # add FDR corrected p-values:
    mr.pval.df <- pval.df[which(pval.df$group == maxrot),]
    
    if (posthocs) {
      for (rotation in mr.pval.df$rotation) {
        
        pval <- mr.pval.df$pval[which(mr.pval.df$rotation == rotation)]
        pval.adj <- mr.pval.df$pval.adj[which(mr.pval.df$rotation == rotation)]
        
        pch <- 1
        col <- '#9900FF33'
        if (pval < 0.05) {
          col <- '#9900FFFF'
        }
        if (pval.adj < 0.05) {
          pch <- 16
        }
        
        points( x = rotation,
                y = -2,
                pch=pch,
                col=col)
        
      }
    }
    
    if (maxrot == 45) {
      axis(side=1,at=c(1,5,10,15,20,25,30,35,40,45),cex.axis=0.75)
    }
    if (maxrot == 60) {
      axis(side=1,at=c(1,5,10,15,20,25,30,40,50,60),cex.axis=0.75)
    }
    if (maxrot == 90) {
      axis(side=1,at=c(1,5,10,15,20,30,40,50,70,90),cex.axis=0.75)
    }
    axis(side=2,at=c(-2,0,2,4),cex.axis=0.75)
    
    
  }
  
  if (target %in% c('pdf','svg','png','tiff')) {
    dev.off()
  }
  
  
}


plotAverageRotIRDs <- function(target='inline') {
  
  width = 8
  if (target %in% c('pdf','svg','png','tiff')) {
    setupFigureFile( target   = target,
                     width    = width,
                     height   = 4,
                     dpi      = 300,
                     filename = sprintf('doc/fig3+_average_rot_IRDS.%s',target))
  }
  
  left  <- 0.75
  right <- 0.05
  
  deg15 <- (width - ((left + right) * 2)) / 7
  
  layout(mat=matrix(c(1:2), ncol=2))
  par(mai=c(0.75,left,0.4,right))
  
  # set up color map:
  
  df <- loadSTLdata(average=median, maxrots=c(45,60))
  all_rots <- sort(unique(df$rotation))
  
  colsT <- c()
  colsO <- c()
  for (rot in c(0:60)) {
    colsT <- c(colsT, Reach::colorMix('#FF000033','#0000FF33',c(rot,max(all_rots)-rot)) )
    colsO <- c(colsO, Reach::colorMix('#FF0000','#0000FF',c(rot,max(all_rots)-rot)) )
  }
  
  
  
  for (maxrot in c(45,60)) {
    
    plot(x=-1000,y=-1000,
         main=sprintf('%d° max rotation',maxrot),ylab='',xlab='',
         # xlim=c(0,maxrot+1),
         xlim=c(-12,18),
         ylim=c(-12,18),
         ax=F,bty='n',asp=1)
    
    title(ylab='arc initial reach deviation [°]',line=2.25)
    title(xlab='dot initial reach deviation [°]',line=2.25)
    
    lines(x=c(-12,18),
          y=c(0,0),
          lty=1,col='#999999')
    lines(x=c(0,0),
          y=c(-12,18),
          lty=1,col='#999999')
    lines(x=c(-12,18),
          y=c(-12,18),
          lty=1,col='#999999')
    
    df <- loadSTLdata(average=median, maxrots=c(maxrot))
    
    rotations <- unique(df$rotation)
    
    for (rot in rotations) {
      
      rot_df <- df[which(df$rotation == rot),]
      
      points(x = rot_df$response[which(rot_df$target=='point')],
             y = rot_df$response[which(rot_df$target=='arc')],
             pch = 16,
             col = colsT[rot],
             cex=1.4)
      
      
    }
    
    dot_df <- aggregate(response ~ participant + rotation, data = df[which(df$target=='point'),], FUN=mean)
    names(dot_df)[which(names(dot_df) == 'response')] <- 'dot'
    arc_df <- aggregate(response ~ participant + rotation, data = df[which(df$target=='arc'),], FUN=mean)
    names(arc_df)[which(names(arc_df) == 'response')] <- 'arc'
    
    dotarc <- merge(dot_df,
                    arc_df,
                    by = c('participant','rotation'))
    

    rot_lm <- lm(arc ~ dot, data=dotarc)
    
    intercept <- rot_lm$coefficients[1]
    slope     <- rot_lm$coefficients[2]
    dotrange <- range(dotarc$dot)
    
    correlation <- cor.test(dotarc$dot, dotarc$arc)
    # print(names(correlation))
    r.stat <- correlation$estimate # the r
    
    text(18,-10,sprintf('r = %0.3f',r.stat),adj=1)
    
    
    odr <- Reach::odregress(x=dotarc$dot,
                            y=dotarc$arc)
    
    
    axis(side=1,at=c(-12,-6,0,6,12,18),cex.axis=0.75)
    axis(side=2,at=c(-12,-6,0,6,12,18),cex.axis=0.75)
    
    if (maxrot == 60) {
      x = seq(-12,-3,len=61*4)
      y = c(15,16)
      z = matrix(seq(0,60,len=61*4),ncol=1)
      image(x,y,z, add=TRUE, col=colsO)
      text(x = c(-12,-7.5,-3),
           y = c(14,14,14),
           labels = sprintf('%d°',c(0,30,60)))
      text(x = c(-7.5),
           y = c(18),
           labels = c('rotation'))
    }
    
  }
  
  
  
  if (target %in% c('pdf','svg','png','tiff')) {
    dev.off()
  }
  
}


plotAverageIRDs <- function(target='inline',variables='og') {
  
  width = 8
  if (target %in% c('pdf','svg','png','tiff')) {
    setupFigureFile( target   = target,
                     width    = width,
                     height   = 4,
                     dpi      = 300,
                     filename = sprintf('doc/fig3+_average_rot_IRDS.%s',target))
  }
  
  left  <- 0.75
  right <- 0.05
  
  deg15 <- (width - ((left + right) * 2)) / 7
  
  layout(mat=matrix(c(1:2), ncol=2))
  par(mai=c(0.75,left,0.4,right))
  
  # set up color map:
  
  df <- loadSTLdata(average=median, maxrots=c(45,60))
  all_rots <- sort(unique(df$rotation))
  
  # colT <- Reach::colorMix('#FF000033','#0000FF33',c(2,3))
  # colO <- Reach::colorMix('#FF0000','#0000FF',c(2,4))
  
  colors <- c('#9900FFFF', '#9900FF33')
  
  for (maxrot in c(45,60)) {
    
    plot(x=-1000,y=-1000,
         main=sprintf('%d° max rotation',maxrot),ylab='',xlab='',
         # xlim=c(0,maxrot+1),
         xlim=c(-3,9),
         ylim=c(-3,9),
         ax=F,bty='n',asp=1)
    
    if (variables == 'og') {
      title(ylab='arc initial reach deviation [°]',line=2.25)
      title(xlab='dot initial reach deviation [°]',line=2.25)
    }
    if (variables == 'flipped') {
      title(ylab='dot initial reach deviation [°]',line=2.25)
      title(xlab='arc initial reach deviation [°]',line=2.25)
    }
    
    lines(x=c(-3,9),
          y=c(0,0),
          lty=1,col='#999999')
    lines(x=c(0,0),
          y=c(-3,9),
          lty=1,col='#999999')
    lines(x=c(-3,9),
          y=c(-3,9),
          lty=1,col='#999999')
    
    df <- loadSTLdata(average=median, maxrots=c(maxrot))
    
    
    rot_df <- aggregate(response ~ participant + target, data=df, FUN=mean)
    
    if (variables == 'og') {
      points(x = rot_df$response[which(rot_df$target=='point')],
             y = rot_df$response[which(rot_df$target=='arc')],
             pch = 1,
             col = colors[1],
             cex=1.4)
    }
    if (variables == 'flipped') {
      points(x = rot_df$response[which(rot_df$target=='arc')],
             y = rot_df$response[which(rot_df$target=='point')],
             pch = 1,
             col = colors[1],
             cex=1.4)
    }
    cat(sprintf("%s group dot & point correlation\n", maxrot))
    
    lm_df <- rot_df[which(rot_df$target == 'point'),]
    lm_df$target <- NULL # remove the column
    names(lm_df) <- c('participant', 'point')
    lm_df$arc <- rot_df$response[which(rot_df$target == 'arc')]
    
    if (variables == 'og') {
      lmod <- lm(arc ~ point, data=lm_df)
    }
    if (variables == 'flipped') {
      lmod <- lm(point ~ arc, data=lm_df)
    }
    print(summary(lmod))
    
    coef <- lmod$coef
    if (variables == 'og') {
      Xrange <- range(rot_df$response[which(rot_df$target=='point')])
    }
    if (variables == 'flipped') {
      Xrange <- range(rot_df$response[which(rot_df$target=='arc')])
    }
    
    X <- seq(Xrange[1],Xrange[2],len=100)
    
    if (variables == 'og') {
      predlm <- predict(lmod, newdata=data.frame(point=X), se.fit = TRUE)
    }
    if (variables == 'flipped') {
      predlm <- predict(lmod, newdata=data.frame(arc=X), se.fit = TRUE)
    }
    
    polygon(x=c(X,rev(X)),
            y=c(predlm$fit + 1.96*predlm$se.fit, rev(predlm$fit - 1.96*predlm$se.fit)),
            col=colors[2],
            border=NA)
    lines(x=X,
          y=predlm$fit,
          col=colors[1])
    
    axis(side=1,at=c(-3,0,3,6,9),cex.axis=0.75)
    axis(side=2,at=c(-3,0,3,6,9),cex.axis=0.75)
    
    
  }
  
  
  
  if (target %in% c('pdf','svg','png','tiff')) {
    dev.off()
  }
  
}




# initial reach deviation over rotation -----

plotData <- function(models=FALSE, target='inline') {
  
  if (models) {
    filename=sprintf('doc/fig4_model_fits.%s',target)
  } else {
    filename=sprintf('doc/fig2_alldata.%s',target)
  }
  
  
  if (target %in% c('pdf','svg','png','tiff')) {
    setupFigureFile( target   = target,
                     width    = 6,
                     height   = 7,
                     dpi      = 300,
                     filename = filename)
  }
  

  layout(mat=matrix(c(1:3), ncol=1))
  par(mar=c(3.25,3.25,2,0.1))
  
  for (maxrot in c(45,60,90)) {
    
    plot(x=-1000,y=-1000,
         main=sprintf('%d° max rotation',maxrot),ylab='',xlab='',
         # xlim=c(0,maxrot+1),
         xlim=c(0,91),
         ylim=c(-0.5,8.5),
         ax=F,bty='n')
    
    title(ylab='reach aftereffect [°]',line=2.25)
    title(xlab='rotation [°]',line=2.25)
    
    
    df <- loadSTLdata(average=median, maxrots=c(maxrot))
    
    # print(sprintf('data/modelFits_%d.csv', maxrot))
    
    if (models) {
      fits <- read.csv(sprintf('data/modelFits_%d.csv', maxrot), stringsAsFactors = FALSE)
      fits <- fits[which(fits$participant=='all'),]
    }
    
    CIdf <- aggregate(response ~ target + rotation, data=df, FUN=Reach::getConfidenceInterval, method='b', resamples=2000)
    avgdf <- aggregate(response ~ target + rotation, data=df, FUN=mean)
    for (tt in c('point','arc')) {
      
      if (tt == 'point') {
        colors <- c('#0066FFFF', '#0066FF33')
      }
      if (tt == 'arc') {
        colors <- c('#FF6600FF', '#FF660033')
      }
      
      tCI <- CIdf[which(CIdf$target == tt),]
      tavg <- avgdf[which(avgdf$target == tt),]
      
      X <- c(tCI$rotation, rev(tCI$rotation))
      Y <- c(tCI$response[,1], rev(tCI$response[,2]))
      polygon(X,Y,border=NA,col=colors[2])
      
      if (models) {
        fit <- fits[which(fits$target == tt),]
        # cross <- fit$c/fit$r
        
        lines(x=c(0,fit$c/fit$r,maxrot),
              y=c(0,fit$c,fit$c),
              lty=3,
              col=colors[1],
              lw=2)
        
        attpar = c('s'=fit$s, 'w'=fit$w)
        rots = seq(0,maxrot,0.5)
        attribution_predictions <- STLPredict(par=attpar,
                                              rotations=rots)
        
        lines(x=rots,
              y=attribution_predictions,
              lty=2,
              col=colors[1],
              lw=2)
      } else {
        lines(x=tavg$rotation, y=tavg$response, col=colors[1], lw=2)
      }
      
    }
    
    if (maxrot == 45) {
      axis(side=1,at=c(1,5,10,15,20,25,30,35,40,45),cex.axis=0.75)
    }
    if (maxrot == 60) {
      axis(side=1,at=c(1,5,10,15,20,25,30,40,50,60),cex.axis=0.75)
    }
    if (maxrot == 90) {
      axis(side=1,at=c(1,5,10,15,20,30,40,50,70,90),cex.axis=0.75)
    }
    axis(side=2,at=c(0,2,4,6,8),cex.axis=0.75)
    
    
    if (maxrot == 45) {
      legend(x=20,y=2,
             legend=c('point target', 'arc target'),
             lwd=c(2,2),
             col=c('#0066FFFF','#FF6600FF'),
             bty='n')
    }
    if (maxrot == 60) {
      if (models) {
        legend(x=20,y=2,
               legend=c('attribution model', 'capped model'),
               lwd=c(2,2),
               lty=c(2,3),
               col=c('#999999','#999999'),
               bty='n')
      }
    }
    
  }
  
  if (target %in% c('pdf','svg','png','tiff')) {
    dev.off()
  }
  
}



plotModelMSEs <- function(maxrots=c(45,60,90), target='inline') {
  
  if (target %in% c('pdf','svg','png','tiff')) {
    setupFigureFile( target   = target,
                     width    = 4,
                     height   = 4,
                     dpi      = 300,
                     filename = sprintf('doc/fig5_mses.%s',target))
  }
  
  df <- NA
  for (maxrot in maxrots) {
    mf <- read.csv(sprintf('data/modelFits_%d.csv', maxrot))
    mf <- mf[which(mf$participant %notin% c('all')),]
    mf$rot <- maxrot
    if (is.data.frame(df)) {
      df <- rbind(df, mf)
    } else {
      df <- mf
    }
  }
  
  plot(-1000,-1000,
       main='model MSEs',xlab='',ylab='',
       xlim = c(78,112),
       ylim = c(78,112),
       ax=F, bty='n',
       asp=1)
  
  title(ylab='capped MSE',       line=2.5)
  title(xlab='attribution MSE',  line=2.5)
  
  lines( x = c(78,112),
         y = c(78,112),
         col='#99999999')
  
  for (rot in unique(df$rot)) {
    # print(rot)
    col <- c('45'='#FF000033', '60'='#9900FF33', '90'='#0099FF33')[sprintf('%d',rot)]
    idx <- which(df$rot == rot)
    points( x = df$attributionMSE[idx], 
            y = df$cappedMSE[idx],
            pch=16, cex=1.5,
            col=col)
    # print(c(df$attributionMSE[idx],df$cappedMSE[idx]))
  }
  
  axis(side=1,at=seq(80,110,10))
  axis(side=2,at=seq(80,110,10))
  
  legend( x=80,
          y=110,
          legend=c('45° group', '60° group', '90° group'),
          col=c('#FF000099', '#9900FF99','#0099FF'),
          pch=c(16,16),
          bty='n')
  
  if (target %in% c('pdf','svg','png','tiff')) {
    dev.off()
  }
  
}

plotLearningPrediction <- function(target='inline') {
  
  if (target %in% c('pdf','svg','png','tiff')) {
    setupFigureFile( target   = target,
                     width    = 7,
                     height   = 8,
                     dpi      = 300,
                     filename = sprintf('doc/fig6_prediction.%s',target))
  }
  
  layout(mat=matrix(c(1:9),ncol=3,byrow=TRUE),widths=c(1.9,1,1))
  
  par(mar=c(3.25,3.25,3.25,0.1))
  
  timecourses <- getLongTimeCourses()
  
  predictions <- getPredictionData()
  
  for (maxrot in c(45,60,90)) {
    
    plot(-1000,-1000,
         main=sprintf('%d° group',maxrot), xlab='', ylab='',
         xlim=c(0,161),ylim=c(-5,25),
         bty='n',ax=F)
    title(xlab='trial',line=2.25)
    title(ylab='reach deviation [°]',line=2.25)
    
    colt <- c('45'='#FF000033', '60'='#9900FF33', '90'='#0099FF33')[sprintf('%d',maxrot)]
    colo <- c('45'='#FF0000ff', '60'='#9900FFff', '90'='#0099FFFF')[sprintf('%d',maxrot)]
    
    df <- timecourses[[sprintf('%d',maxrot)]]
    participants <- unique(df$participant)
    
    CIdf <- aggregate(reachdeviation_deg ~ trial, data=df, FUN=Reach::getConfidenceInterval, method='t', resamples=1000)
    
    X <- c(c(1:160), rev(c(1:160)))
    Y <- c(CIdf$reachdeviation_deg[,1], rev(CIdf$reachdeviation_deg[,2]))
    polygon(X,Y,border=NA,col=colt)
    
    avg <- aggregate(reachdeviation_deg ~ trial, data=df, FUN=mean, na.rm=TRUE)
    lines(avg, col=colo)
    
    axis(side=1,at=c(1,40,80,120,160),cex.axis=0.8)
    axis(side=2,at=c(0,10,20),cex.axis=0.8)
    
    preds <- predictions[[sprintf('%d',maxrot)]]
    ppreds <- preds[which(preds$target == 'point'),]
    
    for (model in c('capped','attribution')) {
      
      plot(-1000,-1000,
           main=sprintf('%d° group',maxrot), xlab='', ylab='',
           xlim=c(0,15),ylim=c(0,15),
           bty='n',ax=F,asp=1)
      title(xlab=sprintf('%s prediction [°]',model),line=2.25)
      title(ylab='exponential fit [°]',line=2.25)
      
      lines(x=c(0,15),y=c(0,15),col='#99999999')
      
      X = ppreds[,model]
      Y = ppreds$exponential
      
      # plot data points:
      points( x = X,
              y = Y,
              pch=16, cex=1.5,
              col=colt)
      
      # plot linear fit with intercept=0
      linmod <- lm(Y~X+0)
      coeff <- linmod$coeff['X']
      lines( x=c(0,15),
             y=c(0,15*coeff),
             col=colo)
      
      axis(side=1,at=seq(0,15,5),cex.axis=0.8)
      axis(side=2,at=seq(0,15,5),cex.axis=0.8)
      
    }
    
  }
  
  if (target %in% c('pdf','svg','png','tiff')) {
    dev.off()
  }
  
}

plotFirstTrialFits <- function(target='inline') {
  
  if (target %in% c('pdf','svg','png','tiff')) {
    setupFigureFile( target   = target,
                     width    = 7,
                     height   = 4,
                     dpi      = 300,
                     filename = sprintf('doc/fig7_exp1st.%s',target))
  }
  
  layout(mat=matrix(c(1:3),ncol=3,byrow=TRUE))
  par(mar=c(3.25,3.25,3.25,0.1))

  firstTrials <- getFirstTrialData()
  
  for (maxrot in c(45,60,90)) {
    
    plot(-1000,-1000,
         main='',xlab='',ylab='',
         xlim=c(0.5,3.5),ylim=c(-35,25),
         bty='n',ax=F)
    
    lines( x=c(0.5,3.5), y=c(0,0), lty=2, col='#999999' )
    lines( x=c(0.5,3.5), y=c(20,20), lty=2, col='#999999' )
    
    title(main=sprintf('second trials - %d°', maxrot))
    title(ylab='(reach) deviation [°]', line=2.25)
    df <- firstTrials[which(firstTrials$maxrot == maxrot),]
    
    col <- c('45'='#FF000033', '60'='#9900FF33', '90'='#0099FF33')[sprintf('%d',maxrot)]
    scol <- c('45'='#FF0000', '60'='#9900FF', '90'='#0099FF')[sprintf('%d',maxrot)]
    
    for (depvar in c(1,2,3)) {
      
      if (depvar == 1) { vals <- df$reachdeviation_deg }
      if (depvar == 2) { vals <- df$prediction }
      if (depvar == 3) { vals <- df$prediction-df$reachdeviation_deg  }

      points(x = seq(depvar-0.3, depvar-0.1, length.out=dim(df)[1]),
             # x=rep(depvar-0.3,dim(df)[1]),
             y=vals,
             col=col)
      
      CI <- Reach::getConfidenceInterval(vals,method='b')
      
      polygon( x = c(depvar+0.1,depvar+0.3,depvar+0.3,depvar+0.1),
               y = c(CI[1],CI[1],CI[2],CI[2]),
               col=col,
               border=NA)
      
      lines(x=c(depvar+0.1,depvar+0.3),
            y=rep(mean(vals),2),
            col=scol)
      
    }
    
    axis(side=1,at=c(1,2,3),labels=c('reach\ndata', 'exponential\nprediction', 'difference\n '), padj=0.5, cex.axis=0.75)
    axis(side=2,at=seq(-20,20,10), cex.axis=0.75)
    

    # print(range(df$reachdeviation_deg))
    
    # print(t.test(df$prediction, df$reachdeviation_deg, paired=TRUE))
    
  }
  
  if (target %in% c('pdf','svg','png','tiff')) {
    dev.off()
  }
  
}

plotTimeOnTask <- function(target='inline') {
  
  if (target %in% c('pdf','svg','png','tiff')) {
    setupFigureFile( target   = target,
                     width    = 7,
                     height   = 7,
                     dpi      = 300,
                     filename = sprintf('doc/fig8_ToT.%s',target))
  }
  
  layout(mat=matrix(c(1:4),ncol=2,byrow=TRUE))
  par(mar=c(3.25,3.25,3.25,0.1))
  
  for (tt in c('point','arc')) {
  
    for (maxrot in c(45,60)) {
      
      
      col <- c('45'='#FF000033', '60'='#9900FF33')[sprintf('%d',maxrot)]
      scol <- c('45'='#FF0000', '60'='#9900FF')[sprintf('%d',maxrot)]
      
      # if (maxrot==45) {
      #   cols <- c('#FF0000','#FF3300','#FF6600','#FF9900','#FFCC00')
      # }
      # if (maxrot==60) {
      #   cols <- c('#9900FF','#9933FF','#9966FF','#9999FF','#99CCFF')
      # }
      
      cols <- c('#9900FF','#9933FF','#9966FF','#9999FF','#99CCFF')
      
      plot(-1000,-1000,
           main='', xlab='', ylab='',
           xlim=c(0,maxrot+6),ylim=c(-1,9),
           bty='n',ax=F)
      title(main=sprintf('%d %s',maxrot, tt))
      title(ylab='reach aftereffect [°]',line=2.25)
      title(xlab='rotation [°]',line=2.25)
      
      lines(x=c(1,maxrot),y=c(0,0),col='#999',lty=2)
      
      sb_avg <- c()
      
      for (superblock in c(1,2,3,4,5)) {
        
        df <- loadSTLdata(targets=c(tt),
                          average=median, 
                          maxrots=c(maxrot), 
                          superblocks=c(superblock))
        
        avg <- aggregate(response ~ rotation, data=df, FUN=mean, na.rm=TRUE)
        
        sb_avg <- c(sb_avg, mean(avg$response))
        
        lines(x=avg$rotation, y=avg$response, col=cols[superblock])
        
      }
      
      points(rep(maxrot+5,5), sb_avg, col=cols)
      
      if (maxrot == 45) {
        axis(side=1, at=c(1,5,10,15,20,25,30,35,40,45), cex.axis=0.8)
      }
      if (maxrot == 60) {
        axis(side=1, at=c(1,5,10,15,20,25,30,40,50,60), cex.axis=0.8)
      }
      axis(side=2, at=c(0,2,4,6,8), cex.axis=0.8)
      
      if (tt == 'point' & maxrot==60) {
        legend(x=40,y=3.5,
               legend=sprintf('block %d',c(1:5)),
               col=cols,lty=1,seg.len=2,bty='n',cex=0.7)
      }
      
    }
    
  }
  
  if (target %in% c('pdf','png','tiff','svg')) {
    dev.off()
  }
  
}