# Presentation-only graphics. Inputs contain retained selections and contributions.
args <- commandArgs(trailingOnly=TRUE)
streamed <- length(args) >= 3 && args[[3]] == '--jsonl'
records <- if(streamed) NULL else jsonlite::fromJSON(args[[1]], simplifyVector=FALSE)
out <- args[[2]]
wrap <- function(x, width=44) {
  lines <- strwrap(x,width=width)
  paste(unlist(lapply(lines,function(line) {
    if(nchar(line)<=width) line else substring(line,seq(1,nchar(line),width),seq(width,nchar(line)+width-1,width))
  })),collapse='\n')
}
value <- function(x, missing='NA') if (is.null(x)) missing else as.character(x)
draw <- function(record) {
  rows <- record$rows
  n <- max(1, length(rows))
  par(mar=c(4, 1, 3, 1), family='sans', xpd=NA)
  if (record$kind == 'tissues') {
    layout(matrix(seq_along(record$panels),ncol=1))
    for(panel in record$panels) {
      rows <- panel$rows; tissues <- unlist(panel$tissues); n <- length(rows); m <- length(tissues)
      par(mar=c(8,19,4,2),family='sans',xpd=NA)
      plot.new(); plot.window(xlim=c(0,max(1,m)),ylim=c(0,max(1,n)))
      if(m) axis(1,at=seq_len(m)-.5,labels=vapply(tissues,wrap,'',width=20),las=2,tick=FALSE,cex.axis=.8)
      for(i in seq_along(rows)) {
        r <- rows[[i]]; y <- n-i
        score <- if(r$N>0) sprintf('%.3f',r$numerator/r$denominator) else 'NA'
        label <- paste(wrap(value(r$name,r$medication_entity_key),28),paste0('S_L=',score,'  P=',r$P,' C=',r$C,' N=',r$N),sep='\n')
        text(-.08,y+.5,label,adj=1,cex=.75)
        for(j in seq_along(tissues)) {
          matches <- Filter(function(x) x$tissue_id==tissues[[j]],r$coverage)
          state <- if(length(matches)) matches[[1]]$state else 'unavailable'
          fill <- switch(state,support='#A6CEE3',conflict='#FDBF6F',unavailable='#D9D9D9','white')
          mark <- switch(state,support='+1',conflict='-1',mixed='0',no_two_direction_evidence='0',pending='Pending','NA')
          rect(j-1,y,j,y+1,col=fill,border='#444444',lty=if(state=='pending') 2 else 1)
          if(state=='mixed') rect(j-1,y,j,y+1,density=14,angle=45,col='#777777',border='#444444')
          text(j-.5,y+.5,mark,cex=if(state=='pending') .65 else .85,font=2)
        }
      }
      title(paste('Tissue expression evidence',panel$gene_id,wrap(panel$scenario_id,80),sep='\n'),cex.main=.8)
      mtext('Support +1 | Conflict -1 | Patterned 0: both | Plain 0: neither | NA: unavailable | Pending: incomplete',side=1,line=6.5,cex=.65)
    }
    layout(1)
  } else if(record$kind == 'tissue_conditions') {
    par(mar=c(5,26,4,2),family='sans',xpd=NA)
    plot.new();plot.window(xlim=c(0,2),ylim=c(0,n))
    axis(1,at=c(.5,1.5),labels=c('rOE','-rKD'),tick=FALSE)
    for(i in seq_along(rows)) {
      r <- rows[[i]]; y <- n-i
      witnesses <- if(length(r$directional_witnesses)) paste(gsub('_',' ',unlist(r$directional_witnesses)),collapse=', ') else 'none'
      label <- paste(wrap(r$pert_iname,40),wrap(paste(r$cell_id,value(r$pert_dose),value(r$pert_dose_unit),
          r$pert_time,r$pert_time_unit,sep=' | '),60),wrap(r$sig_id,60),wrap(paste('Witnesses:',witnesses),60),sep='\n')
      text(-.06,y+.5,label,adj=1,cex=.7)
      vals <- list(r$r_oe,r$kd_reversal)
      for(j in 1:2) {
        v <- vals[[j]]
        fill <- if(is.null(v)) '#D9D9D9' else if(v>0) '#A6CEE3' else if(v<0) '#FDBF6F' else 'white'
        rect(j-1,y,j,y+1,col=fill,border='#444444')
        text(j-.5,y+.5,if(is.null(v)) 'NA' else if(v==0) '0.000' else sprintf('%+.3f',v),cex=.9)
      }
    }
    title(paste('Individual drug conditions',record$gene_id,record$tissue_id,sep=' / '),cex.main=.85)
    mtext('Positive: restoration-direction agreement | Negative: opposite direction | 0: neither sign | NA: undefined or unavailable',side=1,line=3,cex=.65)
  } else if(record$kind == 'tissue_contributions') {
    par(mar=c(6,15,3,2),family='sans',xpd=NA)
    lower <- min(0,vapply(rows,function(r) if(is.null(r$scores$S_L)) 0 else r$scores$S_L,0))
    upper <- max(1,vapply(rows,function(r) max(0,if(is.null(r$scores$S_L)) 0 else r$scores$S_L)+r$scores$S_M,0))
    plot.new();plot.window(xlim=c(lower-.2,upper+1),ylim=c(.3,n+.7));axis(1);abline(v=0,col='#444444',xpd=FALSE)
    for(i in seq_along(rows)) {
      r<-rows[[i]];y<-n-i+1;v<-if(is.null(r$scores$S_L)) 0 else r$scores$S_L
      text(lower-.25,y,wrap(r$name,25),adj=1,cex=.75)
      rect(min(0,v),y-.25,max(0,v),y+.25,col=if(v>=0) '#A6CEE3' else '#FDBF6F',border='#444444')
      x<-max(0,v);rect(x,y-.25,x+r$scores$S_M,y+.25,col='#B2DF8A',border='#444444')
      text(upper+.9,y,paste0('S=',sprintf('%.3f',r$scores$S)),adj=1,cex=.75)
      text(.05,y+.35,paste0('S_L=',if(is.null(r$scores$S_L)) 'NA' else sprintf('%.3f',v),'  S_M=',r$scores$S_M),adj=0,cex=.65)
    }
    title('Figure 2. Tissue and mechanism contributions',cex.main=1)
    mtext('Descriptive score points',side=1,line=2.5,cex=.8)
    mtext('Blue/orange: signed S_L | Green: S_M | Unavailable S_L contributes zero only to S',side=1,line=4.5,cex=.7)
  } else if (record$kind == 'context') {
    plot.new(); plot.window(xlim=c(0, 4), ylim=c(0, n))
    for (i in seq_along(rows)) {
      r <- rows[[i]]; y <- n-i
      labels <- c(paste(r$hypothesis, paste(vapply(r$variants, function(v) paste(v$label, v$description), ''),collapse='\n'),sep='\n'),
        paste('Interpreted effect',r$effect,'Evidence status',r$evidence_status,sep='\n'),
        paste('Affected process',r$process,sep='\n'),
        paste('Objective',if(length(r$objectives)) paste(vapply(r$objectives,function(o) o$desired_modulation,''),collapse='\n') else 'Unavailable',sep='\n'))
      for (j in seq_len(4)) {
        rect(j-1+.06,y+.06,j-.06,y+.94,col='#F4F6F7',border='#263238')
        text(j-.5,y+.5,wrap(labels[[j]],30),cex=.65)
        if(j<4) arrows(j-.05,y+.5,j+.05,y+.5,length=.06)
      }
    }
    title('Figure 1. Retained variant hypotheses',cex.main=1)
  } else if (record$kind == 'contributions') {
    terms <- c('I_L','2B','0.5U','I_M','2M','-0.5C')
    colors <- c('#0072B2','#56B4E9','#009E73','#E69F00','#F0E442','#D55E00')
    contributions <- lapply(rows,function(r) c(r$components$I_L,2*r$components$B,.5*r$components$U,r$components$I_M,2*r$components$M,-.5*r$components$C))
    lower <- min(0,vapply(contributions,function(x) x[[6]],0))
    upper <- max(1,vapply(contributions,function(x) sum(x[1:5]),0))
    par(mar=c(5,14,3,2))
    plot.new(); plot.window(xlim=c(lower-1,upper+1),ylim=c(.3,n+.7))
    axis(1); abline(v=0,col='#333333',xpd=FALSE)
    axis(2,at=seq_len(n),labels=vapply(rev(rows),function(r) wrap(r$name,25),''),las=1,tick=FALSE,cex.axis=.8)
    for (i in seq_along(rows)) {
      y <- n-i+1; vals <- contributions[[i]]; x <- 0
      for (j in seq_len(5)) {
        if (vals[[j]]>0) {
          rect(x,y-.26,x+vals[[j]],y+.26,col=colors[[j]],border='#444444')
          text(x+vals[[j]]/2,y,vals[[j]],cex=.7);x <- x+vals[[j]]
        }
      }
      if(vals[[6]]<0) {
        rect(vals[[6]],y-.26,0,y+.26,col=colors[[6]],border='#444444')
        text(vals[[6]]/2,y,vals[[6]],cex=.7)
      }
      text(upper+.8,y,paste('S =',rows[[i]]$scores$S),adj=1,cex=.7)
    }
    title('Figure 2. Combined-score contributions',xlab='Descriptive score points',cex.main=1)
    legend('bottom',inset=c(0,-.32),legend=terms,fill=colors,horiz=TRUE,bty='n',cex=.7)
  } else {
    plot.new(); plot.window(xlim=c(0,1),ylim=c(0,n+1))
    text(.02,n+.7,'Medication / cell / method / query size',adj=0,cex=.8)
    text(.92,n+.7,'State',cex=.8)
    for(i in seq_along(rows)) {
      r <- rows[[i]]; y <- n-i+.5
      mark <- paste(unlist(r$symbols),collapse=' + ')
      if(isTRUE(r$conflict)) mark <- paste(mark,'!')
      color <- if(r$state=='pass') '#D6EAF8' else if(r$state=='fail') '#FADBD8' else if(r$state=='mixed') '#FFF2CC' else '#EEEEEE'
      rect(.01,y-.42,.99,y+.42,col=color,border='#999999')
      text(.02,y,paste(r$name,r$cell_id,r$method,paste0('n=',value(r$query_gene_count)),sep=' / '),adj=0,cex=.7)
      text(.92,y,mark,cex=.8)
    }
    title('Figure 3. Cell agreement and missing evidence',cex.main=1)
  }
}
render <- function(record) {
  height <- if(record$kind=='context') max(3,2.3*length(record$rows)) else if(record$kind=='tissues')
      sum(vapply(record$panels,function(p) max(4,.7*length(p$rows)+3.4),0)) else if(record$kind=='tissue_conditions')
      max(4,1.05*length(record$rows)+2.5) else max(3.5,.65*length(record$rows)+1.5)
  width <- if(record$kind=='tissues') max(12,6+.7*max(vapply(record$panels,function(p) length(p$tissues),0))) else 12
  for(format in c('png','svg')) {
    path <- file.path(out,paste0(record$id,'.',format))
    if(format=='png') png(path,width=width,height=height,units='in',res=120) else svg(path,width=width,height=height)
    draw(record)
    dev.off()
  }
}
if(streamed) {
  connection <- file(args[[1]],'r')
  repeat {
    line <- readLines(connection,n=1,warn=FALSE)
    if(!length(line)) break
    render(jsonlite::fromJSON(line,simplifyVector=FALSE))
  }
  close(connection)
} else for(record in records) render(record)
