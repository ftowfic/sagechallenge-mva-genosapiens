# Deterministic scientific plots from the same frozen records used in the tables.
suppressPackageStartupMessages(library(jsonlite))
args <- commandArgs(trailingOnly=TRUE)
streamed <- length(args) >= 3 && args[3] == '--jsonl'
if (!streamed) figures <- fromJSON(args[1], simplifyVector=FALSE)
out <- args[2]
draw <- function(f) {
  par(mar=c(4,4,4,2), family='sans', xpd=FALSE)
  if (f$kind == 'selection') {
    q <- f$query; g <- q$genes
    selected <- Filter(function(x) isTRUE(x$selected),g)
    up <- sum(vapply(selected,function(x) x$desired_direction=='up',logical(1)))
    if (!is.null(q$method)) {
      plot.new();plot.window(xlim=c(0,10),ylim=c(0,10))
      title(main=paste(q$gene_id,q$cell_id,sep=' | '))
      text(5,8,q$method,cex=1.2)
      text(5,6,paste(length(g),'measured features.',length(selected),'selected features'))
      text(5,4,paste('Desired direction:',up,'up.',length(selected)-up,'down'))
      text(5,2,if(is.null(q$absent_direction)) 'Both genetic directions available' else paste('Missing genetic direction:',q$absent_direction))
      text(5,.5,'Expression matching does not establish rescue.',cex=.9)
      return(invisible(NULL))
    }
    plot.new();plot.window(xlim=c(0,10),ylim=c(0,10))
    title(main=paste(f$cell_id, '| How this response pattern was chosen'),cex.main=1.3)
    boxes <- list(c(1,7.4,9,9),c(1,4.5,4.6,6.5),c(5.4,4.5,9,6.5),c(1,1,9,3.2))
    cols <- c('#E8F1F8','#DBE9F5','#FBE9D0','#E9F2EB')
    for(i in seq_along(boxes)) {b<-boxes[[i]];rect(b[1],b[2],b[3],b[4],col=cols[i],border='#526170',lwd=1.5)}
    text(5,8.25,paste(q$target_id,':',q$required_modulation,'function\nRestoration objective. Drug rescue remains unproven.'),cex=1.05)
    text(2.8,5.5,'KD: less target expression\nLocal shRNA consensus',cex=1.05)
    text(7.2,5.5,'OE: extra target expression\nLocal overexpression',cex=1.05)
    arrows(c(3,7),7.4,c(2.8,7.2),6.6,length=.12)
    arrows(c(2.8,7.2),4.5,c(4,6),3.3,length=.12)
    text(5,2.1,paste0(length(g),' measured landmarks -> ',length(selected),' selected\n',
       'Opposite signs AND |MODZ| >= 1 in both\nDesired: ',up,' up / ',length(selected)-up,' down. Minimum: 3 usable genes'),cex=1.05)
    text(5,.3,paste('Query:',q$state,if(!is.null(q$reason)) q$reason else '| Selected opposition does not establish validation'),cex=.85)
  } else if (f$kind == 'heatmap') {
    genes<-f$genes;conditions<-f$conditions;ng<-length(genes);nc<-length(conditions)+2
    labels<-c('KD\nlocal genetic','OE\nlocal genetic',vapply(conditions,function(x)
      paste(x$pert_iname,paste(x$pert_dose,x$pert_dose_unit),paste(x$pert_time,'h'),
            sub('_.*','',x$sig_id),sep='\n'),character(1)))
    par(mar=c(8,8,6,4))
    plot.new();plot.window(xlim=c(-.6,nc+1),ylim=c(.4,ng+.6),xaxs='i',yaxs='i')
    pal<-colorRampPalette(unlist(f$palette))(257)
    for(i in seq_len(ng)) {
      y<-ng-i+1;g<-genes[[i]]
      vals<-c(if(is.null(g$kd)) NA_real_ else g$kd,if(is.null(g$oe)) NA_real_ else g$oe,vapply(conditions,function(x) {
        v<-x$values[[f$gene_offset+i]];if(is.null(v)) NA_real_ else v
      },numeric(1)))
      text(0,y,if(g$desired_direction=='up') '\u2191' else '\u2193',cex=1.2)
      for(j in seq_len(nc)) {
        v<-vals[j];color<-if(is.na(v)) '#BDBDBD' else pal[round((v/f$expression_limit+1)*128)+1]
        rect(j-.47,y-.47,j+.47,y+.47,col=color,border='white')
        if(is.na(v))text(j,y,'NA',cex=.7)
      }
    }
    axis(2,at=seq_len(ng),labels=rev(vapply(genes,function(g)paste(g$gene_symbol,g$gene_id,sep=' / '),character(1))),las=1,tick=FALSE,cex.axis=.85)
    text(seq_len(nc),par('usr')[3]-.03*diff(par('usr')[3:4]),labels,
         adj=c(.5,1),xpd=NA,cex=.65)
    text(0,par('usr')[3]-.03*diff(par('usr')[3:4]),'Desired\ndirection',
         adj=c(.5,1),xpd=NA,cex=.7)
    title(main=paste(f$cell_id,'| Measured gene responses'),cex.main=1.3)
    mtext(c(sprintf('Lower: %.2f',-f$expression_limit),'Neutral: 0',
            sprintf('Higher: +%.2f MODZ',f$expression_limit),'Missing: NA'),side=3,line=1,
          at=seq(.5,nc,length.out=4),col=c(pal[1],'#526170',pal[257],'#737373'),cex=.8)
    mtext('Shared zero-centred MODZ scale. No clipping. Arrows are categorical, not measured values.',side=1,line=6.4,cex=.8)
  } else {
    rows<-f$medications
    if(!length(rows)) {
      plot.new();title(main=paste(f$cell_id,'| No local nomination'))
      text(.5,.6,paste(f$n_tested,'tested medications.',f$n_fail,'failed.',f$n_unassessed,'unassessed'),cex=1.2)
      text(.5,.4,'No winning drug is inferred from missing or failed local evidence.\nSee full query, exclusion and condition tables.',cex=1)
    } else {
      par(mar=c(5,15,6,2))
      n<-length(rows);plot.new();plot.window(xlim=c(-1.05,1.30),ylim=c(.4,n+.6))
      abline(v=c(-1,-.5,0,.5,1),col='#DDDDDD');abline(v=0,col='#526170',lwd=1.4)
      colors<-c('#2166AC','#D97706','#23836B')
      for(i in seq_len(n)) {
        r<-rows[[i]];y<-n-i+1
        paired<-is.null(r$method)||r$method=='paired_oe_kd'
        vals<-if(paired) c(r$r_desired,-r$r_opposite,r$joint) else r$r_desired
        yy<-if(paired) y+c(.22,0,-.22) else y
        segments(0,yy,vals,yy,col=colors[seq_along(vals)],lwd=4);points(vals,yy,pch=c(16,17,15)[seq_along(vals)],col=colors[seq_along(vals)],cex=1.2)
        text(1.15,yy,sprintf('%+.3f',vals),cex=.75,col=colors[seq_along(vals)])
        label<-paste0('#',r$rank,' ',r$name,'\n',r$n_agree,'/',r$n_conditions,' conditions, ',length(r$times_hours),' time(s)',
                      if(isTRUE(r$conflict)) ' | CONFLICT' else '', if(paired) paste0('\nrOpp=',sprintf('%+.3f',r$r_opposite)) else '\nOpposite direction unavailable')
        axis(2,at=y,labels=label,las=1,tick=FALSE,cex.axis=.75)
      }
      axis(1,at=seq(-1,1,.5))
      axis_label<-if(is.null(f$method)||f$method=='paired_oe_kd') 'Signed correlation / joint summary (-1 to +1)' else 'Desired Spearman correlation (-1 to +1)'
      title(main=paste(f$cell_id,'| Local expression ranks'),xlab=axis_label,cex.main=1.3)
      if (is.null(f$method)||f$method=='paired_oe_kd') mtext(c('Desired r (circle)','Negated opposite r (triangle)','Joint score (square)'),
            side=3,line=1,at=c(-.65,.1,.85),col=colors,cex=.8)
      else mtext('Desired correlation only. Opposite direction and joint score unavailable.',side=3,line=1,cex=.8)
    }
  }
}
render <- function(f) {
  if (streamed && !identical(f$schema, 'genosapiens.track2-native-figure-page/v1')) stop('Unsupported native figure contract')
  height<-if(f$kind=='heatmap')max(6,length(f$genes)*.3+3) else 7
  for(ext in c('png','svg')) {
    path<-file.path(out,paste0(f$id,'.',ext))
    if(file.exists(path))stop('Refuse to overwrite an existing figure')
    if(ext=='png')png(path,width=12,height=height,units='in',res=130,type='cairo') else svg(path,width=12,height=height)
    draw(f);dev.off()
  }
}
if (streamed) {
  source <- file(args[1], open='r')
  tryCatch({
    repeat {
      line <- readLines(source, n=1, warn=FALSE)
      if (!length(line)) break
      render(fromJSON(line, simplifyVector=FALSE))
    }
  }, finally=close(source))
} else {
  for (f in figures) render(f)
}
