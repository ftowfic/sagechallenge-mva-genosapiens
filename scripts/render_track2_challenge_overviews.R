# Versioned challenge graphics. Legacy renderer bytes remain separate.
args <- commandArgs(trailingOnly=TRUE)
records <- jsonlite::fromJSON(args[[1]], simplifyVector=FALSE)
out <- args[[2]]
wrap <- function(x, width=28) {
  lines <- strsplit(paste(x,collapse='\n'),'\n',fixed=TRUE)[[1]]
  wrapped <- lapply(lines,function(line) {
    tokens <- strsplit(line,' ',fixed=TRUE)[[1]]
    tokens <- vapply(tokens,function(token) {
      if(nchar(token)<=width) return(token)
      starts <- seq(1,nchar(token),by=width)
      paste(substring(token,starts,pmin(starts+width-1,nchar(token))),collapse=' ')
    },'')
    result <- strwrap(paste(tokens,collapse=' '),width=width)
    if(length(result)) result else ''
  })
  paste(unlist(wrapped),collapse='\n')
}
value <- function(x, missing='Unavailable') if(is.null(x)) missing else as.character(x)
method <- function(x) switch(x,paired_oe_kd='Paired OE/KD',kd_only_reversal='KD-only',
  oe_only_matching='OE-only',legacy_cell='Legacy cells',x)
state <- function(x) switch(x,positive='Candidates',negative='No candidates',outside_scope='Outside scope',
  incomplete='Unfinished',not_applicable='Not applicable',support='Support',conflict='Conflict',mixed='Mixed',
  no_two_direction_evidence='No paired direction',single_direction_nomination='Single-direction nomination',
  single_direction_no_nomination='No single-direction nomination',pass='Pass',fail='Fail',
  pending='Pending',unavailable='Unavailable',unassessed='Unassessed',failed='Failed',x)
fill <- function(x) if(x %in% c('positive','support','pass','single_direction_nomination')) '#D9EAF7' else
  if(x %in% c('conflict','fail','failed')) '#F8DDC4' else if(x %in% c('mixed','pending','incomplete')) '#FFF0BB' else '#ECEFF1'
route <- function(x) switch(x,lincs_cell_specific_expression='LINCS expression',
  qualified_target_drug_linkage='Reported mechanism',gsub('_',' ',x))
node <- function(values) {
  if(!length(values)) return('Unavailable\nNo retained assertion')
  paste(vapply(values,function(v) paste(value(v$label),value(v$origin),value(v$applicability),sep='\n'),''),collapse='\n\n')
}
tie_labels <- function(record) {
  labels <- character()
  for(metric in names(record$historical_ties)) for(tie in record$historical_ties[[metric]])
    if(tie$rank<=5) labels <- c(labels,paste('Historical',metric,':',tie$count,'tied at rank',tie$rank))
  for(group in record$derived_methods) for(tie in group$ties)
    if(tie$rank<=5) labels <- c(labels,paste(method(group$method),':',tie$count,'tied at rank',tie$rank))
  labels
}
render <- function(record) {
  n <- max(1,length(record$rows));width <- 13
  labels <- if(record$kind=='coverage') tie_labels(record) else character()
  height <- max(4,1.05*n+2.8+.25*length(labels))
  if(record$kind=='evidence') {
    columns <- length(unique(unlist(lapply(record$rows,function(r) vapply(r$states,function(x) x$label,'')))))
    width <- max(13,5+1.5*columns)
  }
  if(record$kind=='mechanism') {
    heights <- vapply(record$rows,function(r) {
      nodes <- lapply(list(r$variants,r$effects,r$functional,r$objectives),function(v) wrap(node(v),27))
      max(2.5,.16*max(vapply(strsplit(unlist(nodes),'\n'),length,0))+1)
    },0)
    height <- sum(heights)+1.5
  }
  for(format in c('png','svg')) {
    path <- file.path(out,paste0(record$id,'.',format))
    if(format=='png') png(path,width=width,height=height,units='in',res=120) else svg(path,width=width,height=height)
    par(family='sans',mar=c(3,1,3,1),xpd=NA)
    if(record$kind=='coverage') {
      routes <- unique(unlist(lapply(record$rows,function(r) vapply(r$routes,function(x) x$route,''))))
      m <- max(1,length(routes))
      par(mar=c(7+.7*length(labels),13,3,1),family='sans',xpd=NA)
      plot.new();plot.window(xlim=c(0,m),ylim=c(0,n))
      if(length(routes)) axis(1,at=seq_along(routes)-.5,labels=vapply(routes,function(x) wrap(route(x),16),''),
        tick=FALSE,las=2,cex.axis=.75)
      for(i in seq_along(record$rows)) {
        r <- record$rows[[i]];y <- n-i;text(-.08,y+.5,wrap(r$hypothesis,22),adj=1,cex=.8)
        for(j in seq_along(routes)) {
          found <- Filter(function(x) x$route==routes[[j]],r$routes)
          s <- if(length(found)) found[[1]]$category else 'unavailable'
          rect(j-1,y,j,y+1,col=fill(s),border='#53616B')
          text(j-.5,y+.5,wrap(state(s),15),cex=.75,font=2)
        }
      }
      if(!length(record$rows)) text(.5,.5,'No objectives reached a search.',cex=1)
      for(i in seq_along(labels)) mtext(labels[[i]],side=1,line=5.3+.7*(i-1),cex=.8,font=2,col='#23485C')
      title(paste('Figure',record$number,'. Search coverage and ties'),cex.main=1.05)
    } else if(record$kind=='mechanism') {
      plot.new();plot.window(xlim=c(0,4),ylim=c(0,sum(heights)))
      top <- sum(heights)
      for(i in seq_along(record$rows)) {
        r<-record$rows[[i]];bottom<-top-heights[[i]]
        text(.02,top-.12,paste(r$hypothesis,paste(unlist(r$gene_labels),collapse=', '),paste0('[S',i,']'),sep=' / '),adj=c(0,1),font=2,cex=.9)
        nodes<-list(r$variants,r$effects,r$functional,r$objectives)
        for(j in seq_len(4)) {
          rect(j-1+.055,bottom+.12,j-.055,top-.5,col=if(j==4) '#FFF0BB' else '#F0F5F7',border='#526C7A',lty=if(j==4) 2 else 1)
          text(j-.5,(bottom+top-.38)/2,wrap(node(nodes[[j]]),27),cex=.72)
          if(j<4) arrows(j-.05,(bottom+top)/2,j+.05,(bottom+top)/2,length=.06,col='#526C7A')
        }
        top<-bottom
      }
      title(paste('Figure',record$number,'. Evidence origin and proposed intervention'),cex.main=1.05)
      mtext('S labels identify source records in the figure data. Arrows do not establish rescue.',side=1,line=1,cex=.8)
    } else if(record$kind=='evidence') {
      par(mar=c(4,23,3,1),family='sans',xpd=NA)
      units<-unique(unlist(lapply(record$rows,function(r) vapply(r$states,function(x) x$label,''))))
      m<-max(1,length(units));plot.new();plot.window(xlim=c(0,m),ylim=c(0,n))
      if(length(units)) axis(1,at=seq_along(units)-.5,labels=vapply(units,wrap,'',width=14),tick=FALSE,cex.axis=.75)
      for(i in seq_along(record$rows)) {
        r<-record$rows[[i]];y<-n-i;e<-r$expression
        score<-if(is.null(e)) 'Original cell evidence' else paste0('Strength ',format(e$strength,digits=6),' | Percentile ',
          format(100*e$numerator/e$denominator,digits=6),'\nRank ',e$method_rank,' | Tied ',e$tie_count,
          '\nCalibration ',substr(e$entry_key,1,12))
        text(-.08,y+.5,paste(wrap(paste(r$name,r$gene,r$hypothesis,sep=' / '),36),method(r$method),score,sep='\n'),adj=1,cex=.72)
        for(j in seq_along(units)) {
          found<-Filter(function(x) x$label==units[[j]],r$states)
          s<-if(length(found)) found[[1]]$state else 'unavailable'
          mark<-state(s)
          if(length(found) && isTRUE(found[[1]]$conflict) && s!='conflict') mark<-paste(mark,'Conflict',sep='\n')
          if(length(found) && r$unit=='cell') mark<-paste(method(found[[1]]$method),mark,sep='\n')
          rect(j-1,y,j,y+1,col=fill(s),border='#53616B',lty=if(s=='pending') 2 else 1)
          if(s=='mixed') rect(j-1,y,j,y+1,density=10,angle=45,col='#B7A36B',border='#53616B')
          text(j-.5,y+.5,wrap(mark,17),cex=.72,font=2)
        }
      }
      title(paste('Figure',record$number,'. Method-specific expression evidence'),cex.main=1.05)
      mtext('Methods have separate ranks. Full identities and condition records remain in the figure data and supplements.',side=1,line=2.3,cex=.7)
    } else stop('Unsupported challenge overview kind.')
    dev.off()
  }
}
for(record in records) render(record)
jsonlite::write_json(list(R=R.version.string,platform=R.version$platform,jsonlite=as.character(packageVersion('jsonlite'))),
  file.path(out,'runtime.json'),auto_unbox=TRUE)
