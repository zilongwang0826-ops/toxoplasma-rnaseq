#!/usr/bin/env Rscript
# Presentation only: use frozen results; no count transformation or model fitting.
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(ggrepel)})
set.seed(645778)
ROOT <- normalizePath('.')
OUT <- 'report/figures/redesigned'
dir.create(file.path(OUT,'plot-data'),recursive=TRUE,showWarnings=FALSE)
INK <- '#243B49'; GREY <- '#B5BDC2'; BLUE <- '#286A8A'; AMBER <- '#C38317'; RUST <- '#B24D36'
STAGE <- c('11 dpi'=AMBER,'33 dpi'=RUST)
DIR <- c('up'=RUST,'down'=BLUE,'not_significant'=GREY)
theme_set(theme_classic(base_size=10, base_family='sans') + theme(
  text=element_text(colour=INK), plot.title=element_text(face='bold',size=11),
  plot.subtitle=element_text(size=9,margin=margin(b=9)), axis.title=element_text(size=10),
  axis.text=element_text(colour=INK,size=9), axis.line=element_line(colour='#71818A',linewidth=.3),
  axis.ticks=element_line(colour='#71818A',linewidth=.3), legend.position='bottom',
  legend.title=element_text(size=9), legend.text=element_text(size=8.5),
  plot.caption=element_text(size=8,colour='#54666E',hjust=0,margin=margin(t=10)),
  plot.margin=margin(9,12,8,8), strip.background=element_blank(),strip.text=element_text(face='bold',size=10)))
manifest <- list(); inputs <- character(); validation <- list()
read_tsv <- function(p) { inputs <<- unique(c(inputs,p)); read.delim(p,check.names=FALSE,stringsAsFactors=FALSE) }
write_tsv <- function(d,p) write.table(d,p,sep='\t',quote=TRUE,qmethod='double',row.names=FALSE,na='NA')
write_data <- function(d,stem) { p<-file.path(OUT,'plot-data',paste0(stem,'.tsv'));write_tsv(d,p);p }
export <- function(p,stem,w,h,source,purpose,data,category='supplementary') {
  # Main figures are reduced to the report's A4 text width. Keep labels readable
  # after that reduction; supplementary figures retain their standalone design.
  if(category=='main') p <- p & theme(
    axis.title=element_text(size=11), axis.text=element_text(size=11),
    axis.text.x=element_text(size=11), axis.text.y=element_text(size=11),
    legend.title=element_text(size=11), legend.text=element_text(size=11),
    plot.title=element_text(size=12), plot.subtitle=element_text(size=11),
    plot.caption=element_text(size=10), strip.text=element_text(size=11))
  for(ext in c('pdf')) {
    dest<-file.path(OUT,paste0(stem,'.',ext))
    if(ext=='pdf') pdf(dest,width=w,height=h,family='Helvetica',useDingbats=FALSE)
    if(ext=='svg') svg(dest,width=w,height=h,family='sans')
    if(ext=='png') png(dest,width=w,height=h,units='in',res=360,type='cairo',bg='white')
    print(p); dev.off()
  }
  manifest[[length(manifest)+1]] <<- data.frame(figure=stem,category=category,width_in=w,height_in=h,dpi=360,
    source_tables=paste(source,collapse=';'),plot_data=paste(data,collapse=';'),purpose=purpose,
    generator='scripts/figures/redesign_figures.R')
  message('Rendered ',stem)
}
panel <- function(p,letter) p+labs(tag=letter)+theme(plot.tag=element_text(face='bold',size=12),plot.tag.position=c(0,1))
params<-read_tsv('example-results/analysis/parameters.tsv')
alpha<-as.numeric(params$value[params$parameter=='alpha']); lfc<-as.numeric(params$value[params$parameter=='lfc'])
stopifnot(alpha==.05,lfc==1.5)
samples<-read_tsv('metadata/samples.tsv')
stopifnot(!anyDuplicated(samples$sample_id),nrow(samples)==12,all(table(samples$day,samples$condition)==3))
comp_all<-lapply(c(11,33),function(day) read_tsv(sprintf('example-results/comparison/day%s-published-transcript-comparison.tsv',day)))
symbols<-unique(do.call(rbind,lapply(comp_all,function(d)d[,c('gene_id','gene_symbol')])))
label_names<-c('Stat1','Irf1','Cxcl9','Cxcl10','H2-Ab1')
labels<-symbols[symbols$gene_symbol %in% label_names,]
stopifnot(!anyDuplicated(labels$gene_id),setequal(labels$gene_symbol,label_names))

place_pca_labels <- function(d) {
  # Deterministic finite candidate search; does not alter plotted coordinates.
  span<-diff(range(d$PC1)); if(span<=0)span<-1
  width<-span*.145; height<-span*.055
  xb<-range(d$PC1)+c(-.17,.2)*span
  yspan<-diff(range(d$PC2));yb<-range(d$PC2)+c(-.26,.26)*yspan
  angles<-seq(0,2*pi,length.out=17)[-17]
  offsets<-do.call(rbind,lapply(c(.08,.13,.19,.25),function(radius)cbind(cos(angles),sin(angles))*radius*span))
  placed<-matrix(NA_real_,nrow(d),2)
  for(i in seq_len(nrow(d))) {
    candidates<-sweep(offsets,2,c(d$PC1[i],d$PC2[i]),'+')
    costs<-apply(candidates,1,function(z){
      outside<-z[1]-width/2<xb[1]||z[1]+width/2>xb[2]||z[2]-height/2<yb[1]||z[2]+height/2>yb[2]
      prior<-if(i>1)seq_len(i-1) else integer()
      collide<-sum(abs(placed[prior,1]-z[1])<width*1.12 & abs(placed[prior,2]-z[2])<height*1.3)
      points<-sum(abs(d$PC1-z[1])<width*.85 & abs(d$PC2-z[2])<height*.85)
      outside*1e6+collide*1e4+points*1e3+sum((z-c(d$PC1[i],d$PC2[i]))^2)
    })
    placed[i,]<-candidates[which.min(costs),]
  }
  d$label_x<-placed[,1];d$label_y<-placed[,2];d
}

qc_panels <- function(level,scope) {
 base<-file.path('example-results/analysis',level,if(scope=='all')'qc' else paste0('day',scope,'/figures'))
 src<-file.path(base,c('pca-coordinates.tsv','pca-variance.tsv','sample-distances.tsv'))
 d<-read_tsv(src[1]);v<-read_tsv(src[2]);df<-read_tsv(src[3]); mat<-as.matrix(df[,-1]);rownames(mat)<-df$sample_id
 stopifnot(identical(d$sample_id,rownames(mat)),identical(colnames(mat),rownames(mat)),isTRUE(all.equal(mat,t(mat))),all(diag(mat)==0))
 idx<-match(d$sample_id,samples$sample_id);stopifnot(!anyNA(idx),identical(as.character(d$condition),as.character(samples$condition[idx])),identical(as.character(d$day),as.character(samples$day[idx])))
 d$stage<-factor(paste(d$day,'dpi'),levels=names(STAGE));d$condition<-factor(d$condition,levels=c('control','infected'));d<-place_pca_labels(d)
 # Equal coordinate units on both PCA axes prevent visual exaggeration.
 p<-ggplot(d,aes(PC1,PC2))+geom_hline(yintercept=0,linewidth=.25,colour='#E2E6E8')+geom_vline(xintercept=0,linewidth=.25,colour='#E2E6E8')+
 geom_point(aes(colour=stage,shape=condition),size=3,stroke=1)+
 geom_segment(aes(xend=label_x,yend=label_y),colour='#A3ADB1',linewidth=.25,show.legend=FALSE)+
 geom_text(aes(x=label_x,y=label_y,label=sample_id,colour=stage),size=if(level=='gene'&&scope=='all')4 else 3.3,show.legend=FALSE)+
 scale_colour_manual(values=STAGE,drop=TRUE,guide=guide_legend(order=1))+scale_shape_manual(values=c(control=1,infected=16),guide=guide_legend(order=2))+
 scale_x_continuous(expand=expansion(mult=c(.17,.2)))+scale_y_continuous(expand=expansion(mult=c(.26,.26)))+
 coord_fixed()+labs(title='Sample relationships',subtitle=sprintf('%s-level PCA | %s libraries',tools::toTitleCase(level),nrow(d)),x=sprintf('PC1 (%.1f%%)',100*v$fraction_variance[1]),y=sprintf('PC2 (%.1f%%)',100*v$fraction_variance[2]),colour=NULL,shape=NULL)
 ord<-hclust(as.dist(mat),method='complete')$order; ids<-rownames(mat)[ord]
 dd<-expand.grid(row=rownames(mat),column=colnames(mat));dd$distance<-as.vector(mat)
 dd$row<-factor(dd$row,levels=rev(ids));dd$column<-factor(dd$column,levels=ids)
 q<-ggplot(dd,aes(column,row,fill=distance))+geom_tile()+scale_fill_gradient(low='#F1F4F6',high='#244C69',name='VST distance',guide=guide_colourbar(display='rectangles'))+
 coord_equal()+labs(title='Expression distances',subtitle='Complete-linkage ordering',x=NULL,y=NULL)+
 theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=.5),legend.key.width=grid::unit(1.5,'cm'))
 stem<-if(level=='gene'&&scope=='all')'F01_sample_structure' else paste0('S_QC_',level,'_',scope)
 dat<-c(write_data(d,paste0(stem,'_pca')),write_data(dd,paste0(stem,'_distances')))
 composite<-panel(p,'A')+panel(q,'B')+plot_layout(widths=c(1.35,1))+plot_annotation(caption='Each point is one reported biological library. Open = control; filled = infected.\nPCA uses saved variance-stabilized coordinates; distance uses all count-filtered features.')
 export(composite,stem,10,4.9,src,'Sample identity and transformed-expression structure',dat,if(stem=='F01_sample_structure')'main' else 'supplementary')
 # Retain every feature in the saved 40-feature heatmap; reorder for readable group blocks.
 hs<-file.path(base,'heatmap-row-zscores.tsv');hm<-read_tsv(hs);stopifnot(!anyDuplicated(hm$feature_id),setequal(names(hm)[-1],d$sample_id))
 m<-as.matrix(hm[,-1]);rownames(m)<-hm$feature_id;rid<-rownames(m)[hclust(dist(m),method='complete')$order]
 sid<-samples$sample_id[samples$sample_id %in% colnames(m)]
 hd<-expand.grid(feature_id=rownames(m),sample_id=colnames(m));hd$z_score<-as.vector(m)
 hd$feature_id<-factor(hd$feature_id,levels=rev(rid));hd$sample_id<-factor(hd$sample_id,levels=sid)
 mm<-match(hd$sample_id,samples$sample_id);hd$group<-factor(paste(samples$day[mm],'dpi',samples$condition[mm]),levels=unique(paste(samples$day,'dpi',samples$condition)))
 lim<-max(abs(hd$z_score))
 hp<-ggplot(hd,aes(sample_id,feature_id,fill=z_score))+geom_tile()+
 facet_grid(~group,scales='free_x',space='free_x')+scale_fill_gradient2(low=BLUE,mid='white',high=RUST,limits=c(-lim,lim),name='Row Z score',guide=guide_colourbar(display='rectangles'))+
 labs(title=paste(tools::toTitleCase(level),'expression |',if(scope=='all')'all libraries' else paste(scope,'dpi')),subtitle='40 highest-variance features; saved row Z scores',x=NULL,y=NULL,caption='Columns follow study groups; complete-linkage row ordering. Feature identifiers are preserved.\nThe symmetric color scale is specific to this panel; row Z scores are not expression fold changes.')+
 theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.y=element_text(size=7.1),axis.text.x=element_text(angle=90,hjust=1,vjust=.5),panel.spacing.x=grid::unit(4,'pt'))
 hstem<-paste0('S_heatmap_',level,'_',scope);hdpath<-write_data(hd,hstem)
 export(hp,hstem,8.8,8.6,hs,'All saved high-variance features with explicit study groups',hdpath)
}
for(level in c('gene','transcript')) for(scope in c('all','11','33')) qc_panels(level,scope)

for(level in c('gene','transcript')) {
 sources<-sprintf('example-results/analysis/%s/day%s/results.tsv',level,c(11,33))
 allres<-do.call(rbind,lapply(seq_along(sources),function(i){d<-read_tsv(sources[i]);d$day<-c(11,33)[i];d}))
 stopifnot(!anyDuplicated(paste(allres$day,allres$feature_id)))
 selected<-!is.na(allres$padj)&allres$padj<alpha&abs(allres$log2FoldChange)>lfc
 stopifnot(identical(selected,allres$significant))
 d<-allres[is.finite(allres$log2FoldChange)&is.finite(allres$padj),]
 d$neglog10_fdr<--log10(pmax(d$padj,.Machine$double.xmin));d$stage<-factor(paste(d$day,'dpi'),levels=names(STAGE))
 d$gene_symbol<-labels$gene_symbol[match(d$feature_id,labels$gene_id)]
 summaries<-do.call(rbind,lapply(c(11,33),function(day) {s<-d[d$day==day,]; data.frame(day=day,stage=paste(day,'dpi'),label=sprintf('%s dpi  |  %s up / %s down\n%s FDR-eligible; %s missing adjusted P',day,sum(s$direction=='up'),sum(s$direction=='down'),nrow(s),sum(allres$day==day&!is.finite(allres$padj))))}))
 names_map<-setNames(summaries$label,summaries$stage)
 p<-ggplot(d,aes(log2FoldChange,neglog10_fdr))+geom_hline(yintercept=-log10(alpha),linetype='dashed',colour='#7B878D',linewidth=.3)+geom_vline(xintercept=c(-lfc,lfc),linetype='dashed',colour='#7B878D',linewidth=.3)+
 geom_point(aes(colour=direction),size=.65,alpha=.65)+scale_colour_manual(values=DIR,breaks=c('up','down','not_significant'),labels=c('Upregulated','Downregulated','Other FDR-eligible'),name=NULL)+
 facet_wrap(~stage,nrow=1,labeller=as_labeller(names_map))+labs(x='Unshrunken log2 fold change (infected / control)',y=expression(-log[10]*' BH-adjusted P'),caption='Primary selection: BH-adjusted P < 0.05 and |log2 fold change| > 1.5. All finite FDR-eligible estimates are shown.\nAdjusted P = 0 uses the machine-positive floor (2.225e-308); common axes preserve stage comparison.')+
 theme(strip.text=element_text(hjust=0,size=9.5),panel.spacing.x=grid::unit(20,'pt'))
 if(level=='gene') p<-p+geom_label_repel(data=d[!is.na(d$gene_symbol),],aes(label=gene_symbol),seed=645778,max.time=Inf,max.iter=10000,size=3.9,min.segment.length=0,box.padding=.5,max.overlaps=Inf,segment.colour=INK,fill='white',linewidth=0)
 stem<-if(level=='gene')'F02_gene_differential_expression' else 'S_transcript_differential_expression'
 export(p,stem,9.5,5.2,sources,'Paired time-specific differential-expression results',write_data(d,stem),if(level=='gene')'main' else 'supplementary')
 # MA coloring deliberately retains the original BH-only convention, separate from volcano selection.
 m<-allres[is.finite(allres$log2FoldChange)&is.finite(allres$baseMean)&allres$baseMean>0,]
 m$fdr_status<-ifelse(!is.finite(m$padj),'Adjusted P unavailable',ifelse(m$padj<alpha,'BH-adjusted P < 0.05','Other tested'))
 m$stage<-factor(paste(m$day,'dpi'),levels=names(STAGE))
 mp<-ggplot(m,aes(baseMean,log2FoldChange))+geom_hline(yintercept=0,colour='#7B878D',linewidth=.35)+geom_hline(yintercept=c(-lfc,lfc),linetype='dashed',colour='#A6AFB4',linewidth=.3)+
 geom_point(aes(colour=fdr_status),size=.6,alpha=.6)+scale_x_log10()+scale_colour_manual(values=c('BH-adjusted P < 0.05'=RUST,'Other tested'=GREY,'Adjusted P unavailable'='#677B84'),name=NULL)+facet_wrap(~stage,nrow=1)+
 labs(x='Mean normalized count (log10 scale)',y='Unshrunken log2 fold change',caption='Color follows the original MA convention: adjusted P only, without an effect-size filter.\nDashed lines mark +/-1.5. All finite effects are shown without the original +/-8 clipping.')
 stem<-paste0('S_MA_',level)
 export(mp,stem,9.5,4.8,sources,'Abundance versus effect with all finite estimates',write_data(m,stem))
}

# Six report-relevant themes, selected by biological role, not by new testing.
themes<-c('GO:0035456','GO:0034341','GO:0019882','GO:0050863','GO:0046651','GO:0002250')
go_sources<-sprintf('example-results/enrichment/day%s/GO/up.tsv',c(11,33))
go<-do.call(rbind,lapply(seq_along(go_sources),function(i){d<-read_tsv(go_sources[i]);d$day<-c(11,33)[i];d}))
sel<-go[go$term_id %in% themes,];stopifnot(nrow(sel)==12)
sel$stage<-factor(paste(sel$day,'dpi'),levels=names(STAGE));sel$key<-paste(sel$term_id,vapply(sel$term_name,function(x)paste(strwrap(x,width=39),collapse='\n'),character(1)),sep='\n')
keys<-sel$key[match(themes,sel$term_id)];sel$key<-factor(sel$key,levels=rev(keys))
sel$label<-sprintf('%s genes; FDR %s',sel$overlap,formatC(sel$padj,format='e',digits=1))
gp<-ggplot(sel,aes(enrichment_ratio,key,colour=stage))+geom_vline(xintercept=1,linetype='dashed',colour='#9EABB1',linewidth=.3)+
 geom_segment(aes(x=1,xend=enrichment_ratio,yend=key),linewidth=.5)+geom_point(aes(size=overlap))+geom_text(aes(x=26,label=label),hjust=0,size=3.8,colour=INK)+
 scale_colour_manual(values=STAGE,guide='none')+scale_size_area(max_size=6,name='Overlapping genes',breaks=c(50,100,200))+
 scale_x_continuous(limits=c(0,44),breaks=c(1,5,10,15,20,25))+facet_wrap(~stage,ncol=1)+
 labs(x='Enrichment ratio (observed / expected selected-gene proportion)',y=NULL,caption='Six biological themes chosen to represent the report narrative; they are not independent pathways.\nFDR uses the full size-eligible GO family per day and direction.\nAcross-day FDR differences are not interaction tests.')+
 theme(strip.text=element_text(hjust=0),axis.text.y=element_text(size=9.5),panel.spacing.y=grid::unit(15,'pt'))
export(gp,'F03_immune_programs',9.6,6.6,go_sources,'Matched biological themes with effect, overlap and adjusted evidence',write_data(sel,'F03_immune_programs'),'main')

for(kind in c('primary','collapsed')) for(day in c(11,33)) {
 src<-if(kind=='primary')sprintf('example-results/enrichment/day%s/GO/up.tsv',day) else sprintf('example-results/sensitivity/lilrb4a/day%s/up/collapsed.tsv',day)
 d<-read_tsv(src); d<-d[!is.na(d$padj)&d$padj<alpha,];d<-head(d[order(d$padj,d$pvalue,d$term_id),],15)
 d$key<-paste(d$term_id,vapply(d$term_name,function(x)paste(strwrap(x,width=50),collapse='\n'),character(1)),sep='\n');d$key<-factor(d$key,levels=rev(d$key));d$neglog10_fdr<--log10(pmax(d$padj,.Machine$double.xmin))
 p<-ggplot(d,aes(enrichment_ratio,key))+geom_vline(xintercept=1,colour='#ACB6BB',linetype='dashed',linewidth=.3)+geom_segment(aes(x=1,xend=enrichment_ratio,yend=key),colour='#CCD3D7',linewidth=.5)+geom_point(aes(size=overlap,colour=neglog10_fdr))+
 scale_size_area(max_size=7,name='Overlap',breaks=c(50,100,200),guide=guide_legend(order=1))+scale_colour_gradient(low='#E7CFA1',high=RUST,name='-log10 FDR',guide=guide_colourbar(display='rectangles',order=2))+scale_x_continuous(limits=c(0,NA))+
 labs(title=sprintf('%s dpi | %s upward GO terms',day,if(kind=='primary')'primary' else 'collapsed-locus'),subtitle='Top 15 by adjusted P, then raw P and GO identifier',x='Enrichment ratio (observed / expected)',y=NULL,caption='BH-adjusted P < 0.05; full term names and identifiers. Related parent/child terms are retained.\nPoint area is proportional to the number of selected overlapping genes.')+
 theme(axis.text.y=element_text(size=8),legend.position='bottom')
 stem<-sprintf('S_GO_%s_%s',kind,day);export(p,stem,10.5,9,src,'Complete top-15 term labels and quantitative legends',write_data(d,stem))
}
# Preserve all four nonsignificant downward results in one explicit status display.
downsrc<-c();down<-list()
for(kind in c('primary','collapsed')) for(day in c(11,33)) {
 src<-if(kind=='primary')sprintf('example-results/enrichment/day%s/GO/down.tsv',day) else sprintf('example-results/sensitivity/lilrb4a/day%s/down/collapsed.tsv',day)
 d<-read_tsv(src);downsrc<-c(downsrc,src)
 down[[length(down)+1]]<-data.frame(day=day,counting=kind,terms=nrow(d),significant=sum(d$padj<alpha,na.rm=TRUE),minimum_fdr=min(d$padj,na.rm=TRUE),selected_genes=unique(d$selected_size),universe=unique(d$universe_size))
}
down<-do.call(rbind,down);stopifnot(all(down$significant==0));down$row<-paste(down$day,'dpi |',down$counting)
down$detail<-sprintf('0 / %s significant GO terms\n%s selected genes / %s universe genes; minimum FDR = %.3f',down$terms,down$selected_genes,down$universe,down$minimum_fdr)
dp<-ggplot(down,aes(0,factor(row,levels=rev(row))))+geom_point(size=4,shape=1,colour=BLUE)+geom_text(aes(x=.08,label=detail),hjust=0,size=3.7,lineheight=1.2)+
 scale_x_continuous(limits=c(-.02,1))+labs(title='Downward gene lists: no significant GO terms',subtitle='Primary and overlapping-locus counting | 11 and 33 dpi',x=NULL,y=NULL,caption='All eligible terms remain in the BH family, including zero-overlap terms. Threshold: adjusted P < 0.05.\nA nonsignificant result does not establish absence of a biological response.')+
 theme(axis.text.x=element_blank(),axis.ticks=element_blank(),axis.line=element_blank(),axis.text.y=element_text(size=10))
export(dp,'S_GO_downward_status',9,4.3,downsrc,'Explicitly retain all nonsignificant downward enrichment results',write_data(down,'S_GO_downward_status'))

summarysrc<-'example-results/comparison/summary.tsv';cs<-read_tsv(summarysrc);cs<-cs[cs$selection=='significant',]
plots<-list();pairs<-list()
for(i in 1:2) {
 d<-comp_all[[i]];day<-c(11,33)[i];s<-cs[cs$day==day,]
 d<-d[is.finite(d$log2FoldChange_published)&is.finite(d$log2FoldChange_reanalysis),];stopifnot(nrow(d)==s$finite_effect_pairs)
 stopifnot(abs(cor(d$log2FoldChange_published,d$log2FoldChange_reanalysis,method='spearman')-s$effect_spearman)<1e-12)
 d$stage<-paste(day,'dpi');pairs[[i]]<-d
 p<-ggplot(d,aes(log2FoldChange_published,log2FoldChange_reanalysis))+geom_hline(yintercept=0,colour='#DCE1E4',linewidth=.3)+geom_vline(xintercept=0,colour='#DCE1E4',linewidth=.3)+geom_abline(slope=1,intercept=0,linetype='dashed',colour='#596D79',linewidth=.4)+geom_point(colour=STAGE[paste(day,'dpi')],alpha=.45,size=.8)+coord_fixed(xlim=c(-14,14),ylim=c(-14,14))+
 labs(title=sprintf('%s dpi | %s finite pairs',day,nrow(d)),subtitle=sprintf('Spearman rho = %.4f | same direction %.2f%%',s$effect_spearman,100*s$same_direction_fraction),x='Published log2 fold change',y='Reanalysis log2 fold change')
 plots[[i]]<-panel(p,LETTERS[i])
}
cs$stage<-factor(paste(cs$day,'dpi'),levels=rev(names(STAGE)))
cs$recovery_label<-sprintf('%s / %s selected transcripts recovered (%.2f%%)',cs$shared_selected,cs$published_selected,100*cs$published_list_recovery)
rp<-ggplot(cs,aes(stage,published_list_recovery,fill=stage))+geom_col(width=.4)+geom_text(aes(y=.02,label=recovery_label),hjust=0,colour='white',size=3.8)+coord_flip()+scale_y_continuous(limits=c(0,1),labels=function(x)paste0(x*100,'%'),breaks=c(0,.5,1))+scale_fill_manual(values=STAGE,guide='none')+
 labs(title='Published-list recovery',x=NULL,y=NULL)+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank())
cp<-(plots[[1]]+plots[[2]])/panel(rp,'C')+plot_layout(heights=c(4,1.35))+plot_annotation(caption='Correlation and direction agreement are conditional on finite effects in the published selected transcript list.\nThe dashed line denotes equal effects. Recovery uses all published selected transcripts;\nit is not genome-wide accuracy.')
export(cp,'F04_published_agreement',9.6,6.8,c(sprintf('example-results/comparison/day%s-published-transcript-comparison.tsv',c(11,33)),summarysrc),'Conditional effect agreement and selected-list recovery',c(write_data(do.call(rbind,pairs),'F04_published_agreement_pairs'),write_data(cs,'F04_published_agreement_recovery')),'main')
write_tsv(do.call(rbind,manifest),file.path(OUT,'figure-manifest.tsv'))
write_tsv(data.frame(path=inputs,md5=unname(tools::md5sum(inputs))),file.path(OUT,'input-checksums.tsv'))
writeLines(capture.output(sessionInfo()),file.path(OUT,'session-info.txt'))
