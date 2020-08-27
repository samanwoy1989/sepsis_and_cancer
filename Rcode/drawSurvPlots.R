# Same barplot tried with ggplot
plot.dat <- data.frame(Cancers= names(percSignifSurv[ids.slc]), Percent_pathways= percSignifSurv[ids.slc])

bar.plot <- ggplot(data = plot.dat, aes(y=Percent_pathways, x= Cancers, fill= Cancers))+
  geom_bar(stat = "identity", width = 0.43, show.legend = F)+ theme_gray()+
  labs(title="Association of Selected Pathways \nwith Survival from Cancer in SLC group", 
  y = "Percentage of pathways \nassociated with survival")

rm(plot.dat)

##########################################
# Draw representative K-M plots
# Fig 4B
##########################################
pathid = "hsa04666"
proj = c("HNSC","KIRC","LIHC")[1]
#graphics.off(); rm(sfit, dat, info, pathX, projES, to.keep, sel.ptids)
graphics.off()
# get the clinical metadata from an external file and add to the expression set
info = read.table(file=paste0("Metadata/",proj,"_surv_clinical.tsv"), header=T, sep="\t")
to.keep = (apply(info[,2:3], 1, function(x) !all(x=="--")))
info = info[to.keep,]
is.na(info) <- info=="--"

rownames(info) = gsub("-",".", info$submitter_id)
sel.ptids = intersect(rownames(info), sampleNames(pathES))
#print.noquote(length(sel.ptids))
projES = pathES[, sel.ptids]
info = info[sel.ptids,]

pathX = as.vector(exprs(projES[pathid, ]))
names(pathX) = sampleNames(projES)
dat = data.frame(info, pathX)
dat$day = as.integer(as.character(dat$days_to_death))
dat$day[dat$outcome==0] = as.integer(as.character(dat$days_to_last_follow_up[dat$outcome==0]))
dat$high_expressing = as.integer(dat$pathX > median(dat$pathX))
dat = dat[!is.na(dat$day),]
#print.noquote(nrow(dat))
sfit <- survfit(Surv(day, outcome)~high_expressing, data=dat)
KM_plot_HNSC = ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE, 
                 legend.labs=c("Low", "High"), 
                 legend.title=paste0(pathid, " Score"),  
                 palette=c("dodgerblue2", "orchid2"), 
                 title=paste0("Kaplan-Meier Curve for survival in ", proj, "\nusing the pathway ", pathid), 
                 risk.table.height=.25)
pval = surv_pvalue(sfit, data = NULL, method = "survdiff",
                   test.for.trend = FALSE, combine = FALSE)$pval
rm(sfit, dat, info, pathX, projES, to.keep, sel.ptids)

#################################
#          Fig 4C
#################################
proj = c("HNSC","KIRC","LIHC")[2]
#graphics.off(); rm(sfit, dat, info, pathX, projES, to.keep, sel.ptids)
graphics.off()
# get the clinical metadata from an external file and add to the expression set
info = read.table(file=paste0("Metadata/",proj,"_surv_clinical.tsv"), header=T, sep="\t")
to.keep = (apply(info[,2:3], 1, function(x) !all(x=="--")))
info = info[to.keep,]
is.na(info) <- info=="--"

rownames(info) = gsub("-",".", info$submitter_id)

sel.ptids = intersect(rownames(info), sampleNames(pathES))
#print.noquote(length(sel.ptids))
projES = pathES[, sel.ptids]
info = info[sel.ptids,]

pathX = as.vector(exprs(projES[pathid, ]))
names(pathX) = sampleNames(projES)
dat = data.frame(info, pathX)
dat$day = as.integer(as.character(dat$days_to_death))
dat$day[dat$outcome==0] = as.integer(as.character(dat$days_to_last_follow_up[dat$outcome==0]))
dat$high_expressing = as.integer(dat$pathX > median(dat$pathX))
dat = dat[!is.na(dat$day),]
#print.noquote(nrow(dat))
sfit <- survfit(Surv(day, outcome)~high_expressing, data=dat)
KM_plot_KIRC = ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE, 
                 legend.labs=c("Low", "High"), 
                 legend.title=paste0(pathid, " Score"),  
                 palette=c("dodgerblue2", "orchid2"), 
                 title=paste0("Kaplan-Meier Curve for survival in ", proj, "\nusing the pathway ", pathid), 
                 risk.table.height=.25)
pval = surv_pvalue(sfit, data = NULL, method = "survdiff",
                   test.for.trend = FALSE, combine = FALSE)$pval
rm(sfit, dat, info, pathX, projES, to.keep, sel.ptids)

#################################
#          Fig 4D
#################################
proj = c("HNSC","KIRC","LIHC")[3]
#graphics.off(); rm(sfit, dat, info, pathX, projES, to.keep, sel.ptids)
graphics.off()
# get the clinical metadata from an external file and add to the expression set
info = read.table(file=paste0("Metadata/",proj,"_surv_clinical.tsv"), header=T, sep="\t")
to.keep = (apply(info[,2:3], 1, function(x) !all(x=="--")))
info = info[to.keep,]
is.na(info) <- info=="--"

rownames(info) = gsub("-",".", info$submitter_id)

sel.ptids = intersect(rownames(info), sampleNames(pathES))
#print.noquote(length(sel.ptids))
projES = pathES[, sel.ptids]
info = info[sel.ptids,]

pathX = as.vector(exprs(projES[pathid, ]))
names(pathX) = sampleNames(projES)
dat = data.frame(info, pathX)
dat$day = as.integer(as.character(dat$days_to_death))
dat$day[dat$outcome==0] = as.integer(as.character(dat$days_to_last_follow_up[dat$outcome==0]))
dat$high_expressing = as.integer(dat$pathX > median(dat$pathX))
dat = dat[!is.na(dat$day),]
#print.noquote(nrow(dat))
sfit <- survfit(Surv(day, outcome)~high_expressing, data=dat)
KM_plot_LIHC = ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE, 
                 legend.labs=c("Low", "High"), 
                 legend.title=paste0(pathid, " Score"),  
                 palette=c("dodgerblue2", "orchid2"), 
                 title=paste0("Kaplan-Meier Curve for survival in ", proj, "\nusing the pathway ", pathid), 
                 risk.table.height=.25)
pval = surv_pvalue(sfit, data = NULL, method = "survdiff",
                   test.for.trend = FALSE, combine = FALSE)$pval
rm(sfit, dat, info, pathX, projES, to.keep, sel.ptids)
fig = ggarrange(bar.plot, KM_plot_HNSC[[1]], KM_plot_KIRC[[1]], KM_plot_LIHC[[1]], 
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)
print(fig)
# Remove variables not needed
rm(fig, bar.plot, KM_plot_HNSC, KM_plot_KIRC, KM_plot_LIHC)