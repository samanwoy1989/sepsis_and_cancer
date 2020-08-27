drawBoxplot3DiseaseGroup = function(pathid=sel, showTitle=T, pathgroup = "66") {
  xall = as.matrix(gsea.list$Zobs[pathid, ])
  xss = rowMeans(xall[,ids.sepsis])
  xca = rowMeans(xall[,ids.cancer])
  xslc = rowMeans(xall[,ids.sepsis.like])
  plotdat = list(xca, xslc, xss)
  titlestr = ""
  if(showTitle) {
    titlestr = paste0("Mean pathway score for the three disease groups\n(using ", pathgroup, " pathways)")
  }
  par(mar=c(4,10,4,5))
  b=boxplot(plotdat, names=c("CA","SLC","SS"), 
            main= titlestr, 
            #ylim=c(-200, 500),
            cex.lab=1.5, cex.axis=1.3, outline=F,
            boxwex=0.4,
            col=c("blue", "cyan", "purple"),
            main="",
            ylab="Pathway Score (Group-mean)")
  abline(h=0, lty=2, lwd=2, col="gray20")
  p_slc.ss = t.test(xslc-xss)$p.value
  p_ca.ss = t.test(xca-xss)$p.value
  p_slc.ca = t.test(xslc-xca)$p.value
  pvalues = c("p_slc.ss"=p_slc.ss,
             "p_ca.ss"=p_ca.ss,
             "p_slc.ca"=p_slc.ca)
  return(pvalues)
}
