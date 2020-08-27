############################################################################################################
# Network analysis
############################################################################################################
# get a network of KEGG pathways (with size > 10)
which.large = sapply(genes.by.pathway, length)>10
ids = names(genes.by.pathway[which.large])

# make the pathway names shorter
shortids = sapply(strsplit(ids, split="hsa"), "[[", 2)
shortids = as.character(as.integer(shortids))
names(shortids) = ids

# net info
keggpairs = t(combn(shortids,2))
overlap = c()
for(i in 1:nrow(keggpairs)) {
  id1 = names(which(shortids==keggpairs[i,1]))
  id2 = names(which(shortids==keggpairs[i,2]))
  genes1 = genes.by.pathway[[id1]]
  genes2 = genes.by.pathway[[id2]]
  common = intersect(genes1, genes2)
  both = union(genes1, genes2)
  overlap[i] = length(common)/length(both)
}
netdat = data.frame(keggpairs, overlap)
colnames(netdat) = c("from","to","Overlap")

# keep those kegg pairs with overlap above a threshold (5%)
overlap.threshold = 0.05
links = netdat[netdat[,3]>overlap.threshold,]
nodes = unique(as.character(unlist(links[,1:2])))

# create a graph object
require(igraph)
g <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
# plot(g, vertex.label.color="black", vertex.label.dist=1, vertex.size=7)
E(g)$width = links[,3]

# Identify vertices with coresness <= threshold
toosmall <- which(graph.coreness(g) <= 2)
# Delete vertices and plot
gp <- delete.vertices(g, toosmall)

deg_dist_gp <- degree_distribution(gp, cumulative = TRUE)

#####################################################################
# Supplementary Figure S3_A
#####################################################################
plotdat <- data.frame(deg= 1:length(deg_dist_gp), cf= deg_dist_gp)
plot1=ggplot(plotdat, aes(y=cf , x=deg))+geom_point()+xlab("Degree")+
  ylab("Cumulative Frequency")+ 
  labs(title = "Degree distribution of \nKEGG pathway overlap network")

#####################################################################
# Supplementary Figure S3_B
#####################################################################
## Boxplot of degrees of selected and other pathways
dg = degree(gp)
names(dg) = paste0("hsa",formatC(as.integer(names(dg)), width=5, flag="0"))
selected = intersect(names(dg), sel)
dg.sel = dg[selected]
dg.other = dg[setdiff(names(dg),selected)]

Pathway.grp = rep("",length(dg))
names(Pathway.grp) = names(dg)
Pathway.grp[selected] = "Selected"
Pathway.grp[Pathway.grp!="Selected"] = "Other"

#plotdat<- rbind(cbind(rep("Other", length(dg.other)), dg.other), 
#                cbind(rep("Selected66", length(dg.sel)),dg.sel))

plotdat = data.frame("Pathways"=Pathway.grp, "val"=dg)

#plotdat<- as.data.frame(plotdat[,1:2], row.names = NULL)  
#colnames(plotdat) <- c("Pathways", "val")
plotdat$val <- as.numeric(plotdat$val)
plotdat= plotdat[!is.na(plotdat$val),]
plot2 = ggplot(plotdat, aes(x=Pathways, y=val, fill= Pathways))+
  geom_boxplot(show.legend = FALSE)+
  ylab("Degree of the nodes in \nthe KEGG overlap network")+ 
  labs(title = "Box plot showing Degree\nof the Selected Pathways")+ theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid"))
print(plot1)
print(plot2)
scan()

#####################################################################
# Figure 3
#####################################################################
# layout
coords = layout_(gp, nicely())

# color the nodes on the network
sel.nodes = c(pathways.c.only, pathways.sl.only, pathways.both)
nodes.66 = shortids[sel.nodes]

nodenames = vertex_attr(gp)$name
mycolors = rep("gray90",length(nodenames))
names(mycolors) = nodenames

mycolors[nodenames%in%nodes.66] = "red"
# resize the vertex size
node.size <- rep(4, length(V(gp)))

V(gp)$color = mycolors

#opar = par()
par(mar=c(1,1,1,1))
plot(gp, layout=coords, vertex.label.color=V(gp)$color, vertex.label.dist=1, vertex.label=NA,
     vertex.label.cex=0.8, vertex.label.dist=1, 
     vertex.size= node.size, vertex.frame.color=0, edge.curved=0.5)
cat("Selected pathways are at the core of the network!")
gc(reset=TRUE)
scan()

#####################################################################
# Betweenness calculation
#####################################################################

bet = betweenness(gp)
names(bet) = paste0("hsa",formatC(as.integer(names(bet)), width=5, flag="0"))
selected = intersect(names(bet), sel)
bet.sel = bet[selected]
bet.other = bet[setdiff(names(bet),selected)]

Pathway.grp = rep("",length(bet))
names(Pathway.grp) = names(bet)
Pathway.grp[selected] = "Selected"
Pathway.grp[Pathway.grp!="Selected"] = "Other"

#plotdat<- rbind(cbind(rep("Other", length(dg.other)), dg.other), 
#                cbind(rep("Selected66", length(dg.sel)),dg.sel))

plotdat = data.frame("Pathways"=Pathway.grp, "val"=bet)

#plotdat<- as.data.frame(plotdat[,1:2], row.names = NULL)  
#colnames(plotdat) <- c("Pathways", "val")
plotdat$val <- as.numeric(plotdat$val)
plotdat= plotdat[!is.na(plotdat$val),]
plot3 = ggplot(plotdat, aes(x=Pathways, y=val, fill= Pathways))+
  geom_boxplot(show.legend = FALSE)+
  ylab("Betweenness of the nodes in \nthe KEGG overlap network")+ 
  labs(title = "Box plot showing Centrality\nof the Selected Pathways")+
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid"))

print(plot3)

rm(plot1,plot2,plot3,bet,bet.other,bet.sel,coords,deg_dist_gp,dg,dg.other,
   dg.sel,g,i,keggpairs,mycolors,node.size,nodenames,nodes.66,
   Pathway.grp,sel.nodes,selected,shortids,which.large, ids,both,common,
   genes1,genes2,id1,id2,links,netdat,nodes,overlap,overlap.threshold,toosmall)
gc(reset=TRUE)
