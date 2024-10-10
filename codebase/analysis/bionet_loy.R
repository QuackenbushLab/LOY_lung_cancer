# run BioNet
library(BioNet)
library(dplyr)
library(igraph)
library(stringr)

load("../data/kate_activesubnet_info.rdata")

top_left = matrix(0,ncol=735,nrow=735)
bottom_right = matrix(0,nrow=16856,ncol=16856)
colnames(network) = paste0(colnames(network),"-tf")
row.names(network) = paste0(row.names(network),"-gene")
adj_mat = cbind(rbind(top_left,network),rbind(t(network),bottom_right))

adj_graph = graph_from_adjacency_matrix(adj_mat,
                                        weighted=T,
                                        mode="undirected",
                                        diag=F)

names(node_info_tfs_pval) = paste0(names(node_info_tfs_pval),"-tf")
names(node_info_genes_pval) = paste0(names(node_info_genes_pval),"-gene")

names(node_info_tfs) = paste0(names(node_info_tfs),"-tf")
names(node_info_genes) = paste0(names(node_info_genes),"-gene")

node_p_vals = c(node_info_tfs_pval,node_info_genes_pval)
node_reg_scores = c(node_info_tfs,node_info_genes)

set.seed(1989)
bum = fitBumModel(node_p_vals)
node_scores = BioNet::scoreNodes(adj_graph, bum, fdr=5e-5)
bionet_mod = BioNet::runFastHeinz(adj_graph,node_scores)
length(V(bionet_mod))
sum(node_scores > 0)
# 337

# take the bionet_mod graph as an edgelist
bionet_edgelist = as_edgelist(bionet_mod) %>% 
  data.frame()

# run a basic over-representation analysis on the set of genes included
active_genes = unique(str_remove(bionet_edgelist$X2,"-gene"))
universe_genes = str_remove(names(node_info_genes),"-gene")

library(fgsea)
# gene set from https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2023.2.Hs/c2.cp.v2023.2.Hs.symbols.gmt
cp_gene_sets = gmtPathways("../ext/c2.cp.v2023.2.Hs.symbols.gmt")
reactome_gene_sets = cp_gene_sets[grep("REACTOME",names(cp_gene_sets))]
active_gene_ora = fora(pathways = cp_gene_sets,
                       genes = active_genes,
                       universe = universe_genes)
  
# reimpose the GRN edge weights
bionet_edgelist$reg_edge_weight=NA
for(i in 1:nrow(bionet_edgelist))
{
  gene_idx = which(rownames(network) == bionet_edgelist[i,]$X2)
  tf_idx = which(colnames(network) == bionet_edgelist[i,]$X1)
  bionet_edgelist[i,]$reg_edge_weight = network[gene_idx,tf_idx] 
}

# add the regulatory network edges back in
bionet_weighted_graph = graph_from_edgelist(as.matrix(bionet_edgelist[,1:2]))
E(bionet_weighted_graph)$reg_edge_weight = bionet_edgelist$reg_edge_weight
E(bionet_weighted_graph)$pos_reg_edge_weight = ifelse(bionet_edgelist$reg_edge_weight > 0,
                                                      bionet_edgelist$reg_edge_weight, 0)
E(bionet_weighted_graph)$neg_reg_edge_weight = ifelse(bionet_edgelist$reg_edge_weight < 0,
                                                      bionet_edgelist$reg_edge_weight, 0)

# calculate centrality measures (maybe not appropriate for bipartite graph)
bionet_measures = data.frame("strength" = strength(bionet_weighted_graph,
                                                   weights =abs(E(bionet_weighted_graph)$reg_edge_weight)),
                             "positive_strength" = strength(bionet_weighted_graph,
                                                            weights = E(bionet_weighted_graph)$pos_reg_edge_weight),
                             "negative_strength" = strength(bionet_weighted_graph,
                                                            weights=abs(E(bionet_weighted_graph)$neg_reg_edge_weight)),
                             "hub_score" = hub_score(bionet_weighted_graph, 
                                                     weights = abs(E(bionet_weighted_graph)$reg_edge_weight))$vector,
                             "betweenness" = betweenness(bionet_weighted_graph, weights = 1/abs(E(bionet_weighted_graph)$reg_edge_weight)))
                             
bionet_measures$node = row.names(bionet_measures)

# add node attributes to convert to bipartite network
V(bionet_weighted_graph)$type = FALSE
V(bionet_weighted_graph)$type[grep("-tf",V(bionet_weighted_graph)$name)] = TRUE

score_df = data.frame("node"=V(bionet_weighted_graph)$name) %>%
 left_join(data.frame("node"=names(node_p_vals),
                      "p_loy_assoc" = node_p_vals,
                      "FDR_loy_assoc"=p.adjust(node_p_vals,method="BH")),by="node") %>%
 left_join(data.frame("node"=names(node_scores),
                      "score_bum_loy"=node_scores),by="node") %>%
 left_join(data.frame("node"=names(node_reg_scores),
                      "node_loy_assoc"=node_reg_scores),by="node")

sum(V(bionet_weighted_graph)$name != score_df$node)
V(bionet_weighted_graph)$score_bum_loy = score_df$score_bum_loy
V(bionet_weighted_graph)$node_loy_assoc = score_df$node_loy_assoc

tfs_graph_idx = grep("-tf",V(bionet_weighted_graph)$name)

jpeg("loy_reg_net.jpeg",width=8,height=10,units="in",res=300)

#par(mar=c(0,0,0,0)+.1)
LO = layout_as_bipartite(bionet_weighted_graph,
                         hgap=10,vgap=10)
LO[tfs_graph_idx,1]=seq(from=3500,to=500,length.out=10)
LO[-tfs_graph_idx,1]=seq(from=4000,to=0,length.out=200)
#LO[-tfs_graph_idx,2] = LO[-tfs_graph_idx,2] + rnorm(n=length(tfs_graph_idx))
plot(bionet_weighted_graph,
     asp=2,
     layout=LO[,2:1],
     vertex.size = ifelse(V(bionet_weighted_graph)$type==T,25,4), 
     edge.arrow.size = 0,
     edge.color = ifelse(E(bionet_weighted_graph)$reg_edge_weight > 0, "red","blue"),
     edge.width = 2*abs(E(bionet_weighted_graph)$reg_edge_weight),
     vertex.label.degree = pi,
     vertex.label.dist = 7,
     #vertex.label.size=3,
     vertex.label.color="black",
     vertex.label = ifelse(V(bionet_weighted_graph)$type == "FALSE",NA,V(bionet_weighted_graph)$name),
     vertex.color = ifelse(V(bionet_weighted_graph)$node_loy_assoc > 0, "red4","royalblue"))

dev.off()

all_tf_nbhds = list()
for(i in 1:length(tfs_graph_idx))
{
  this_graph = make_ego_graph(bionet_weighted_graph, 1,V(bionet_weighted_graph)$name[tfs_graph_idx[i]])
  all_tf_nbhds[[i]] = this_graph
  names(all_tf_nbhds)[[i]] = V(bionet_weighted_graph)$name[tfs_graph_idx[i]]
}

circle_plot = function(mygraph,mytitle="My Graph")
{
  radian.rescale <- function(x, start=0, direction=1) {
    c.rotate <- function(x) (x + start) %% (2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  n = length(V(mygraph))
  V(mygraph)$type = F
  V(mygraph)$type[grep("-gene",V(mygraph)$name)] = T
  print(table(V(mygraph)$type))
  lab.locs <- radian.rescale(x=1:n, direction=-1, start=0)
  LO = layout_as_bipartite(mygraph,hgap=20,vgap=20)
  plot(mygraph,layout = LO[,2:1], #layout_as_bipartite,
       asp=2,
     #vertex.label.degree=lab.locs,
     vertex.size = ifelse(V(mygraph)$score_bum_loy>0,2*V(mygraph)$score_bum_loy,0), 
     edge.arrow.size = 0.5,
     edge.color = ifelse(E(mygraph)$reg_edge_weight > 0, "red","blue"),
     edge.width = 3*abs(E(mygraph)$reg_edge_weight),
     vertex.label.dist = ifelse(V(mygraph)$type==T,3,-2),
     vertex.label.color="black",
     vertex.label.cex=ifelse(V(mygraph)$type==T,0.5,1),
     vertex.label.degree=ifelse(V(mygraph)$type==T,pi,pi/2),
     #vertex.label.size=3,
     vertex.label = str_remove(V(mygraph)$name,"-gene"), # ifelse(V(mygraph)$type == "FALSE",NA,V(bionet_weighted_graph)$name),
     vertex.color = ifelse(V(mygraph)$node_loy_assoc > 0, "red4","royalblue"))
  title(mytitle)
  }

#ZHX2 and GATA3
pdf("tf_nbhds.pdf",width=5,height=10)
par(mfrow=c(1,1))
for(i in 1:length(all_tf_nbhds))
{
  print(i)
  circle_plot(all_tf_nbhds[[i]][[1]],mytitle=names(all_tf_nbhds)[[i]])
}
dev.off()

table(V(bionet_weighted_graph)$node_loy_assoc > 0)
# 201     9
sum(E(bionet_weighted_graph)$reg_edge_weight > 0)
#215
sum(E(bionet_weighted_graph)$reg_edge_weight < 0)
#18
length(E(bionet_weighted_graph))
#233

measures_and_scores = bionet_measures %>% left_join(score_df,by="node")
table_paper = measures_and_scores %>% filter(grepl("tf",node)) %>%
  select(node,strength,positive_strength, 
         negative_strength,node_loy_assoc,p_loy_assoc,FDR_loy_assoc) %>% 
  arrange(-strength)
table_paper[,2:5] = round(table_paper[,2:5],2)
table_paper$p_loy_assoc = formatC(table_paper$p_loy_assoc, format = "e", digits = 2)
table_paper$FDR_loy_assoc = formatC(table_paper$FDR_loy_assoc, format = "e", digits = 2)

names(table_paper)=c("Node","Absolute Weighted Degree","Positive Degree","Negative Degree","LOY Effect","LOY p","LOY FDR")
write.table(table_paper,file="table_paper.asv",
            sep="&",quote=F)

jpeg("cytochrome_p450_paths.jpeg",width=12,height=4,units="in",res=300)
par(mfrow=c(1,3))
# pull out subnets for some of the enriched pathways
kegg_drug_metabolism_cytochrome_p450 = c("ADH1B", 
                                         "ADH1C", 
                                         "ADH5", 
                                         "CYP2A13", 
                                         "FMO5", 
                                         "GSTA3", 
                                         "MAOA")

wp_cytochrome_p450 = c("CYP2A13", 
                       "CYP4A11", 
                       "CYP4A22", 
                       "CYP4B1", 
                       "CYP4V2", 
                       "CYP4X1")

reactome_cytochrome_p450 = c("ARNT2",
                             "CYP2A13", 
                             "CYP4A11", 
                             "CYP4A22", 
                             "CYP4B1", 
                             "CYP4V2")

this_graph = make_ego_graph(bionet_weighted_graph, 
                            order=1,
                            nodes = paste(kegg_drug_metabolism_cytochrome_p450,"-gene",sep=""))

whole_graph = do.call(union, this_graph)
V(whole_graph)$type = F
V(whole_graph)$type[grep("-tf",V(whole_graph)$name)] = T
my_layout = layout_as_bipartite(whole_graph)
plot(whole_graph,
     layout = my_layout[,2:1],
     vertex.color="white",
     vertex.size=4,vertex.label.dist=2,
     vertex.label.degree=-90)
title("KEGG Drug Metabolism Cytochrome p450")

this_graph = make_ego_graph(bionet_weighted_graph, 
                            order=1,
                            nodes = paste(wp_cytochrome_p450,"-gene",sep=""))

whole_graph = do.call(union, this_graph)
V(whole_graph)$type = F
V(whole_graph)$type[grep("-tf",V(whole_graph)$name)] = T
my_layout = layout_as_bipartite(whole_graph)
plot(whole_graph,
     layout = my_layout[,2:1],
     vertex.color="white",
     vertex.size=4,vertex.label.dist=2,
     vertex.label.degree=-90)
title("WikiPathways Oxidation by Cytochrome p450")


this_graph = make_ego_graph(bionet_weighted_graph, 
                            order=1,
                            nodes = paste(reactome_cytochrome_p450,"-gene",sep=""))

whole_graph = do.call(union, this_graph)
V(whole_graph)$type = F
V(whole_graph)$type[grep("-tf",V(whole_graph)$name)] = T
my_layout = layout_as_bipartite(whole_graph)
plot(whole_graph,
     layout = my_layout[,2:1],
     vertex.color="white",
     vertex.size=4,vertex.label.dist=2,
     vertex.label.degree=-90)
title("Reactome Cytochrome p450")
dev.off()
