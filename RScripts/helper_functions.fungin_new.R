
fungin_initialize <- function(fungin_data_dir=file.path("data","fungin"), sourcedb="reactomeFI",
                              rds_file=NULL, fromScratch=F, genome="hg38", string_score_threshold=NULL) {
  require(igraph)

  sourcedb=tolower(sourcedb)
  if (! sourcedb %in% c("reactomefi","stringdb")) {
    stop("'sourcedb' must be one of 'rectomeFI' or 'stringDB'.")
  }

  fungin_data_dir <- file.path(fungin_data_dir,sourcedb)
  if (is.null(rds_file)) {
    rds_file <- ifelse(sourcedb=="stringdb", file.path(fungin_data_dir,paste0("fungin.",sourcedb,".",genome,".Rds")),
                                               file.path(fungin_data_dir,paste0("fungin.",sourcedb,".Rds")))
  }
  if (!file.exists(rds_file) | fromScratch) {
    warning(paste0("Building network data from scratch (",sourcedb,")...\n"))
    if (sourcedb=="reactomefi") {
      rds_file <- gather_FI_network_data(fungin_data_dir,rds_file=rds_file, fromScratch=fromScratch)


    } else if (sourcedb=="stringdb") {
      rds_file <- gather_stringdb_data(fungin_data_dir,rds_file=rds_file, fromScratch=fromScratch, genome=genome)
    }
  }
  load(rds_file)
  if (sourcedb=="stringdb" & !is.null(string_score_threshold)) {
    # browser()
    # get.edge.attribute(full_interaction_network,"combined_score")
    # full_interaction_network <- delete.edges(full_interaction_network, which(E(full_interaction_network)$combined_score >= string_score_threshold)-1)
    # string_score_threshold=800
    # load(rds_file)
    # dim(get.edgelist(full_interaction_network))
    # range(edge_attr(full_interaction_network,"combined_score"))
    # # full_interaction_network <- full_interaction_network %>% delete_edges(which(edge_attr(full_interaction_network,"combined_score") >= string_score_threshold)-1)
    full_interaction_network <- delete.edges(full_interaction_network, E(full_interaction_network)[!edge_attr(full_interaction_network,"combined_score") >= string_score_threshold])
    # dim(get.edgelist(full_interaction_network))
    # range(edge_attr(full_interaction_network,"combined_score"))
    # # isolated_nodes <- which(strength(full_interaction_network, mode="all") < 1)
    isolated_nodes <- which(igraph::degree(full_interaction_network) < 1)
    nodes_to_remove <- isolated_nodes[!is.na(isolated_nodes)]
    full_interaction_network <- igraph::simplify(igraph::delete_vertices(full_interaction_network, nodes_to_remove), edge.attr.comb="first")
  }
  full_interaction_network <- fungin_add_attr(full_interaction_network)


  return(full_interaction_network)
}

fungin_add_attr <- function(igraph_obj) {
  all_node_attrs <- names(edge_attr(igraph_obj))
  all_edge_attrs <- names(edge_attr(igraph_obj))
  if (! "logFC" %in% all_node_attrs) {
    igraph_obj <- set_vertex_attr(igraph_obj,name="logFC",value=NA)
  }
  if (! "pval" %in% all_node_attrs) {
    igraph_obj <- set_vertex_attr(igraph_obj,name="pval",value=NA)
  }
  if (! "pval_binary" %in% all_node_attrs) {
    igraph_obj <- set_vertex_attr(igraph_obj,name="pval_binary",value=NA)
  }
  if (! "mutated" %in% all_node_attrs) {
    igraph_obj <- set_vertex_attr(igraph_obj,name="mutated",value=NA)
  }
  if (! "changetype" %in% all_edge_attrs) {
    igraph_obj <- set_edge_attr(igraph_obj,name="changetype",value="Other")
  }
  return(igraph_obj)
}

fungin_trim <- function(full_interaction_network, query_genes=NULL, get_neighbors=FALSE,get_largest_connected_subnet=FALSE) {
  require(igraph)
  query_genes <- query_genes[query_genes %in% names(V(full_interaction_network))]
  if (length(query_genes)==0) {
    stop("Query genes not found in interaction data")
    # stop(paste0("Top ", n_genes_cutoff, " not found in interaction data"))
    # return(NA)
  }

  if (get_neighbors) {
    neighborhood <- ego(full_interaction_network, nodes=V(full_interaction_network)[query_genes])
    query_nodes <- unlist(neighborhood)
  } else {
    query_nodes <- V(full_interaction_network)[query_genes]
  }

  # my_interaction_graph <- simplify(induced_subgraph(full_interaction_network,query_nodes))#, edge.attr.comb="sum")
  my_interaction_graph <- igraph::simplify(induced_subgraph(full_interaction_network,query_nodes), edge.attr.comb="first")

  full_strengths <- strength(my_interaction_graph, mode="all")
  my_interaction_graph <- set_vertex_attr(my_interaction_graph,"strength",V(my_interaction_graph)[names(full_strengths)],full_strengths)

  isolated_nodes <- V(my_interaction_graph)$name[V(my_interaction_graph)$strength < 1]

  if (length(isolated_nodes)/length(V(my_interaction_graph)) > 0.8) {
    stop(paste0("Not enough interactions left to plot.",ifelse(get_neighbors,""," Try setting 'get_neighbors' to TRUE.")))
    # warning(paste0("Not enough interactions to plot.",ifelse(get_neighbors,""," Try setting 'get_neighbors' to TRUE.")))
    # return(NA)
  }
  # nodes_to_remove <- union(nonsig_nodes, isolated_nodes)
  nodes_to_remove <- isolated_nodes[!is.na(isolated_nodes)]

  # browser()
  # curr_graph <- simplify(delete_vertices(curr_graph, nodes_to_remove), edge.attr.comb="sum")
  # my_interaction_graph <- simplify(delete_vertices(my_interaction_graph, nodes_to_remove))
  my_interaction_graph <- igraph::simplify(igraph::delete_vertices(my_interaction_graph, nodes_to_remove), edge.attr.comb="first")

  if(get_largest_connected_subnet)
  {
    g<-components(my_interaction_graph)
    largestCluster<-names(sort(table(g$membership),decreasing = T)[1])
    genesInCluster<-igraph::groups(g)[[largestCluster]]
    genesNotInSubnetwork<-setdiff(V(my_interaction_graph)$name,genesInCluster)
    my_interaction_graph <- igraph::simplify(igraph::delete_vertices(my_interaction_graph, genesNotInSubnetwork), edge.attr.comb="first")
  }

  return(my_interaction_graph)
}


fungin_annotate <- function(fungin_graph,
                            diff_exp_results=NULL, pval_cutoff=0.05, n_genes_cutoff=10,
                            gene_column="gene",fc_column="logFC",pval_column="adj.P.Val",expression_column="exp",isTFFactor_column="isTF",
                            isDrugTarget="isdrugTarget",maf=NULL, min_mutated_samples = 1,
                            fillVar="mutated") {
  require(igraph)
  # browser()
  all_vertex_attr <- data.frame(do.call(cbind,vertex_attr(fungin_graph)), stringsAsFactors = F)
  if (is.null(diff_exp_results)) {
    # diff_exp_results <- as.data.frame(matrix(0, nrow=1, ncol=3))
    # colnames(diff_exp_results) <- c(gene_column, fc_column, pval_column)
    diff_exp_results <- data.frame(name=names(V(fungin_graph)),
                                   all_vertex_attr[,c("logFC","pval")],
                                   stringsAsFactors = F)
  } else {
    if (any(! c(gene_column, fc_column, pval_column,expression_column) %in% colnames(diff_exp_results))) {
      stop(paste0("Diff exp results are missing these columns: ", paste0(setdiff(c(gene_column, fc_column, pval_column,expression_column), colnames(diff_exp_results)), collapse=",")))
    }
    diff_exp_results <- diff_exp_results[,c(gene_column, fc_column, pval_column,expression_column,isTFFactor_column,isDrugTarget)]
    colnames(diff_exp_results) <- c("name","logFC","pval","Expression","isTF","isdrugTarget")
  }

  diff_exp <- diff_exp_results[order(diff_exp_results$pval, decreasing = F),]

  mut_status <- data.frame(name=names(V(fungin_graph)),
                           mutated=all_vertex_attr[,c("mutated")],
                           stringsAsFactors = F)
  if (! is.null(maf)) {
    # browser()
    print("In MAF function")
    mut_summary <- make_mut_summary(maf)
    # mut_summary$mutated <- ifelse(mut_summary$AlteredSamples >= min_mutated_samples,
    #                                paste0("Mutated in ", min_mutated_samples, " or more samples"),
    #                                NA)
    mut_summary$mutated <- unlist(lapply(mut_summary$AlteredSamples,FUN=function(x){
      ifelse(x >= 10, paste0("Mutated in 10 or more samples"),
             ifelse(x >= 5, paste0("Mutated in 5 to 9 samples"),
                    ifelse(x > 1,paste0("Mutated in 2 to 4 samples"),
                           ifelse(x == 1,paste0("Mutated in 1 sample"), NA))))
    }))
    mut_status <- data.frame(name=mut_summary$Hugo_Symbol,
                             mutated=mut_summary$mutated,
                             Nmutated=mut_summary$AlteredSamples,
                             stringsAsFactors = F)
  }

  curr_graph <- fungin_graph
  my_attrs <- data.frame(vname=V(curr_graph)$name, stringsAsFactors = F)

  my_attrs$logFC <- diff_exp$logFC[match(my_attrs$vname, diff_exp$name)]
  my_attrs$exp <- diff_exp$Expression[match(my_attrs$vname, diff_exp$name)]
  my_attrs$isTF <- diff_exp$isTF[match(my_attrs$vname, diff_exp$name)]
  my_attrs$isDrugTarget <- diff_exp$isdrugTarget[match(my_attrs$vname, diff_exp$name)]
  #print(head(diff_exp$isdrugTarget[match(my_attrs$vname, diff_exp$name)]))

  my_attrs$pval <- diff_exp$pval[match(my_attrs$vname, diff_exp$name)]
  # my_attrs$pval_binary <- ifelse(my_attrs$pval < pval_cutoff, paste0("p < ",pval_cutoff),
  #                                ifelse(my_attrs$pval >= pval_cutoff, "ns",NA))
  my_attrs$pval_binary <- my_attrs$pval
  my_attrs$pval_binary[is.na(my_attrs$pval_binary)] <- "No data"

  my_attrs$mutated <- mut_status$mutated[match(my_attrs$vname, mut_status$name)]
  my_attrs$mutated[is.na(my_attrs$mutated)] <- "No data"
  my_attrs$Nmutated <- mut_status$Nmutated[match(my_attrs$vname, mut_status$name)]
  # print(table(my_attrs$mutated))
  ###Community Clustering
  # curr_graph.communities <- cluster_walktrap(curr_graph, weights=NULL)#edge.betweenness.community(curr_graph, weights=NULL, directed=FALSE)
  # crossingObj<-crossing(curr_graph.communities,curr_graph)
  # my_attrs$membership<-membership(curr_graph.communities)
  #curr_graph.clustering <- make_clusters(curr_graph, membership=curr_graph.communities$membership)
  # browser()
  #print(head(my_attrs))
  curr_graph <- set_vertex_attr(curr_graph,"logFC",V(curr_graph)[my_attrs$vname],my_attrs$logFC)
  curr_graph <- set_vertex_attr(curr_graph,"pval",V(curr_graph)[my_attrs$vname],my_attrs$pval)
  curr_graph <- set_vertex_attr(curr_graph,"Expression",V(curr_graph)[my_attrs$vname],my_attrs$exp)
  curr_graph <- set_vertex_attr(curr_graph,"pval_binary",V(curr_graph)[my_attrs$vname],my_attrs$pval_binary)
  curr_graph <- set_vertex_attr(curr_graph,"mutated",V(curr_graph)[my_attrs$vname],my_attrs$mutated)
  curr_graph <- set_vertex_attr(curr_graph,"Nmutated",V(curr_graph)[my_attrs$vname],my_attrs$Nmutated)
  curr_graph <- set_vertex_attr(curr_graph,"isTF",V(curr_graph)[my_attrs$vname],my_attrs$isTF)
  curr_graph <- set_vertex_attr(curr_graph,"isDrugTarget",V(curr_graph)[my_attrs$vname],my_attrs$isDrugTarget)
  #curr_graph <- set_vertex_attr(curr_graph,"membership",V(curr_graph)[my_attrs$vname],my_attrs$membership)
  # nonsig_nodes <- V(curr_graph)$name[V(curr_graph)$pval > pval_cutoff]
  #curr_graph <- set_edge_attr(curr_graph,"same.community",E(curr_graph),ifelse(crossingObj, 2, 1))
  return(curr_graph)
}


fungin_make_plotdata <- function(fungin_graph) {
  require(igraph)
  L <- graph_attr(fungin_graph %>%
                    add_layout_(nicely(), component_wise()),
                  "layout")


  vis.nodes <- igraph::as_data_frame(fungin_graph,what = "vertices")
  vis.nodes$label  <- vis.nodes$name
  vis.nodes$Xn <- L[,1]
  vis.nodes$Yn <- L[,2]
  vis.nodes$color <- "lightblue"


  vis.nodes$shape <- ifelse(vis.nodes$mutated=="No data","circle",
                            ifelse(grepl("Mutated",vis.nodes$mutated), "star","square"))

  vis.nodes$color.background <- "grey90"
  # browser()
  maxFC<-NA
  if (!all(is.na(vis.nodes$logFC))) {
    # suppressWarnings(maxFC<-max(abs(vis.nodes$logFC), na.rm = T))
    maxFC<-max(abs(vis.nodes$logFC), na.rm = T)
  }else if(!all(is.na(vis.nodes$Expression))){
    maxFC<-max(abs(vis.nodes$Expression), na.rm = T)
  }
  if (!is.na(maxFC)) {
    require(circlize)
    fc_color_breaks <- seq(-maxFC, maxFC, length.out=6)
    # maxFC=4
    fc_color_breaks <- seq(-maxFC, maxFC, length.out=6)
    color_func <- colorRamp2(fc_color_breaks, brewer.pal(name = "PiYG",n = 6))
    vis.nodes$color.background <- color_func(vis.nodes$logFC)
  }
  vis.nodes$color.border <- "grey50"
  vis.nodes$color.border[vis.nodes$id %in% query_genes] <- "red"
  # browser()

  vis.links <- igraph::as_data_frame(fungin_graph,what = "edges")
  # vis.links <- merge.data.frame(vis.links, gene_interactions, by.x=c("from","to"), by.y=c("Gene1","Gene2"), all.x=T)
  # vis.links$changetype <- ifelse(vis.links$change == 0, "Other",ifelse(vis.links$change < 0, "Inhibition","Activation"))
  # vis.links$width <- 4
  # vis.links$arrows <- "middle"

  interaction_colors <- c("Other"="grey70", "Inhibition"="blue2","Activation"="gold")
  # vis.links$color.background <- interaction_colors[vis.links$changetype]
  vis.links$color <- interaction_colors[vis.links$changetype]

  vis.links$width <- ifelse(vis.links$changetype == "Other", 1,3)
  vis.links$arrowType <- ifelse(vis.links$changetype == "Other",0,ifelse(vis.links$changetype=="Inhibits",7,2))
  # browser()
  # marker_offset <- marker_size/100
  vis.links$startx <- vis.nodes$Xn[match(vis.links$from, vis.nodes$name)]#+marker_offset
  vis.links$starty <- vis.nodes$Yn[match(vis.links$from, vis.nodes$name)]#+marker_offset
  vis.links$endx <- vis.nodes$Xn[match(vis.links$to, vis.nodes$name)]#-marker_offset
  vis.links$endy <- vis.nodes$Yn[match(vis.links$to, vis.nodes$name)]#-marker_offset


  return(list(nodes=vis.nodes, edges=vis.links))
}



fungin_plot <- function(fungin_graph, marker_size=40, fillVar="logFC",savename=NULL) {


  # browser()
  require(ggnetwork)
  require(ggnewscale)
  layout_area_param=NULL
  plotdata <- ggnetwork(fungin_graph,scale = F, stringsAsFactors=F)
  edgesize <- 1.5
  shape_vals <- c(21, 22, 23, 24, 25)
  names(shape_vals)<-c("Mutated in 1 sample","Mutated in 2 to 4 samples","Mutated in 5 to 9 samples","Mutated in 10 or more samples","No data")
  shape_vals<-shape_vals[unique(plotdata$mutated)]
  #names(shape_vals) <- sort(unique(plotdata$mutated))
  #print(unique(plotdata$mutated))

  interaction_colors <- c("Other"="grey70", "Inhibition"="blue2","Activation"="gold")
  #plotdata$pval_binary[is.na(plotdata$pval_binary)] <- "NA"
  alpha_vals <- rep(1,length(unique(plotdata$pval_binary)))#c(1, 0.25, 0.1)
  names(alpha_vals) <- rev(sort(unique(plotdata$pval_binary)))

  outline_vals <- c("blue2","grey", "green","orange","yellow")
  names(outline_vals) <- c("SNVs","No data","CNVs-DEL","CNVs-AMP","CNVs-SNVs")
  outline_vals <- outline_vals[names(alpha_vals)]
  #names(outline_vals) <- names(alpha_vals)

  nedges <- length(E(fungin_graph))
  edgesize <- ifelse(nedges > 1000, 0.25,
                     ifelse(nedges > 100, 0.5,
                            ifelse(nedges > 10, 2, edgesize)))

  nnodes <- length(V(fungin_graph))
  nodesize <- ifelse(nnodes > 1000, 2,
                     ifelse(nnodes > 100, 5,
                            ifelse(nnodes > 10, 10, nnodes)))
  # noOfMutations<- plotdata$Nmutated
  # nodesize <- unlist(lapply(noOfMutations,FUN=function(x){
  #   ifelse(x >= 10, 3,ifelse(x >= 5, 3,ifelse(x > 2, 2, ifelse(x == 1, 1, x))))
  # }))

  if (!fillVar %in% colnames(plotdata)) {
    warning("fillVar not found, using logFC for fill colors...")
    fillVar="logFC"
  }
  plotdata$fillvar <- plotdata[,fillVar]
  mypal <- "PiYG"
  maxFillValue<-max(abs(plotdata$fillvar), na.rm = T)
  color_limits <- c(-maxFillValue, maxFillValue)
  showFClegend=T
  if (all(is.na(plotdata$fillvar))) {
    plotdata$logFC <- 0
    showFClegend=F
    # suppressWarnings(maxFC<-max(abs(vis.nodes$logFC), na.rm = T))
  } else {
    if (min(plotdata$fillvar, na.rm=T) >= 0) {
      mypal <- "Reds"
      # color_limits <- range(plotdata$fillvar, na.rm = T)
      color_limits <- c(0, max(plotdata$fillvar, na.rm = T))
    }
  }


  # gene_name_colors <- c("#ff3300","#cc0000","grey50")
  # gene_name_colors <- c("#ff3300","#cc0000","grey50")
  require(randomcoloR)
  gene_name_colors <- distinctColorPalette(length(unique(plotdata$mutated)))
  names(gene_name_colors) <- sort(unique(plotdata$mutated))
  gene_name_colors["No data"] <- "black"
  print(length(((plotdata$mutated))))
  #print(head(plotdata$isTF))
  plotdata$isTFColor<-unlist(lapply(plotdata$isTF,FUN=function(x){if (x) "Yes" else "No"}))
  labelNamesColors<-c("red","black")
  names(labelNamesColors)<-c("Yes","No")
  # browser()
  ##Adding drug bank information
  plotdata$name[plotdata$isDrugTarget == "Y"]<-paste0(plotdata$name[plotdata$isDrugTarget == "Y"],"*")
  require("ggforce")
  require("concaveman")
  edgeplot <- ggplot(plotdata,
                   # layout = "fruchtermanreingold", cell.jitter = 2, niter=1000, area=layout_area_param,
                   # layout = "mds",
                   # layout = "spring", repulse=T, mass=0.1, k=0.001, kfr=0.01, repeqdis=0.5,
                   # layout = "kamadakawai", niter=10000,initemp=1000,cool.exp=0.1,
                   # layout = "mds", niter=10, var="geodist", dist="maximum",
                   # layout = "hall", niter=100,
                   aes(x = x, y = y, xend = xend, yend = yend))+#,
    # layout = "target", niter=1000) +
    # geom_edges(color= "grey50", alpha=0.5,size=0.1,curvature = 0.2,
    #ggalt::geom_encircle(expand=0,aes(color=as.factor(membership))) +
    # geom_edges(data = plotdata,aes(linetype = as.factor(same.community)),colour= "grey70", alpha=0.5,size=edgesize,curvature = 0.2,
    #            arrow = arrow(length = unit(0, "pt"), type = "closed", angle=90)) +
    geom_edges(data = plotdata,colour= "grey70", alpha=0.5,size=edgesize,curvature = 0,
               arrow = arrow(length = unit(0, "pt"), type = "closed", angle=90)) +
    guides(linetype = "none")

    # geom_edges(data = subset(plotdata, changetype %in% "Other"),
    #            aes(color= changetype), alpha=0.5,size=edgesize,curvature = 0.2,
    #            arrow = arrow(length = unit(0, "pt"), type = "closed", angle=90)) +
    # geom_edges(data = subset(plotdata, changetype %in% "Activation"),
    #            aes(color= changetype), alpha=0.5,size=edgesize,curvature = 0.2,
    #            arrow = arrow(length = unit(3, "pt"), type = "closed")) +
    # geom_edges(data = subset(plotdata, changetype %in% "Inhibition"),
    #            aes(color= changetype), alpha=0.5,size=edgesize,curvature = 0.2,
    #            arrow = arrow(length = unit(3, "pt"), type = "closed", angle=90)) +
    # scale_color_manual(values=interaction_colors)+
    # labs(color = "Gene Interaction")

  #print(table(plotdata$mutated))
  nodeplot <- edgeplot +
    new_scale_color()+
    # geom_nodes(aes(fill=logFC, alpha=pval_binary),
    # geom_nodes(data = subset(plotdata, mutated %in% "Inhibition"),
    geom_nodes(aes(fill=fillvar, shape=mutated, alpha=pval_binary, color=pval_binary),stroke = 2,#,size = strength),
               # show.legend = showFClegend) +
               size = nodesize,show.legend = showFClegend) +
    guides(size = "none") +
    # color = "grey50", size = nodesize, lwd=2) +
    # shape = 21, color = "grey50") +#, size = nodesize, lwd=2) +
    scale_color_manual(values=outline_vals)+
    scale_shape_manual(values=shape_vals) +
    labs(color = "Mutation Type\n(outline color)",
         fill = "Gene Expression\n(Z-Score)",
         shape = "Mutations") +
    guides(#title.theme = element_text(size=14,face = "bold"),label.theme = element_text(size=12),
           color=guide_legend(title.theme = element_text(size=14,face = "bold"),label.theme =  element_text(size=12),order = 1),
           shape=guide_legend(title.theme = element_text(size=14,face = "bold"),label.theme =  element_text(size=12),order = 2,keyheight=0.5,default.unit="inch"),
           fill=guide_colorbar(title.theme = element_text(size=14,face = "bold"),label.theme =  element_text(size=12),order = 3)
           ) +
    scale_fill_distiller(palette = mypal,na.value = "grey80", limit=color_limits, n.breaks=5, direction = -1) +
    # scale_fill_gradientn(colors=fc_colors,values=fc_colors.values,na.value = "grey80") +
    scale_alpha_manual(values=alpha_vals, guide="none")
  # nodeplot
  #print(labelNamesColors)
  annotatedplot <- nodeplot +
    new_scale_color()+
    #geom_mark_hull(expand=0.01,aes(color=as.factor(membership),fill=(membership)))+#
    #ggalt::geom_encircle(expand=0.01,aes(color=as.factor(membership))) +
    geom_nodetext_repel(aes( label = name,color = isTFColor),#, colour=isTFColor,
                        # geom_nodetext(aes( label = vertex.names, color=logFC ),
                        # fontface = "bold",size=1, color="steelblue") +
                        # fontface = "bold",size=1, show.legend = F) +
                        fontface = "bold",size=nodesize/2, show.legend = TRUE,box.padding = unit(1, "lines")) +
                        # fontface = "bold",size=1, show.legend = T) +
    # scale_colour_gradient2(low="grey70",mid="black",high="grey70") +
    # geom_nodetext(aes(label = name,color=isTFColor))+
    scale_colour_manual(values=labelNamesColors) +
    labs(color = "Is TF\n(label color)") +
    guides(color=guide_legend(title.theme = element_text(size=14,face = "bold"),label.theme = element_text(size=12),order = 3)) +
    #guides(color="none") +
    #theme(legend.title = element_text(size=14,face = "bold"),legend.text = element_text(size=12))+
    #ggtitle(paste0(length(V(fungin_graph))," genes plotted")) +
    theme_blank()

  # annotatedplot
  if (! is.null(savename)) {
    ggsave(annotatedplot, filename = savename, width=8, height=8)
  } #else {
  # print(myplot)
  # }


  return(annotatedplot)

}
fungin_plot_interactive <- function(fungin_graph,fill_var="mutated",marker_size=40) {
  # browser()
  plotdata <- fungin_make_plotdata(fungin_graph)
  nodedata <- plotdata$nodes
  edgedata <- plotdata$edges
  require(plotly)

  node_plot <- plot_ly(data=nodedata,
                       x = ~Xn, y = ~Yn,
                       mode = "markers",
                       text = nodedata$name, hoverinfo = "text",
                       marker = list(
                         color = nodedata$color.background,
                         shape = nodedata$shape,
                         size = marker_size,
                         opacity=0.5,
                         line = list(
                           # color = 'red2',
                           color = nodedata$color.background,
                           width = 5
                         )#,
                         # selected = list(maker=list(opacity=0.1, color="red",size=20))
                       )
  )



  network <- node_plot %>%
    add_annotations(
      data=edgedata,
      x = ~startx,
      y = ~starty,
      ax = ~endx,
      ay = ~endy,
      text = "",
      showarrow = TRUE,
      arrowcolor = ~color,
      arrowhead = ~arrowType,
      arrowsize = 1,
      arrowwidth = ~width,
      standoff = marker_size*0.9,
      startstandoff = marker_size*0.9,
      xref = "x",
      yref = "y",
      axref="x",
      ayref="y"
    )

  a <- list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)
  mynetwork <- network %>% layout(xaxis = a, yaxis = a, dragmode="select",hovermode="closest")

  return(mynetwork)

}


make_mut_summary <- function(maf) {
  require(maftools)
  if (class(maf)=="character") {
    if (file.exists(maf)) {
      mymaf <- read.maf(maf)
    } else {
      stop(paste0("MAF file does not exist: ",maf))
    }
  } else if (class(maf)=="MAF") {
    mymaf <- maf
  } else {
    stop("Argument 'maf' must be a path to a MAF file or a maftools object.")
  }
  # browser()
  gene_mut_summary <- mymaf@gene.summary
  # flag_genes <- c("TTN","MUC16","OBSCN","AHNAK2","SYNE1","FLG","MUC5B","DNAH17","PLEC","DST","SYNE2","NEB","HSPG2","LAMA5","AHNAK","HMCN1","USH2A","DNAH11","MACF1","MUC17","DNAH5","GPR98","FAT1","PKD1","MDN1","RNF213","RYR1","DNAH2","DNAH3","DNAH8","DNAH1","DNAH9","ABCA13","SRRM2","CUBN","SPTBN5","PKHD1","LRP2","FBN3","CDH23","DNAH10","FAT4","RYR3","PKHD1L1","FAT2","CSMD1","PCNT","COL6A3","FRAS1","FCGBP","RYR2","HYDIN","XIRP2","LAMA1")
  # gene_mut_summary <- gene_mut_summary[!gene_mut_summary$Hugo_Symbol %in% flag_genes,]
  num_samples <- as.integer(mymaf@summary$summary[mymaf@summary$ID=="Samples"])
  gene_mut_summary$AlteredFraction <- gene_mut_summary$AlteredSamples/num_samples

  # match_idx <- match(V(my_interaction_graph)$name,gene_mut_summary$Hugo_Symbol,nomatch = 0)
  # mut_status <- gene_mut_summary[match_idx,c("Hugo_Symbol",	"AlteredSamples",	"AlteredFraction")]
  # gene_mut_summary$mutated <- ifelse(gene_mut_summary$AlteredSamples >= min_mutated_samples,
  #                                 paste0("Mutated in ", min_mutated_samples, " or more samples"),
  #                                 NA)

  # return_df <- data.frame(name=gene_mut_summary$Hugo_Symbol, mutated=gene_mut_summary$mutated, stringsAsFactors = F)
  return(gene_mut_summary)
}


gather_FI_network_data <- function(dataDir=file.path("data","fungin"), rds_file=NULL, fromScratch=F) {
  require(igraph)
  if (is.null(rds_file)) {
    rds_file <- file.path(dataDir, "fungin.reactomefi.Rds")
  }

  if (!file.exists(rds_file) | fromScratch ) {
    require(openxlsx)
    if (! dir.exists(dataDir)) {
      dir.create(dataDir, recursive = T)
    }

    print(paste0("Getting FI data..."))
    reactome_FI_url="https://reactome.org/download/tools/ReatomeFIs/FIsInGene_020720_with_annotations.txt.zip"
    reactome_FI_file <- file.path(dataDir,basename(reactome_FI_url))
    if (!file.exists(reactome_FI_file)) {
      download.file(reactome_FI_url,reactome_FI_file)
    }
    # browser()
    print(paste0("Reading FI data..."))
    reactome_FIs <- read.table(unz(reactome_FI_file, gsub(".zip","",basename(reactome_FI_file))), sep="\t", header=T, stringsAsFactors = F)
    reactome_FIs$predicted <- grepl("predicted",reactome_FIs$Annotation)

    change_factor <- list(
      "|-|"= c(-1,-1),
      "<-|" = c( 1,-1),
      "|->" = c(-1, 1),
      "-|"  = c( 0,-1),
      "|-"  = c(-1, 0),
      "<->" = c( 1, 1),
      "->"  = c( 0, 1),
      "<-"  = c( 1, 0),
      "-"   = c( 0, 0)
    )

    include_columns=c("Annotation","predicted")

    print(paste0("Processing FI data..."))
    ## Here we expand all the interactions so that interaction direction is from Column 1 to Column 2
    ## This takes a minute or two
    ## Gene interaction list
    gene_int_list <- apply(reactome_FIs, 1, function(curr_data){
      genes=curr_data[c("Gene1","Gene2")]
      change_vec=change_factor[[curr_data["Direction"]]]

      expanded_df <- cbind(
        do.call(rbind,lapply(change_vec, function(x){
          return_val=genes
          if(x<0){return_val=rev(genes)}
          return(return_val)
        })),
        change=change_vec,
        matrix(rep(curr_data[include_columns],each=2),
               nrow = length(change_vec),
               ncol=length(include_columns),
               dimnames = list(c(1:length(change_vec)), include_columns))
      )

      return(data.frame(expanded_df, stringsAsFactors = F))
    })

    ## Bind into data frame
    gene_int_mat <- do.call(rbind, gene_int_list)
    gene_interactions <- unique(gene_int_mat)





    print(paste0("Geting miRNA target data..."))
    ## Add miRNA target data
    mirna_data_url="http://mirtarbase.cuhk.edu.cn/cache/download/8.0/miRTarBase_MTI.xlsx"
    mirna_data_file=file.path(dataDir, basename(mirna_data_url))
    if (!file.exists(mirna_data_file)) {
      download.file(url=mirna_data_url,destfile = mirna_data_file)
    }

    print(paste0("Reading miRNA target data..."))
    mirna_data.raw <- read.xlsx(mirna_data_file)
    mirna_data <- mirna_data.raw[mirna_data.raw$`Species.(Target.Gene)`=="Homo sapiens",]
    mirna_data <- mirna_data[mirna_data$Support.Type=="Functional MTI",]

    print(paste0("Processing miRNA target data..."))
    mirna_interactions <- data.frame(Gene1=mirna_data$miRNA,
                                     Gene2=mirna_data$Target.Gene,
                                     change=-1,
                                     Annotation=mirna_data$Support.Type,
                                     stringsAsFactors=F)
    mirna_interactions <- unique(mirna_interactions)
    mirna_interactions$predicted <- F
    # mirna_interactions$Annotation <- mirna_interactions$Support.Type
    mirna_interactions$data_source <- "mirTarBase"

    if (! "data_source" %in% colnames(gene_interactions)) {
      gene_interactions$data_source="Reactome_FI"
    }
    gene_interactions <- rbind(gene_interactions, mirna_interactions)

    print(paste0("Combining interactions data..."))
    gene_interactions$change <- as.numeric(gene_interactions$change)
    gene_interactions$changetype <- ifelse(gene_interactions$change == 0, "Other",ifelse(gene_interactions$change < 0, "Inhibition","Activation"))

    full_interaction_network <- igraph::graph_from_data_frame(gene_interactions)

    print(paste0("Saving data..."))
    save(full_interaction_network, file = rds_file)
  } else {
    # load(rds_file)
  }

  return(rds_file)
}




gather_stringdb_data <- function(dataDir=file.path("data","fungin"), rds_file=NULL,
                                 fromScratch=F, genome="hg38",
                                 stringdb_version="11") {

  species_taxids <- c(hg19=9606, hg38=9606, mm10=10090)
  if (! genome %in% names(species_taxids)) {
    stop(paste0("'genome' must be one of '", paste0(names(species_taxids), collapse=","),"'"))
  }
  taxid=species_taxids[genome]
  if (is.null(rds_file)) {
    rds_file <- file.path(dataDir, paste0("fungin.stringdb.",taxid,".Rds"))
  }

  if (!file.exists(rds_file) | fromScratch ) {
    require(STRINGdb)
    if (! dir.exists(dataDir)) {
      dir.create(dataDir, recursive = T)
    }

    print(paste0("Getting STRINGdb data..."))
    # browser()
    # mygenomeinfo <- detect_maf_genome(mafobj_strict)
    string_db <- STRINGdb$new(version=stringdb_version, species=taxid, score_threshold=0, input_directory=dataDir)# ,
    #STRINGdb$new(version="11", species=9606,score_threshold=ppiEdgeThreshold, input_directory="")
    # backgroundV = unique(c(mafobj_strict@maf.silent$Hugo_Symbol, mafobj_strict@data$Hugo_Symbol)))
    string_ppi_graph <- string_db$get_graph()

    #require(org.Hs.eg.db)
    #allgenes <- unique(keys(org.Hs.eg.db,keytype = "SYMBOL"))
    string_gene_file<-file.path(dataDir,"human.name_2_string.tsv.gz")
    print(string_gene_file)
    download.file("https://string-db.org/mapping_files/STRING_display_names/human.name_2_string.tsv.gz",destfile = string_gene_file)
    allgenes <- read.table(string_gene_file, sep="\t", header = F, stringsAsFactors = F)
    string_id_mapfile <- file.path(dataDir,"stringdb_id_map.txt")
    if (!file.exists(string_id_mapfile)) {
      string_id_map <- string_db$map(data.frame(gene=allgenes$V2, stringsAsFactors = F), "gene")
      string_id_map <- string_id_map[!is.na(string_id_map$STRING_id),]
      write.table(string_id_map, file=string_id_mapfile, sep="\t", row.names=F,quote = F)
    } else {
      string_id_map <- read.table(string_id_mapfile, sep="\t", header = T, stringsAsFactors = F)
    }

    nodenames <- names(V(string_ppi_graph))
    id2symbol <- string_id_map$gene[match(nodenames, string_id_map$STRING_id)]
    id2symbol <- ifelse(is.na(id2symbol), gsub(paste0(taxid,"\\."),"",names(V(string_ppi_graph))), id2symbol)
    full_interaction_network <- set.vertex.attribute(string_ppi_graph, "name", value=id2symbol)
    save(full_interaction_network, file = rds_file)
  } else {
    # load(rds_file)
  }
  return(rds_file)
}
