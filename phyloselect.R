#!/usr/bin/env Rscript
 
### Author: Ludovic V. Mallet, PhD, Frank Cerutti
### 2016.06.22
### licence: GPLv3
### Version: Alpha.1
### Garanteed with misatkes. <- Including this one.
### Important notes:

### dependencies:
# library(gtools)
library(ape)
library(getopt) #install.packages("getopt")   #maybe launch it as admin or root if you want it for other users too.
# ### phyloligo.py in the exec PATH 
### EMBOSS in the exec PATH 
### samtools in the exec PATH 

# Functions :


dev.new <- function(width = 7, height = 10) {
  platform <- sessionInfo()$platform
  if (grepl("linux",platform)) {
    x11(width=width, height=height) 
  } else if (grepl("pc",platform)) {
    windows(width=width, height=height)
  } else if (grepl("apple", platform)) {
    quartz(width=width, height=height)
  }
}

get_children<-function(tree,node){
  # Get all children of a given node
  tree$edge[,2][tree$edge[,1]==node]
}

get_all_children<-function(tree,node){
  # Retrieve all nodes (internal and leaves) from a given node
  edges<-tree$edge
  file=c()
  nt=c()
  file=c(file,node)
  while(length(file)>0){
    tmp=edges[,2][edges[,1]==file[1]]
    file=c(file,tmp)
    file<-file[-1]
    nt=c(nt,file[1])
  }
  return(sort(nt))
}

get_all_edges<-function(phy,node){
  # Retrieve all edges from a given node
  nodes=get_all_children(phy,node)
  res=lapply(nodes,function(x,edges){
   which((x==edges[,1] | x==edges[,2]))
  },edges=phy$edge)
  return(sort(unique(unlist(res))))
}

get_parent<-function(tree,node){
  # Get parent of a given node
  tree$edge[,1][tree$edge[,2]==node]
}

get_all_leaves<-function(tree,node){
  # Return all leaf considering a given tree and node
  tips=tree$tip.label
  edges=tree$edge
  
  # retrieve specific node-corresponding leaves
  file=c()
  nt=c()
  file=c(file,node)  
  
  while(length(file)>0){
    tmp=edges[,2][edges[,1]==file[1]]
    if(length(tmp)>0){
      file=c(file,tmp)
      file<-file[-1]
    }else{
      nt=c(nt,tips[file[1]])
      file<-file[-1]
    }
  }
  sort(nt)
}

clade_select<-function(phy, title = TRUE, subbg = "white", return.tree = FALSE,edge_size=rep(0.5,nrow(phy$edge)),edge_label=rep("",nrow(phy$edge)),lim_tip=125,...){
  # Allow user to select specific clades(s) on contigs tree, corresponding sequences could be exported in a fasta file
  lastPP <- get("last_plot.phylo", envir= .PlotPhyloEnv)
  devmain <-dev.cur()
  restore <- function(){
    dev.set(devmain)
    assign("last_plot.phylo", lastPP, envir = .PlotPhyloEnv)
  }
  on.exit(restore())
  
  #NEW<- TRUE
  cat("Click close to a node. Right-click to exit.\n")
  
  subtree=NaN
  
  edge_size_tree<-edge_size
  edge_size<-log(edge_size+1,base=1.5)
  
  user_decision=1
  subset=1
  while (user_decision != "E"){
    x <- identify.phylo(phy, quiet=TRUE)
    if (is.null(x))
      return(invisible(NULL))
    else{
      x <- x$nodes
      
      orig_edges_id=get_all_edges(phy,x)
      
      dev.set(devmain)
      
      edge_selected_color=rep("black",length(phy$edge))
      edge_selected_color[orig_edges_id]<-"red"
      
      plot(phy,use.edge.length=lastPP$use.edge.length,type=lastPP$type,show.tip.label=lastPP$show.tip.label,edge.width=edge_size_tree,edge.color=edge_selected_color)
      edgelabels(text=edge_label,adj=c(0.5,-0.5),frame="none",font=2,cex=0.5,col = edge_selected_color)
      
      if (is.null(x))
        cat("Try again!\n")
      else{
        #if (NEW){
        dev.new()
        par(bg=subbg)
        devsub <- dev.cur()
        #NEW <- FALSE
        #}
        #else dev.set(devsub)
        tr <<- extract.clade(phy,x)
        
        cex_tip=(1/(1+2*log(base=10,tr$Nnode))*1+log(base=10,2-(tr$Nnode/phy$Nnode)))
        
        if(!is.null(edge_size)){
          v<-0.1+edge_size[orig_edges_id]
          k<-as.integer(tree$Nnode)/length(get_all_children(phy,x))
          k<-1
          es<-v*k
          edge.width<-es
          if(tr$Nnode<=lim_tip){
            plot(tr,edge.width=es,cex=cex_tip,...) 
          }else{
            plot(tr,edge.width=es,show.tip.label=F,...)
          }
        }else{
          if(tr$Nnode<=lim_tip){
            plot(tr,cex=cex_tip,...) 
          }else{
            plot(tr,show.tip.label=F,...)
          }
        }
        if(!is.null(edge_label)){
          edgelabels(text=edge_label[orig_edges_id],adj=c(0.5,-1),cex=0.5,frame="none",font=2)
        }
        if (is.character(title))
          title(title)
        else if (title) {
          tl <- if (is.null(phy$node.label)) 
            paste("From node #", x, sep = "")
          else paste("From", phy$node.label[x - Ntip(phy)])
          title(tl)
        }
        if(!is.null(tr) & return.tree==TRUE)
          subtree=tr
          cat("Do you want to export this contig subset as fasta? (y)es, export it. (n)o, discard this selection and let me select another node. (E)xit, I'm done with the selection and exports.\n")
          user_decision <- readLines(file("stdin"),n=1L)
          # print(line)
          dev.off(devsub)
          if (user_decision == "y")
            #return(subtree)
            cmd=paste("samtools faidx ",assembly," '",paste(paste(tr$tip.label,collapse="' '",sep=""),"'",sep=""), " > ",outfile,"_",subset,".fa",sep="" )
            print(cmd)
            system(cmd)
#             plots à join colored mapping  to each selected clade:
            pdf(file=paste(outfile,"_",subset,".pdf",sep=""),width=16,height=9)
                plot(phy,use.edge.length=lastPP$use.edge.length,type=lastPP$type,show.tip.label=lastPP$show.tip.label,edge.width=edge_size_tree,edge.color=edge_selected_color)
                edgelabels(text=edge_label,adj=c(0.5,-0.5),frame="none",font=2,cex=0.5,col = edge_selected_color)
            dev.off()
            if (opt[["verbose"]]) print(paste(date(), "User exported a clade stored in ",paste(outfile,"_",subset,".fa",sep="")))
            subset=subset+1
            
          if (user_decision != "E")
            restore()
            plot(phy,use.edge.length=lastPP$use.edge.length,type=lastPP$type,show.tip.label=lastPP$show.tip.label,edge.width=edge_size_tree)
            edgelabels(text=edge_label,adj=c(0.5,-0.5),frame="none",font=2,cex=0.5)
      }
    }
  }
}
#  interactive_mode()

# print("toto")

# clade_select = function (){
#     user_decision=1
#     subset=1
#     while (substr(user_decision, 1, 1) != "E"){
#         g=tree_zoom(tree,return.tree=TRUE,type='c',use.edge.length=branchlength,cex=0.5,font=2,edge_size=edge_size,edge_label=edge_lab)
#         user_decision = readline("Do you want to export this contig subset as fasta? (y)es, export it. (n)o, discard this selection and let me select another node. (E)xit, I'm done with the selection and exports.")
#         if (substr(user_decision, 1, 1) == "y"){
#             system(paste("samtools faidx ",assembly," ",paste(g$tip.label, collapse=" "), " > ",outfile,"_",subset,".fa",sep="" ))
#             subset=subset+1
#         }
#     }
# 
#     
# }

outersect <- function(x, y, ...) {
  # Extraire les valeurs qui diffèrent entre plusieurs listes 
  big.vec <- c(x, y, ...)
  duplicates <- big.vec[duplicated(big.vec)]
  setdiff(big.vec, unique(duplicates))
}


interactive_mode=function(){
    library(ape)
    X11(width=12,height=10)
    edge_size=edge_size/sum(edge_size)*100*40
    plot.phylo(tree,use.edge.length=branchlength,type="c",show.tip.label=FALSE,edge.width=edge_size)
    edgelabels(text=edge_lab,adj=c(0.5,-0.5),frame="none",font=2,cex=0.5)
    clade_select(tree,return.tree=TRUE,type='c',use.edge.length=branchlength,font=2,edge_size=edge_size,edge_label=edge_lab)
}

spec <- matrix(c(
  'matrix'         , 'i', 1, "character", "all-by-all contig distance matrix, tab separated (required)",
  'assembly'         , 'a', 1, "character", "multifasta file of the contigs (required)",
  'tree_method'     , 't', 2, "character", "tree building method (optional), by default NJ. No other option implemented ATM",
  'dump_R_session'     , 'd', 0, "logical", "Should the R environment saved for later exploration? the filename will be generated from the outfile parameter or its default value",
  'max_perc'     , 'g', 2, "double", "max edge assembly length percentage displayed (%)",
  'min_perc'     , 'l', 2, "double", "min edge assembly length percentage displayed (%)",
  'keep_perc'           , 'k', 2, "double",   "ratio of out-of-range percentages to display (%)",
  'outfile'     , 'o', 2, "character", "outfile name, default:phyloligo.out",
  'branchlength'           , 'b', 0, "logical",   "display branch length",
  'verbose'           , 'v', 0, "logical",   "say what the program do. Not implemented yet.",
  'help'           , 'h', 0, "logical",   "this help"
),ncol=5,byrow=T)


opt = getopt(spec);
# [[""]]
if (!is.null(opt[["help"]]) || is.null(opt[["matrix"]])) {
  cat(paste(getopt(spec, usage=T),"\n"));
}

if (is.null(opt[["verbose"]])) {
  opt[["verbose"]]=FALSE;
}

outfile = ifelse(is.null(opt[["outfile"]]), outfile <- "phyloligo.out" , outfile <-opt[["outfile"]])
branchlength = ifelse(is.null(opt[["branchlength"]]), branchlength <- "FALSE" , branchlength <-"TRUE")
keep_perc = ifelse(is.null(opt[["keep_perc"]]), keep_perc <- 5 , keep_perc <-opt[["keep_perc"]])
min_perc = ifelse(is.null(opt[["min_perc"]]), min_perc <- 0.5 , min_perc <-opt[["min_perc"]])
max_perc = ifelse(is.null(opt[["max_perc"]]), max_perc <- 30 , max_perc <-opt[["max_perc"]])

if ( system("which infoseq",intern=TRUE) == ""){
  print("Please Install EMBOSS")
}
stopifnot(system("which infoseq",intern=TRUE,ignore.stdout = TRUE, ignore.stderr = TRUE) == "")

# print("parameters:")
# print(paste(names(opt),opt[names(opt)]))


# cat(paste(getopt(spec, usage=T),"\n"));

# dist_matrix_file="~/Documents/Ludovic/data/mat_th.out"
# assembly="~/Documents/Ludovic/data/TH12-rn_prefix.fna"

dist_matrix_file = opt[["matrix"]]
tree_method = opt[["tree_method"]]
assembly = opt[["assembly"]]
branchlength = opt[["branchlength"]]







if (opt[["verbose"]]) print(paste(date(), "Reads matrix file"))
dist_matrix = as.matrix(read.delim(dist_matrix_file,sep="\t", header = FALSE))




if (opt[["verbose"]]) print(paste(date(), "Format data"))
# labels = system(paste("infoseq -auto -nocolumns -only -noheading -name ",assembly),intern=TRUE)
labels = system(paste("grep -Po \">[^ ]+\" ",assembly,"|sed 's/>//' "),intern=TRUE)
lengths = as.numeric(system(paste("infoseq -auto -nocolumns -only -noheading -length ",assembly),intern=TRUE))
colnames(dist_matrix) = labels


if (opt[["verbose"]]) print(paste(date(), "Building NJ tree"))
tree = nj(dist_matrix)

# Vector containing contigs size with corresponding name

vlen=lengths
names(vlen)<-labels

# Cumative size for tree edges

edge_size<-unlist(apply(as.matrix(tree$edge),1,function(x){
  sum(vlen[get_all_leaves(tree,x[2])])
}))

# Percentage of cumulative contigs length for tree edges

edge_perc<-unlist(apply(as.matrix(tree$edge),1,function(x){
  #res=round(sum(vlen[get_all_leaves(tree,x[2])])/sum(vlen)*100,digits = 0)
  res=sum(vlen[get_all_leaves(tree,x[2])])/sum(vlen)*100
}))

# Edge labels (filtering percentages to display)

edge_lab=unlist(lapply(edge_perc,function(x){
  paste(round(x,digits=0),"%",sep="")
}))

avoidable_edge_perc=outersect(
  which(!(min_perc<=edge_perc & edge_perc<=max_perc)), # edge_perc avoidable
  sample(edge_perc, round((keep_perc/100)*length(edge_perc)), replace=FALSE) # edge_perc kept
)

edge_lab[avoidable_edge_perc]<-""

#edge_lab[which(tree$edge[,2] %in% 1:length(tree$tip.label))]<-""

# for(i in 1:nrow(tree$edge)){
#   current_node=tree$edge[i,2]
#   
#   root_id=length(tree$tip.label)+1
#   
#   if(tree$edge[i,1]!=root_id){
#     parent=get_parent(tree,current_node)
#     parent_edge_id=which(tree$edge[,2]==parent)
#     
#     calc=as.numeric(edge_perc[parent_edge_id])-as.numeric(edge_perc[i])
#     print(calc)
#     if(calc<1.0){
#       edge_lab[i]<-""
#     }
#     
#   }
#   
#   old_node=current_node
# }

# Plot the tree with cumulative contig length (width of edges) and cumalative contig length percentage (edge label)

#colfunc<-colorRampPalette(c("red","orange","yellow","green"))

# Remove label for terminal edges


if (opt[["verbose"]]) print(paste(date(), "Computing Tree plot"))
X11(width=12,height=10) # external display when script is launched with Rscript command
# par(ljoin = 1, lend = 1)

edge_size=edge_size/sum(edge_size)*100*40

plot(tree,use.edge.length=branchlength,type="c",show.tip.label=FALSE,edge.width=edge_size)
#plot(tree,use.edge.length=FALSE,type="c",show.tip.label=FALSE,edge.width=edge_size/sum(edge_size)*100*15,edge.color = colfunc(100)[round(edge_perc)])
edgelabels(text=edge_lab,adj=c(0.5,-0.5),frame="none",font=2,cex=0.5)

if (opt[["verbose"]]) print(paste(date(), "Indexing fasta file"))
###manually selecting the subtree of the putative contaminant:
system(paste("samtools faidx",assembly))


if (opt[["verbose"]]) print(paste(date(), "Dumping the R environment"))

if (!is.null(opt[["dump_R_session"]])) {
  save.image(file=paste(outfile,"_phyloligo_dump.RData",sep=""));
}

if (opt[["verbose"]]) print(paste(date(), "Entering interactive mode of clade selection to constitute learning sets for contalocate"))
clade_select(tree,return.tree=TRUE,type='c',use.edge.length=branchlength,font=2,edge_size=edge_size,edge_label=edge_lab)

# Export clade-corresponding contig in fasta format

# paste((g$tip.label),collapse= " ")

if (opt[["verbose"]]) print(paste(date(), "Dumping the R environment. To reuse it in the future, open R, load this file with load(\"file\") and call the function interactive_mode(). Beware that files and former selections might get overwritten."))
if (!is.null(opt[["dump_R_session"]])) {
  save.image(file=paste(outfile,"_phyloligo_dump.RData",sep=""));
}


# interactive_mode()
