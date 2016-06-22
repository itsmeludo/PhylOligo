#!/usr/bin/Rscript

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

node_all_leaves<-function(tree,node){
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

spec <- matrix(c(
        'matrix'         , 'i', 1, "character", "all-by-all contig distance matrix, tab separated (required)",
        'assembly'         , 'a', 1, "character", "multifasta file of the contigs (required)",
        'tree_method'     , 't', 2, "character", "tree building method (optional), by default NJ",
        'outfile'     , 't', 2, "character", "tree building method (optional), by default NJ",
        'branchlength'           , 'h', 0, "logical",   "this help"
        'help'           , 'h', 0, "logical",   "this help"
),ncol=5,byrow=T)

#         

opt = getopt(spec);
# [[""]]
if (!is.null(opt[["help"]]) || is.null(opt[["matrix"]])) {
    cat(paste(getopt(spec, usage=T),"\n"));
}

outfile = ifelse(is.null(opt[["outfile"]]), outfile <- "phyloligo.out" , outfile <-opt[["outfile"]])

branchlength = ifelse(is.null(opt[["branchlength"]]), branchlength <- "FALSE" , branchlength <-opt[["branchlength"]])



if ( system("which infoseq",intern=TRUE) == ""){
    print("Please Install EMBOSS")
}
stopifnot(system("which infoseq",intern=TRUE,ignore.stdout = TRUE, ignore.stderr = TRUE) == "")

# print("parameters:")
# print(paste(names(opt),opt[names(opt)]))


# cat(paste(getopt(spec, usage=T),"\n"));
dist_matrix_file = opt[["matrix"]]
tree_method = opt[["tree_method"]]
assembly = opt[["assembly"]]
branchlength = opt[["branchlength"]]

dist_matrix = as.matrix(read.delim(dist_matrix_file,sep="\t", header = FALSE))

labels = system(paste("infoseq -auto -nocolumns -only -noheading -name ",assembly),intern=TRUE)
lengths = as.numeric(system(paste("infoseq -auto -nocolumns -only -noheading -length ",assembly),intern=TRUE))

colnames(dist_matrix) = labels
tree = nj(dist_matrix)

# Vector containing contigs size with corresponding name

vlen=lengths
names(vlen)<-labels

# Cumative size for tree edges

edge_size<-unlist(apply(as.matrix(tree$edge),1,function(x){
  sum(vlen[node_all_leaves(tree,x[2])])
}))

# Percentage of cumulative contigs length for tree edges

edge_perc<-unlist(apply(as.matrix(tree$edge),1,function(x){
  res=round(sum(vlen[node_all_leaves(tree,x[2])])/sum(vlen)*100,digits = 0)
}))

# Plot the tree with cumulative contig length (width of edges) and cumalative contig length percentage (edge label)

colfunc<-colorRampPalette(c("red","orange","yellow","green"))
edge_lab=lapply(edge_perc,function(x){paste(x,"%",sep="")})

X11() # external display when script is launched with Rscript command
# par(ljoin = 1, lend = 1)
plot(tree,use.edge.length=FALSE,type="c",show.tip.label=FALSE,edge.width=edge_size/sum(edge_size)*100*40)
# plot(tree,use.edge.length=FALSE,type="c",show.tip.label=FALSE,edge.width=edge_size/sum(edge_size)*100*15,edge.color = colfunc(100)[round(edge_perc)])
edgelabels(text=edge_lab,adj=c(0.5,-0.5),frame="none",font=2)



###manually selecting the subtree of the putative contaminant:
g=trex(tree,return.tree=TRUE)
plot(g)

# Wait the plot of g

message("Press Return To Continue")
invisible(readLines("stdin", n=1))

system(paste("samtools faidx",assembly))

# paste((g$tip.label),collapse= " ")

system(paste("samtools faidx ",assembly, paste(g$tip.label, collapse=" "), "> ",outfile ))
