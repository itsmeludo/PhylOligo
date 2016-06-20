#!/usr/bin/Rscript

### Author: Ludovic V. Mallet, PhD
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


spec <- matrix(c(
        'matrix'         , 'i', 1, "character", "all-by-all contig distance matrix, tab separated (required)",
        'assembly'         , 'a', 1, "character", "multifasta file of the contigs (required)",
        'tree_method'     , 't', 2, "character", "tree building method (optional), by default NJ",
        'help'           , 'h', 0, "logical",   "this help"
),ncol=5,byrow=T)

#         

opt = getopt(spec);
# [[""]]
if (!is.null(opt[["help"]]) || is.null(opt[["matrix"]])) {
    cat(paste(getopt(spec, usage=T),"\n"));
}


if ( system("which infoseq",intern=TRUE) == ""){
    print("Please Install EMBOSS")
}
stopifnot(system("which infoseq",intern=TRUE) == "")

# print("parameters:")
# print(paste(names(opt),opt[names(opt)]))


# cat(paste(getopt(spec, usage=T),"\n"));
dist_matrix_file = opt[["matrix"]]
tree_method = opt[["tree_method"]]
assembly = opt[["assembly"]]



dist_matrix = as.matrix(read.delim(dist_matrix_file,sep="\t", header = FALSE))



labels = system(paste("infoseq -auto -nocolumns -only -noheading -name ",opt[["assembly"]]),intern=TRUE)
lengths = as.numeric(system(paste("infoseq -auto -nocolumns -only -noheading -length ",opt[["assembly"]]),intern=TRUE))


colnames(dist_matrix) = labels

tree = nj(dist_matrix)


plot(tree,show.tip.label=FALSE,type="p")


###selection of the contamiant subtree:

g=trex(tree,return.tree=TRUE )


system(paste("samtools faidx",opt[["assembly"]] ))

# paste((g$tip.label),collapse= " ")

system(paste("samtools faidx ",opt[["assembly"]], paste(g$tip.label, collapse=" "), "> out.test" ))
