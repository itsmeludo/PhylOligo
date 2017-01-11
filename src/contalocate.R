#!/usr/bin/Rscript

### Author: Ludovic V. Mallet, PhD
### 2016.03.22
### licence: GPLv3
### Version: Alpha.2
### Garanteed with misatkes. <- Including this one.
### Important notes:


### dependencies:
library(gtools)
library(getopt) #install.packages("getopt")   #maybe launch it as admin or root if you want it for other users too.
### Kount.py in the exec PATH or where this script is launched


### Standard parameters
win_step=500 #bp: sliding step
win_size=5000 #bp:window size for composition computation
dist="KL"
# 
# genome_fasta="/home/ludovic/taff/projects/inra/sequences/TH12_prefix.fna"
# host_sample_fasta="/home/ludovic/taff/projects/inra/sequences/TH12_prefix.fna"
# conta_sample_fasta="TH12_positions_Burk.GFF.gdna.fasta"


spec <- matrix(c(
        'genome'         , 'i', 1, "character", "file from fastq-stats -x (required)",
        'host_learn'     , 'r', 2, "character", "input gc content file (optional)",
        'conta_learn'    , 'c', 1, "character", "output filename (optional)",
        'win_step'       , 't', 2, "int", "output filename (optional)",
        'win_size'       , 'w', 2, "int", "output filename (optional)",
        'outputdir'      , 'W', 1, "character", "path to outputdir directory",
        'dist'           , 'd', 2, "character", "Divergence metric used to compare profiles: (KL), JSD or Eucl",
        'manual_threshold' , 'm', 0, "logical", "You will be asked to manually set the thresholds",
        'help'           , 'h', 0, "logical",   "this help"
),ncol=5,byrow=T)

#         

opt = getopt(spec);
# [[""]]
if (!is.null(opt[["help"]]) || is.null(opt[["genome"]])) {
    cat(paste(getopt(spec, usage=T),"\n"));
    quit(save="no", status = 1)
}
# print("parameters:")
# print(paste(names(opt),opt[names(opt)]))


# cat(paste(getopt(spec, usage=T),"\n"));
genome_fasta = opt[["genome"]]
working_dir = opt[["outputdir"]]
conta_sample_fasta = opt[["conta_learn"]]
host_sample_fasta = ifelse(is.null(opt[["host_learn"]]), test <- "" , test <-opt[["host_learn"]])
dist = ifelse(is.null(opt[["dist"]]), test <- dist , test <-opt[["dist"]])
win_step = ifelse(is.null(opt[["win_step"]]), test <- win_step , test <-opt[["win_step"]])
win_size = ifelse(is.null(opt[["win_size"]]), test <- win_size , test <-opt[["win_size"]])


### compute profiles:
data=list()
if ( is.null(opt[["host_learn"]])) {
  cmd=paste("bash -lic \"ionice -c2 -n3 Kount.py -i ",genome_fasta ," -W ",working_dir," -w ",win_size," -t ",win_step," -c ",conta_sample_fasta," -d  ",dist," -n 0.5\"",sep="")
  print(cmd)
#   system(cmd)
  input_path=paste(working_dir, basename(genome_fasta),".mcp_hostwindows_vs_wholegenome_",dist,".dist",sep="")
  data[["host"]]=read.delim(file=input_path, header=F)
}else{
  cmd=paste("bash -lic \"ionice -c2 -n3 Kount.py -i ",genome_fasta ," -W ",working_dir," -w ",win_size," -t ",win_step," -c ",conta_sample_fasta," -r ",host_sample_fasta," -d ",dist," -n 0.5\"",sep="")
  print(cmd)
#   system(cmd)
  input_path=paste(working_dir,basename(genome_fasta),".mcp_hostwindows_vs_host_",basename(host_sample_fasta),"_",dist,".dist",sep="")
  data[["host"]]=read.delim(file=input_path, header=F)
}
input_path=paste(working_dir, basename(genome_fasta),".mcp_hostwindows_vs_conta_",basename(conta_sample_fasta),"_",dist,".dist",sep="")
data[["conta"]]=read.delim(file=input_path, header=F)

print("end comp")

if (! is.null(opt[["manual_threshold"]])) {
  
  ### Ask the trusty human to set the thresholds:

  threshold_conta=0
  repeat{
    plot(density(data[["conta"]][,4],na.rm=TRUE),xlim=c(0,5000),lwd=2)
    abline(v=threshold_conta,col="red")
    new_threshold= ask("Give a different threshold value for the contaminant threshold. Give the same value to confirm it.\n")
    new_threshold <- as.numeric(new_threshold)
    
    if(new_threshold==threshold_conta){
      break
    }
    threshold_conta=new_threshold
  }


  threshold_host=0
  repeat{
    plot(density(data[["host"]][,4],na.rm=TRUE),xlim=c(0,5000),lwd=2)
    abline(v=threshold_host,col="red")
    new_threshold= ask("Give a different threshold value for the host threshold. Give the same value to confirm it.\n")
    new_threshold <- as.numeric(new_threshold)
    
    if(new_threshold==threshold_host){
      break
    }
    threshold_host=new_threshold
  }

}else{

  ### Humans are not worthy to set the thresholds, stats will guess it:


  des_conta=density(data[["conta"]][which(!is.nan(data[["conta"]][,4]) ),4] )
  plot(des_conta,lwd=2)
#   points(x= des_conta[["x"]][which.max(des_conta[["y"]])],y= des_conta[["y"]][which.max(des_conta[["y"]])])
  steep=des_conta[["y"]][seq(which.max(des_conta[["y"]]),0)]
  i=1
  while(steep[i+1]<steep[i]){i=i+1}
  conta_min=(which.max(des_conta[["y"]])-i)
  abline(v=des_conta[["x"]][conta_min],col="blue",lwd=2)
  threshold_conta=des_conta[["x"]][conta_min]
  ask("Please inspect that the automatic threshold for the contaminant was set-up properly.")


  des_host=density(data[["host"]][which(!is.nan(data[["host"]][,4]) ),4] )
  plot(des_host,lwd=2)
#   points(x= des_host[["x"]][which.max(des_host[["y"]])],y= des_host[["y"]][which.max(des_host[["y"]])])
  steep=des_host[["y"]][seq(which.max(des_host[["y"]]),length(des_host[["y"]]))]
  i=1
  while(steep[i+1]<steep[i]){i=i+1}
  host_min=(which.max(des_host[["y"]])+i)
  abline(v=des_host[["x"]][host_min],col="blue",lwd=2)
  threshold_host=des_host[["x"]][host_min]
  ask("Please inspect that the automatic threshold for the host was set-up properly.")
}


### Perform the split over the double threshold

data[["Select_conta"]]=which((data[["conta"]][,4]<= threshold_conta )* (data[["host"]][,4]>= threshold_host)>=1)
data[["windows_conta"]]=data[["conta"]][data[["Select_conta"]],]
# write(file=paste(f,"fenetres_Burk.txt",sep=""),as.character(data[["Fenetres_Burk"]]))


# a=cbind((data[["conta"]][,4]<= threshold_conta ),(data[["host"]][,4]>= threshold_host),((data[["conta"]][,4]<= threshold_conta )* (data[["host"]][,4]>= threshold_host)),(data[["conta"]][,4]<= threshold_conta )* (data[["host"]][,4]>= threshold_host)>=1)
# b=data[["conta"]][data[["Select_conta"]],]



### Regroups contiguous windows into islands and write a GFF file of the positions of the targeted species
# regroup_struct=data[["Select_conta"]]
write(x="##gff-version 2", file = paste(working_dir, basename(genome_fasta),"_contaminant_",basename(conta_sample_fasta),".gff",sep=""),ncolumns=1)
start_index=data[["Select_conta"]][1]
for(i in seq(length(data[["Select_conta"]]))){
#   if(ifelse(is.na(data[["Select_conta"]][i+1]), test <- 0 , test <- data[["Select_conta"]][i+1]) == data[["Select_conta"]][i]+1 | data[["conta"]][i,1] == ifelse(is.na(data[["conta"]][i+1,1]), test <- -1 , test <- data[["conta"]][i+1,1])){ #2 windows have a consecutive index number, so they are physically touching. data[["conta"]][i,1] describe the contig, and regions are only spanning one contig max. we assert that ssuccessive indexes are in the same contig to group them. If everything of this if fine, we will then group the windows in a region, by etending the end of it.
  if(ifelse(is.na(data[["Select_conta"]][i+1]), test <- 0 , test <- data[["Select_conta"]][i+1]) == data[["Select_conta"]][i]+1 ){ #2 windows have a consecutive index number, so they are physically touching. data[["conta"]][i,1] describe the contig, and regions are only spanning one contig max. we assert that ssuccessive indexes are in the same contig to group them. If everything of this if fine, we will then group the windows in a region, by etending the end of it.
    end_index=data[["Select_conta"]][i+1]
  }else{ #well, the window index i+1 is not touching the window i, so i is the last of its island, and the island can be written. i+1 is the start of a new island.
    end_index=data[["Select_conta"]][i]
    line=paste(sep="\t", as.character(data[["conta"]][start_index,1]),"SignatureGohtam\tregion",data[["conta"]][start_index,2],data[["conta"]][end_index,c(3)],"\t.\t.\t.")#,c(paste(sep="/", data[["conta"]][start_index:end_index,4])))
#     print(line)
    write(x=line,append=TRUE, file = paste(working_dir, basename(genome_fasta),"_contaminant_",basename(conta_sample_fasta),".gff",sep=""),ncolumns=1)
    start_index=data[["Select_conta"]][i+1]
  }
}

print("Done")

### Booom Done

# image.save(file="test.Rdata")





