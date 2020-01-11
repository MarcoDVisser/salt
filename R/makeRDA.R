setwd("/home/mvisser/symlinks/git/salt/")
bands <- read.csv("./src/bands.csv")
save(bands,file="./data/bands.rda")
