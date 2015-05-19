# This script will:

#Loop through columns in Taxa output file or OTU table to 
#perform kruskal wallis on each and 
#print the taxa names, p values and adjusted p values in a dataframe
#within the function and a .csv file in your working directory

#Then..

#Subset the OTUs from your original table into a new .csv file
#containing only OTUs with significant differences between groups
#And create a set of boxplots displaying the relative abundance of
#each different OTU in your experimental groups

#Read .csv file into the environment as a dataframe
getwd()
setwd("/Users/marymoore/Desktop/MilkAnalysis/RawMilkTaxa")
taxa<-read.csv("Genera_All_Raw_Milks_no_low_seq_samp_1%_forR.csv", header=TRUE)

#Function
kruskal.taxa<-function(taxa, ...){
        df<-as.data.frame(taxa)
        n<-ncol(df)
        v<- n-2
        dat<-as.list(df[,3:n]) #dataframe into a list of vectors named by taxa
        subject<-df[,c(1)] #vector denoting which groups to compare

        #Calculate the kruskal wallis statistics and put them into a list
        result<-vector(mode = "list", length = v)
        for (i in 1:v) {
                result[[i]]<-kruskal.test(dat[[i]],subject)
        }
        
        taxa.name<-names(df[,3:n])
        p.vals<-lapply(result, '[', 'p.value')
        p<-as.numeric(c(do.call("rbind",p.vals)))
        p.corrected<-as.numeric(p.adjust(p, method="fdr"))
        kruskal.vals<-as.data.frame(cbind(taxa.name,p,p.corrected))
        write.csv(kruskal.vals,file="kruskal_wallis_values.csv")
        kruskal.vals
        
        #subset the dataframe just generated
        #       to contain only significant p-values
        kruskal.vals[,3]<-as.numeric(as.character(kruskal.vals[,3]))
        sig<-subset(kruskal.vals, p.corrected<=0.05)
        
        #Use this subset to create a new OTU table containing only
        #       OTUs with significant differences
        sig.otus<-cbind(taxa[,1:2],taxa[,c(colnames(taxa) %in% sig$taxa.name), 
                                        drop=FALSE])
        write.csv(sig.otus,file="otus_with_significant_changes.csv")
        
        #Create a file named "boxplots" within the working directory
        #       write boxplots in pdf for all of the significantly
        #       different taxa
        dir.create(paste(getwd(),"/boxplots", sep=""))
        setwd(paste(getwd(),"/boxplots", sep=""))
        for (i in 3:ncol(sig.otus)){
                pdf(paste(names(sig.otus[i]),".pdf"))
                boxplot(sig.otus[,i] ~ sig.otus[,1], main=names(sig.otus[i]),
                        las=2, cex.axis=0.5)
                dev.off()
        }
        
}