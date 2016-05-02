Granges.2.df <- function(grange, midpoint = TRUE){
        #
        # Function that coerces a Granges object into a data frame. If midpoint is TRUE
        # then the midpoint for each interval is added between the start and end columns
        #
        if(!require("dplyr")){
                install.packages("dplyr")
        }
        library(dplyr)
        #   
        df <- as.data.frame(grange)
        #
        #Add midpoint value
        #
        if(midpoint == TRUE){
                temp <- as.numeric(ncol(df))
                df <- df %>%
                        dplyr::mutate( midpoint = (start + end)/2) %>%
                        dplyr::select(c(1:start,midpoint,end:temp))
        }
        df
}


which.index <- function(df,namesvec){
        #
        # Function which takes a dataframe and a character vector of names and returns
        # a vector of the indices of these names
        #
        index <- vector("numeric",length = length(namesvec))
        for(i in 1:length(namesvec)){
                index[i] <- which(names(df)==as.character(namesvec[i]))
        }
        index
}

threshold.classes <- function(df, length, index = c(1,2,3),
                              threshold = rep(0.05, times = 3),
                              greater = rep(TRUE, times = 3),
                              eq = rep(TRUE, times = 3),
                              head = c("class1","class2","class3")){
        #
        # Function which takes a data frame df and creates additional columns 
        # representing a classification of existing variables as being either above or 
        # below a given threshold. The new columns are placed after the column with 
        # index given by length. The columns to be classified are given by a numeric
        # vector index. The thresholds to compare these to are given as threshold. The
        # argument greater sets is the comparison should be greater than or less than,
        # with TRUE meaning greater than for the given comparison. The argument eq
        # sets if the comparison should be strict, or aloow equality, with TRUE 
        # represnting equality being allowed. Head sets the names for the new columns.
        # 
        #
        if(!require("dplyr")){
                install.packages("dplyr")
        }
        library(dplyr)
        #
        if(!is.logical(greater)){
                stop("greater must be logical")
        }
        if(!is.logical(eq)){
                stop("eq must be logical")
        }
        #        
        temp <- as.numeric(ncol(df))
        #
        for (i in 1:length(index)){
                vect <- if(greater[i]){
                        if(eq[i]){
                                df[,index[i]] >= threshold[i]
                        } else {
                                df[,index[i]] > threshold[i]
                        }
                } else {if(eq[i]){
                        df[,index[i]] <= threshold[i]
                } else {
                        df[,index[i]] < threshold[i]
                }
                }
                vect[is.na(df[,index[i]])] <- NA
                df <- df %>%
                        dplyr::mutate(temp = vect)
                names(df)[temp + i] <- head[i]
        }
        #
        df <- df %>%
                dplyr::select(c(1:length,
                                (temp + 1):(temp + length(index)),
                                (length + 1):temp))
        df
}


col.ratio <- function(df,index,split,title){
        #
        # Funtion to take a data frame df and add the ratio of the columns given by the
        # indices in index to the dataframe, adding it after the column with the index
        # given by split, and giving the new column the name given by title
        #
        if(!require("dplyr")){
                install.packages("dplyr")
        }
        library(dplyr)
        #
        temp <- as.numeric(ncol(df))
        data <- df[,index[1]]/df[,index[2]]
        data[which(df[,index[2]]==0)] <- NA
        data[is.na(df[,index[2]])] <- NA
        data[is.na(df[,index[1]])] <- NA
        #
        df <- df %>%
                dplyr::mutate(data_ratio = data) %>%
                dplyr::select(1:split,
                              (temp + 1),
                              (split + 1):temp)
        names(df)[split + 1] <- title
        df
}


df.overlap <- function(df, compare, split, title, strand = TRUE, type = "any",
                       relabel = FALSE){
        #
        # Function that takes a dataframe, df, (with columns for start and end of an 
        # Irange, as well a column for strand) and a Granges object, compare and adds a 
        # column of any overlaps between the ranges. This overlap is added to the data 
        # frame as a new column after the column with index split, and is given the
        # title as given by title. Strand sets whether the strand has to match for an
        # overlap to count. Type specifies that the mechanism inherited from Granges for
        # counting a range in df as overlapping one in compare. If relabel is TRUE, the
        # labels for chromosomes in compare are strpped down to simply be numerical
        #
        #
        if(!require("plyr")){
                install.packages("plyr")
        }
        library(plyr)
        #
        if(!require("dplyr")){
                install.packages("dplyr")
        }
        library(dplyr)     
        #
        if(!require("GenomicRanges")){
                source("https://bioconductor.org/biocLite.R")
                biocLite("GenomicRanges")
        }
        library(GenomicRanges)  
        #
        if(!require("stringr")){
                install.packages("stringr")
        }
        library(stringr) 
        #
        df.gr <- GenomicRanges::makeGRangesFromDataFrame(df)
        temp <- as.numeric(ncol(df))
        #
        if(relabel){
                compare.df <- as.data.frame(compare)        
                regexp <- "[[:digit:]]+"
                from <- unique(grep(regexp,compare.df[,1],value = TRUE))
                to <- stringr::str_extract(from,regexp)
                compare.df[,1] <- plyr::mapvalues(compare.df$seqnames, from = from, to = to)
                compare <- GenomicRanges::makeGRangesFromDataFrame(compare.df,
                                                                   keep.extra.columns = TRUE)
        }
        #
        if(strand){
                df <- df %>%
                        dplyr::mutate(
                                overlap = GenomicRanges::overlapsAny(df.gr,compare,
                                                                     ignore.strand = TRUE,
                                                                     type = type)
                        ) %>%
                        dplyr::select(c(1:split,
                                        overlap,
                                        (split +1):temp))
        } else {
                df <- df %>%
                        dplyr::mutate(
                                overlap = GenomicRanges::overlapsAny(df.gr,compare,
                                                                     ignore.strand = FALSE,
                                                                     type = type)
                        ) %>%
                        dplyr::select(c(1:split,
                                        overlap,
                                        (split +1):temp)) 
        }
        names(df)[split + 1] <- title
        df
}

classes.overlap <- function(df, compare, indexName, split, strand = TRUE, 
                            type = "any", append = FALSE, addName = "",
                            relabel = FALSE){
        #
        # Function that takes a dataframe, df, and a granges object, compare. It then
        # adds columns to df, starting after column given by split. These contain a
        # logical vector indicating whether each bin overlaps the features contained in
        # the column given by indexName of the Granges object. If append is TRUE, then 
        # the names of the columns are modfied by adding the string given as addName 
        # after first adding an underscore. If relabel is TRUE, the labels for 
        # chromosomes in compare are strpped down to simply be numerical
        #
        #
        if(!require("dplyr")){
                install.packages("dplyr")
        }
        library(dplyr)     
        #
        if(!require("GenomicRanges")){
                source("https://bioconductor.org/biocLite.R")
                biocLite("GenomicRanges")
        }
        library(GenomicRanges)  
        #
        comp.df <- as.data.frame(compare)
        types <- as.vector(levels(factor(comp.df[,indexName])))
        ntypes <- types
        if(append == TRUE){
                ntypes <- paste(types,rep(addName, times = length(types)),sep="_")
        }
        
        #
        i <- 1
        for (data in types){
                filt <- comp.df[which(comp.df[,indexName] == data),] %>%
                        GenomicRanges::makeGRangesFromDataFrame()
                tmp.spl <- split + i - 1
                tmp.ttl <- ntypes[i]
                
                df <- df.overlap(df, filt, tmp.spl, tmp.ttl, strand = strand, 
                                 type = type, relabel = relabel)
                i <- i + 1
        }
        df
}

Granges.score.merge <- function(df, compare, title, split, score.index = 6, 
                                bin.width = 200, cutoff.s = FALSE,
                                cutoff.e = FALSE, strand = TRUE){
        #
        # Function that takes a data frame from a granges object, and a Granges object,
        # compare, that is on the same ranges, with ranges containing a score for some
        # measure on that range. The function maps these scores onto df, and calculates
        # the average score for each range in df. It add these scores as a new column
        # called per the argument title, and places it after column given by split. 
        # Score index tells the function where the score will be in the coerced form of
        # compare. cutoff.s and cutoff.e, when set to true, makes the cutoff for ranges
        # in compare overlapping boundaries of ranges in df hard (using interpolation), 
        # for the start and end respectively. Strand decides whther or not to ignore
        # strand (TRUE will ignore strand)
        #
        if(!require("dplyr")){
                install.packages("dplyr")
        }
        library(dplyr)     
        #
        if(!require("GenomicRanges")){
                source("https://bioconductor.org/biocLite.R")
                biocLite("GenomicRanges")
        }
        library(GenomicRanges)
        #
        comp.df <- as.data.frame(compare)
        df.gr <- GenomicRanges::makeGRangesFromDataFrame(df)
        #
        olaps <- if(strand){
                GenomicRanges::findOverlaps(df.gr, compare,
                                            ignore.strand = TRUE) %>%
                        as.data.frame()
        } else {
                GenomicRanges::findOverlaps(df.gr, compare,
                                            ignore.strand = FALSE) %>%
                        as.data.frame()        
        }
        #
        olaps <- olaps %>%
                dplyr::mutate(
                        score = comp.df[olaps[,2],score.index]*comp.df$width[olaps[,2]]
                )
        
        if(cutoff.s){
                ind1 <- comp.df[olaps[,2],"start"] < df[olaps[,1],"start"]
                diff1 <- df[olaps[ind1,1],"start"] - comp.df[olaps[ind1,2],"start"]
                adj1 <- diff1*comp.df[olaps[ind1,2],score.index]
                olaps[ind1,3] <- olaps[ind1,3] - adj1
                rm(ind1,diff1,adj1)
        }
        
        if(cutoff.e){
                ind2 <- comp.df[olaps[,2],"end"] > df[olaps[,1],"end"]
                diff2 <- comp.df[olaps[ind2,2],"end"] - df[olaps[ind2,1],"end"]
                adj2 <- diff2*comp.df[olaps[ind2,2],score.index]
                olaps[ind2,3] <- olaps[ind2,3] - adj2
                rm(ind2,diff2,adj2)
        }
        
        olaps <- olaps %>%
                dplyr::group_by(queryHits) %>%
                dplyr::summarise(sum(score)/bin.width) %>%
                dplyr::ungroup() %>%
                as.data.frame()
        
        temp.r <- nrow(df)
        df <- df %>%
                mutate(temp = rep(NA, times = temp.r))
        #
        temp.c <- ncol(df)
        df[olaps[,1],temp.c] <- olaps[,2]
        colnames(df)[temp.c] <- title
        df <- dplyr::select(df, c(1:split, temp.c, (split + 1):(temp.c - 1)))
        df
}


col.gr.mean <- function(df,index,split,title){
        #
        # Funtion to take a data frame df and add the mean of the columns given by the
        # indices in index to the dataframe, adding it after the column with the index
        # given by split, and giving the new column the name given by title
        #
        if(!require("dplyr")){
                install.packages("dplyr")
        }
        library(dplyr)
        #
        temp <- as.numeric(ncol(df))
        data <- rowMeans(df[,index], na.rm = TRUE)
        #
        df <- df %>%
                dplyr::mutate(temp_mean = data) %>%
                dplyr::select(1:split,
                              (temp + 1),
                              (split + 1):temp)
        names(df)[split + 1] <- title
        df
}

col.prod <- function(df,index,split,title){
        #
        # Funtion to take a data frame df and add the product of the columns given by 
        # the indices in index to the dataframe, adding it after the column with the 
        # index given by split, and giving the new column the name given by title
        #
        if(!require("dplyr")){
                install.packages("dplyr")
        }
        library(dplyr)
        #
        temp <- as.numeric(ncol(df))
        data <- df[,index[1]]*df[,index[2]]
        data[is.na(df[,index[2]])] <- NA
        data[is.na(df[,index[1]])] <- NA
        #
        df <- df %>%
                dplyr::mutate(data_prod = data) %>%
                dplyr::select(1:split,
                              (temp + 1),
                              (split + 1):temp)
        names(df)[split + 1] <- title
        df
}


annotation.builder <- function(data.base,granges,features,index, names, strand = TRUE){
        #
        # df data frame, files = granges files, feautres - name of thing to be overlaid,
        # index - whihc column is this given in
        #
        if(!require("dplyr")){
                install.packages("dplyr")
        }
        library(dplyr)     
        #
        if(!require("GenomicRanges")){
                source("https://bioconductor.org/biocLite.R")
                biocLite("GenomicRanges")
        }
        library(GenomicRanges)
        if(!require("rtracklayer")){
                source("https://bioconductor.org/biocLite.R")
                biocLite("rtracklayer")
        }
        library(rtracklayer)
        #
        labels <- list()
        base <- data.base %>%
                GenomicRanges::reduce()
        for (i in length(features)){
                data <-  granges[[i]]
                if (is.na(features[i]) == FALSE){
                        data <- as.data.frame((data))
                        data <- data[which(data[,index[i]] == features[i]),] %>%
                                GenomicRanges::makeGRangesFromDataFrame()
                }
                #
                extend <- annot.diff(base,data,strand)
                base <- sort(c(base,extend))
        }
        labels[[1]] <- base
        #
        for (i in 1:length(features)){
                data <-  granges[[i]]
                if (is.na(features[i]) == FALSE){
                        data <- as.data.frame((data))
                        data <- data[which(data[,index[i]] == features[i]),] %>%
                                GenomicRanges::makeGRangesFromDataFrame()
                }
                #
                labels[[i+1]] <- annot.diff(data,labels[[i]],strand)
                labels[[i]] <- annot.intersect(data,labels[[i]],strand)
                mcols(labels[[i]]) <- as.factor(names[i])
                rm(data)
        }
        mcols(labels[[length(names)]]) <- as.factor(names[length(names)])
        annot <-  labels[[1]]
        
        for (i in 2:length(names)){
                annot <- sort(c(annot,labels[[i]]))
        }
        names(mcols(annot)) <- "annotation"
        annot
}

annot.intersect <- function(data,base,strand){
        int <- if(strand){
                GenomicRanges::intersect(data,base,ignore.strand=TRUE)
        } else {
                GenomicRanges::intersect(data,base,ignore.strand=FALSE)
        }
        int
}

annot.diff <- function(data,base,strand){
        diff <- if(strand){
                GenomicRanges::setdiff(base,data,ignore.strand =TRUE)
        } else {
                GenomicRanges::setdiff(base,data,ignore.strand =FALSE)   
        }
        diff
}

annot.2.bin <- function(df,annot,split,strand = TRUE){
        #
        # Function that appends the data frame (of a Granges object) with information
        # from a matching genome annotation. This must produce a unique annotation for 
        # each nucleotide. The annotation is added to each window in the data frame df.
        # In cases where the annotation is not unique, the returned column will display
        # non-unique. This appended column is then added after the column given by 
        # split.
        #
        if(!require("dplyr")){
                install.packages("dplyr")
        }
        library(dplyr)     
        #
        if(!require("GenomicRanges")){
                source("https://bioconductor.org/biocLite.R")
                biocLite("GenomicRanges")
        }
        #
        df.gr <- GenomicRanges::makeGRangesFromDataFrame(df)
        temp <- if(strand){GenomicRanges::findOverlaps(df.gr,annot,ignore.strand=TRUE)}
        else {GenomicRanges::findOverlaps(df.gr,annot,ignore.strand=FALSE)}
        rm(df.gr)
        temp <- as.data.frame(temp)
        
        non.unique <- temp %>%
                dplyr::group_by(queryHits) %>%
                dplyr::summarise(overlaps = n()) %>%
                dplyr::filter(overlaps != 1) %>%
                as.data.frame()
        
        Annotation <- vector("character",length = nrow(df))
        df <- dplyr::mutate(df, Annotation = Annotation)
        rm(Annotation)
        temp.l <- as.numeric(ncol(df))
        
        df[non.unique[,1],temp.l] <- "non-unique"
        
        unique <- dplyr::filter(temp, queryHits %in% non.unique[,1] == FALSE)
        df[unique[,1],temp.l] <- as.character(mcols(annot)[unique[,2],1])
        df[,temp.l] <- as.factor(df[,temp.l])
        
        df<- dplyr::select(df,c(1:split,temp.l,(split + 1):(temp.l - 1)))
}

train.partition <- function(df,ratio){
        #
        # Function that takes a data frame and appends it with a column indicating 
        # training and test sets. These are picked at random in proportions given by 
        # ratio
        #
        if(!require("dplyr")){
                install.packages("dplyr")
        }
        library(dplyr)
        #
        df <- df %>% 
                dplyr::mutate(act_id = 1:as.numeric(nrow(.)),
                              rand = runif(nrow(df)))
        cut <- ceiling(nrow(df)*ratio)
        df <- df[order(df[,"rand"]),] %>%
                dplyr::mutate(
                        set = c(rep("training",times = cut),
                                rep("test",times = (nrow(df)-cut)))
                )
        df <- df[order(df[,"act_id"]),]
        df <- df[,-which.index(df,("rand"))]
        rownames(df) <- 1:nrow(df)
        df
}

quant.bin <- function(df,index = 1, split, title, class.number = 5,
                      zero.exclude = FALSE){
        #
        # Function for annotating a dataframe, df, with the bins for a set number of
        # quantiles. The number of bins is set by class.number. Each row is binned
        # according to the data in the column given by index. The new column of bin 
        # labels (1 lowest) is added after the column given by split. The new column
        # is given the name title. The option zero.exclude gived the option to ignore
        # zero values in calculating the quantiles (only for positive data). Thus all
        # zero values will default to bin 1
        #
        if(!require("dplyr")){
                install.packages("dplyr")
        }
        library(dplyr)
        #
        if(zero.exclude){
                mod.data <- df[which(df[,index]>0),]
        } else {
                mod.data <- df
        }
        #
        classes <- stats::quantile(mod.data[,index], 
                                   probs = seq(0, 1, length.out = (class.number + 1)))
        temp <- matrix(nrow = nrow(df),ncol = class.number)
        #
        for (i in 1:class.number){
                temp[,i] <- df[,index] <= rep(classes[(i+1)],nrow(df))       
        }
        length <- ncol(df)
        df <- df %>% 
                dplyr::mutate(
                        temp.class = rep(class.number+1, nrow(df)) - rowSums(temp)) %>%
                dplyr::select(1:split, length + 1, (split + 1):length)
        names(df)[split + 1] <- title
        df
}


select.narm <- function(df,index,select = TRUE){
        #
        # function to modify dplyr select by also removing all rows in the columns to be
        # selected where there is a missing value
        #
        if(!require("dplyr")){
                install.packages("dplyr")
        }
        library(dplyr)
        #
        for (i in index){
                df <- df %>%
                        dplyr::filter(is.na(df[,i]) == FALSE)
        }
        if(select){
                df <- dplyr::select(df,index)
        }
        df
}


quant.bin.plot <- function(df, index = c(1,2), class.number = 5, varname, binname, 
                           zero.exclude = FALSE, binwidth = NULL){
        #
        # function that plots a lattice of histograms for the data in df, in the 
        # column given by index[2]. The histograms are binned by the quantiles of
        # the column given by index[1]. The number of plots is given by class.number
        # if zero.exclude is true then zero values are ignored in calculating quantiles
        #
        if(!require("ggplot2")){
                install.packages("ggplot2")
        }
        library(ggplot2)        
        #
        df <- select.narm(df,index)
        #
        df <- quant.bin(df,index = 1, split = 2, title = "class_act",
                        class.number = class.number, zero.exclude = zero.exclude)
        df[,3] <- factor(df[,3])
        names(df) <- c("a", "b", "c")
        g <- ggplot(df, aes(b, col = c))
        g <- g + facet_grid(~ c)
        g <- g + stat_bin(binwidth = binwidth)
        g <- g + labs(title = "Quantile Binned Histograms", x = varname, col = binname)
        g
}



classifier.test <- function(df, index = c(1,2), top.class = TRUE, top.num = 5,
                            bot.class = TRUE, bot.num = 1){
        #
        # Function that can be used to evaluate a classifier. Given a data frame df,
        # the function compares the actual classes with those predicted (given as index
        # in this order). The result is the accuracy of the predictor across all 
        # classes. Setting top.class to True will evaluate the classifier on the top 
        # It will report the precision, recall and f measure. bot.class does the same
        # only for the bottom class (in both cases the class number for each must be
        # given specifically)
        #
        scores <- vector("numeric", length = 1)
        scores[1] <- mean(df[,index[1]]==df[,index[2]], na.rm = TRUE)
        names(scores) <- "accuracy"
        #
        if(top.class){
                t.pos <- df[which(df[,index[1]] == top.num & df[,index[2]]== top.num),]
                f.pos <- df[which(df[,index[1]] != top.num & df[,index[2]]== top.num),]
                f.neg <- df[which(df[,index[1]] == top.num & df[,index[2]]!= top.num),]
                top_recall <- nrow(t.pos)/(nrow(t.pos) + nrow(f.neg))
                top_precision <- nrow(t.pos)/(nrow(t.pos) + nrow(f.pos))
                f_top <- 2*top_recall*top_precision/(top_recall + top_precision)
                scores <- c(scores,top_recall,top_precision,f_top)
                tmp <- length(scores)
                names(scores)[(tmp-2):tmp] <- c("top_recall","top_precision","top_f")
        }
        if(bot.class){
                t.pos <- df[which(df[,index[1]] == bot.num & df[,index[2]]== bot.num),]
                f.pos <- df[which(df[,index[1]] != bot.num & df[,index[2]]== bot.num),]
                f.neg <- df[which(df[,index[1]] == bot.num & df[,index[2]]!= bot.num),]
                bot_recall <- nrow(t.pos)/(nrow(t.pos) + nrow(f.neg))
                bot_precision <- nrow(t.pos)/(nrow(t.pos) + nrow(f.pos))
                f_bot <- 2*bot_recall*bot_precision/(bot_recall + bot_precision)
                scores <- c(scores,bot_recall,bot_precision,f_bot)
                tmp <- length(scores)
                names(scores)[(tmp-2):tmp] <- c("bottom_recall","bottom_precision",
                                                "bottom_f")
        }
        scores
}

neg.bin.prob.add <- function(df,index,class.number,param,method = "cdf"){
        #
        # Function that takes a data frame and returns it appended with probabilities
        # for a given variable belonging to a set of negative binomial distributions,
        # given via fitted parameters param. The variable on which this is done is 
        # given by index. The number of distributions is given by class.number. The 
        # returned probabilities are given for each distribution in order, either the
        # probability density, or the cdf depending on the value of method.
        #
        prob.est <- matrix(nrow = nrow(df), ncol = class.number)
        #
        if(method == "density"){
                for(i in 1:class.number){
                        prob.est[,i] <- dnbinom(df[,index], size = param[[i]]$estimate[1],
                                                mu = param[[i]]$estimate[2])
                }
        }
        
        if(method == "cdf"){
                for(i in 1:class.number){
                        prob.est[,i] <- pnbinom(df[,index], size = param[[i]]$estimate[1],
                                                mu = param[[i]]$estimate[2])
                }
        }
        #
        df <- cbind(df,prob.est)
        df
}

neg.bin.prob <- function(df, index = c(1,2,3), class.number = 5,
                         method = "density", zero.exclude = FALSE){
        #
        # Function that splits a data frame based on the quantiles for the variable
        # given by the column as index[1], using only values in the training set as 
        # given by index[3]. For each group, a negative binomial distribution is fitted 
        # to the variable in index[2] via the MLE. The function returns the 
        # probabilities of each data point in both training and test sets for each class
        # This is either the pdf or cdf depending on the value of method
        #
        if(!require("dplyr")){
                install.packages("dplyr")
        }
        library(dplyr)
        #
        if(!require("fitdistrplus")){
                install.packages("fitdistrplus")
        }
        library(fitdistrplus)
        #
        df <- select.narm(df,index) %>% 
                quant.bin(index = 1, split = 3, title = "class_act",
                          class.number = class.number, zero.exclude = zero.exclude)
        names(df)[3] <- "set"
        #
        train.set <- dplyr::filter(df,set == "training")
        train.set <- quant.bin(train.set, index = 1, split = 4, title = "class_train",
                               class.number = class.number, zero.exclude = zero.exclude)
        #
        param <- NULL
        for (i in 1:class.number){
                param[[i]] <- dplyr::filter(train.set, class_train == i)[,2] %>%
                        fitdistrplus::fitdist("nbinom",method = "mle")
        }
        rm(train.set)
        #
        df <- neg.bin.prob.add(df,2,class.number,param,method)
        df
}


simple.quant.count <- function(df, index = c(1,2,3), classes = c(5,5,5),
                               title = c('A','B','C')){
        #
        # Function that returns a table of counts for the variables in the
        # indices given by index, seperated out by quantiles of these 
        # variables. The number of quantiles is given by classes for each
        # variable. The columns of the output are given the titles from title
        #
        if(!require("dplyr")){
                install.packages("dplyr")
        }
        library(dplyr)     
        
        tmp.l <- ncol(df)
        
        data <- quant.bin(df, index = index[1], split = tmp.l, 
                          title = 'class_1', class.number = classes[1]) 
        data <- quant.bin(data, index = index[2], split = tmp.l + 1, 
                          title = 'class_2', class.number = classes[2])
        data <- quant.bin(data, index = index[3], split = tmp.l + 2, 
                          title = 'class_3', class.number = classes[3])   
        data <- dplyr::group_by(data,class_1,class_2,class_3) %>%
                dplyr::summarise(count = n())
        names(data)[1:3] <- title
        data
}

cv.partition <-function(df,ratio_train,ratio_CV){
        #
        # Function that takes a data frame and appends it with a column indicating 
        # training, cross-validation and test sets. These are picked at random in
        # proportions given by ratio_train, ratio_CV (out of total of 1)
        #
        if(!require("dplyr")){
                install.packages("dplyr")
        }
        library(dplyr)
        #
        df <- df %>% 
                dplyr::mutate(act_id = 1:as.numeric(nrow(.)),
                              rand = runif(nrow(df)))
        cut_train <- ceiling(nrow(df)*ratio_train)
        cut_CV <- ceiling(nrow(df)*(ratio_train + ratio_CV))
        df <- df[order(df[,"rand"]),] %>%
                dplyr::mutate(
                        set = c(rep("training",times = cut_train),
                                rep("CV",times = cut_CV - cut_train),
                                rep("test",times = (nrow(df)-cut_CV)))
                )
        df <- df[order(df[,"act_id"]),]
        df <- df[,-which.index(df,("rand"))]
        rownames(df) <- 1:nrow(df)
        df
}

kern.density <- function(x,data){
        #
        # Function that takes a set of data, fits a gaussian kernal density 
        # function and then returns the density for a specified data point x
        #
        func <- stats::approxfun(stats::density(data))
        func(x)
}