### Joerg Heintz, April 2017
### Health Care Engineering Systems Center 


# Funtion Library
Ingestion <- function(alpha, ColumnsOfInterest, fileLineOffset, mypath, lowlimit, uplimit){
        allfiles<-list.files(mypath)
        helper <- data.frame()
        mydata<-data.frame()
        number_of_loops<-(uplimit-lowlimit+1)/512
        for (i in allfiles){
                
                myfile<-paste0(mypath,"/", i)
                helper<-read.table(myfile, skip = fileLineOffset)[lowlimit:uplimit, ColumnsOfInterest]
                colnames(helper)<-c("x", "y_trace", "y_retrace")
                #####   introduce variables and format table
                
                #### extracting infos from filename
                fnc<-strsplit(i, "\\.")[[1]][1]
                fnc<-unlist(strsplit(fnc, " "))
                
                #### generating the data.frame 
                helper$f_on_samp<-fnc[2]
                helper$date<-as.Date(fnc[4], "%d %m %Y")
                helper$alpha<-alpha
                helper$ID<-1:length(helper$x) # ID will be overwritten, when points are tranformed into nodes.
                helper$tip<-fnc[3]
                helper$system<-fnc[1]
                #print(paste("file = ", i, "; #rows = ", length(helper$x) ))
                helper$y_trace<-helper$y_trace*alpha
                helper$y_retrace<-helper$y_retrace*alpha
                # Loop will be overwritten, with the next ingested file. 
                helper$loop<-rep(seq(1, number_of_loops), each = 512)
                # Calculating speed, retrieving data from AFT raw.txt file head   
                mycon<-file(myfile)
                open(mycon)
                lg<-read.table(mycon,skip=6,nrow=1, comment.char="?",check.names=FALSE)
                fr<-read.table(mycon,skip=1,nrow=1, comment.char="?",check.names=FALSE)
                close(mycon)
                helper$speed <- as.numeric(lg[[2]] * fr[[2]] * 2) # um/s
                helper$line_rate<-fr[[2]]
                helper$length<-lg[[2]]
                
                # collecting all the data one data frame
                mydata<-rbind(helper, mydata)
                #print(i)
        }
        
        # assignes unique loop numbers for all ingested files
        mydata$loop <- rep(seq(1, length(allfiles) * number_of_loops), each = 512)
        
        write.csv(mydata, "CellMembrane_Friction_20170324.csv", row.names = FALSE)
        mydata
}



#####################################################################################################
#####################################################################################################


PreprocessingIngestedRawData<-function(myfilename, y_trace__y_retrace = "y_trace"){
        setwd("/Users/joergheintz/Documents/06_Projects/FrictionCellMembranes")
        mypath<-""
        getwd()
        df<-read.csv(paste0(mypath, myfilename))[, c("x", y_trace__y_retrace,"loop", "f_on_samp", "system", "ID", "line_rate", "length", "tip", "alpha", "date")]
        #df<-read.xlsx("df_jump.xlsx", 2)
        mydata<-df
        # columnname loop changes from loop to half_loop, because only tracing or retracing is read. 
        colnames(mydata)<-c("x", "y", "half_loop", "f_on_samp", "system", "ID", "line_rate", "length", "tip", "alpha", "date")
        if (colnames(df)[2] == "y_trace") mydata$trace_d<-1 else  mydata$trace_d<-0
        mydata
}




#####################################################################################################
#####################################################################################################


#This function marks positive and negative slope and converts points into nodes, but 
# copying points to realize vertexes. 
StickSlip <- function(mydata, grp){
        colnames(mydata)[which(names(mydata) == grp)] <- "bygrp"
        mytemp<-mydata
        myData_mSL<-data.frame()
        for (i in unique(mytemp$bygrp)){
                #print(paste(grp, " = ", (i)))
                mydata<-mytemp[mytemp$bygrp == i, ]
                for (i in 1:length(mydata[mydata$bygrp == i,]$y)){
                        y1 <- mydata[i,]$y
                        y2 <- mydata[i+1,]$y
                        if (!is.na(y2) & (mydata$bygrp[i] == mydata$bygrp[i+1])) {
                                if (y1>y2) {mydata$mSL[i] = 0 ; a = 1; b = 0}
                                if (y1<y2) {mydata$mSL[i] = 1 ; a = 0; b = 1}
                        }
                }
                #for (i in 1:length(mydata$y)){
                k = length(mydata$y)
                i=1
                b=1
                while ((i < k) & !is.na(b)) {
                        a <- mydata[i,]$mSL
                        b <- mydata[i+1,]$mSL
                        if ( a==1 & !b == a) {
                                mydata[i+1, ]$mSL = 1
                                helper<-mydata[i+1, ]
                                helper$mSL<-0
                                mydata <- rbind(mydata, helper)
                                mydata<-mydata[order(mydata$ID),]
                                k = length(mydata$y)
                                i=i+1
                        }
                        if ( a==0 & !b == a) {
                                mydata[i+1, ]$mSL = 0
                                helper<-mydata[i+1, ]
                                helper$mSL<-1
                                mydata <- rbind(mydata, helper)
                                mydata<-mydata[order(mydata$ID),]
                                k = length(mydata$y)
                                i=i+1
                        }
                        i = i + 1
                        
                }
                myData_mSL<-rbind(myData_mSL, mydata )
        }
        myData_mSL
        ## correction on boarder between different applied normal forces f_on_samp
        i=0
        myerror<-data.frame()
        for (i in 2:length(myData_mSL$y)){
                if ((myData_mSL$bygrp[i] != myData_mSL$bygrp[i+1]) 
                    & (!is.na(myData_mSL$bygrp[i+1])) 
                    & (myData_mSL$y[i] == myData_mSL$y[i-1]) ){
                        myerror<-rbind(myData_mSL[i, ], myerror)
                        myData_mSL<-myData_mSL[-i, ]
                }
                myData_mSL$y[5]
        }
        mydata<-myData_mSL
        # marking sticks and slip 
        mydata$st_sl <- 1
        mydata[mydata$mSL < 1,"st_sl"]<- -1
colnames(mydata)[which(names(mydata) == "bygrp")] <- grp
mydata$ID <-1:length(mydata$x)
print("Stick - Slip function executed")
mydata
}


#####################################################################################################
#####################################################################################################


BuildCluster<-function(mydata, BinVariable){
        colnames(mydata)[which(names(mydata) == BinVariable)] <- "byBinVariable"
        # cluster
        n=0
        flipflop = 1
        for (i in 1:length(mydata$y)){
                if (!is.na(mydata$y[i+1]) & (mydata$half_loop[i] == mydata$half_loop[i+1])){
                        if (mydata$byBinVariable[i] == 1)  {
                                if (flipflop == 1) {n = n + 1}
                                mydata$cluster[i]<-n
                                flipflop=0
                        }
                        if (mydata$byBinVariable[i] == -1) {
                                if (flipflop == 0) {n = n + 1}
                                mydata$cluster[i]<-n
                                flipflop=1
                        }
                }
                else{
                        #catches the last [i+1] value and avoids NA. 
                        mydata[i,"cluster"] <- n
                        n=n+1
                        #print(i)
                        #print(paste("half-loop border = ", mydata$half_loop))
                }
        }
        colnames(mydata)[which(names(mydata) == "byBinVariable")] <- BinVariable

        mydata
}


#####################################################################################################
#####################################################################################################


MinMaxSlope <- function(mydata){
        # determining min, and max values by group and cluster
        mydata<-mydata %>% group_by(half_loop, cluster) %>% mutate(
                xmin = min(x),
                xmax = max(x),
                ymin = min(y),
                ymax = max(y),
                x.fitted.min = min(x.fitted),
                x.fitted.max = max(x.fitted),
                y.resid.min = min(y.resid),
                y.resid.max = max(y.resid)
                
        )
        # calculate dxmax, max, slope
        mydata <- mydata %>% mutate(
                dxmax = xmax - xmin,
                dymax = (ymax - ymin) * st_sl,
                slope = dymax/dxmax,
                dx.fitted.max = x.fitted.max - x.fitted.min,
                dy.resid.max = (y.resid.max - y.resid.min) * st_sl,
                slope.fit = dy.resid.max / dx.fitted.max
        )
        
        mydata <- mydata %>% select(x:y, dxmax, dymax, slope, reg.slope, reg.intercept,  xmin:ymax, x.fitted, y.resid, x.fitted.min, x.fitted.max, 
                                    y.resid.min, y.resid.max, dx.fitted.max, dy.resid.max, slope.fit,
                                    half_loop, cluster, st_sl, ID, f_on_samp, 
                                    system, trace_d, mSL, line_rate, length, tip, date, alpha)
        mydata[mydata$slope < 0, c("ymin", "ymax")] <- mydata[mydata$slope < 0, c("ymax", "ymin")]
        mydata <- ungroup(mydata)
        mydata
}


#####################################################################################################
#####################################################################################################


Segments <- function(mydata=mydata){
        mydata$dy<-0
        mydata$dx<-0
        for (i in 1:(length(mydata$cluster)-1)){
                mydata$dy[i+1]=abs((mydata$y[i]-mydata$y[i+1])) * mydata$st_sl[i+1]
                mydata$dx[i+1]=mydata$x[i]-mydata$x[i+1]
                if (mydata$half_loop[i+1] != mydata$half_loop[i]) {
                        mydata$dy[i+1]<-0
                        mydata$dx[i+1]<-0
                        
                }
        }
        mydata
}

#####################################################################################################
#####################################################################################################


SegmentsFit <- function(mydata=mydata){
        mydata$dy.resid<-0
        mydata$dx.fitted<-0
        for (i in 1:(length(mydata$cluster)-1)){
                mydata$dy.resid[i+1]=abs((mydata$y.resid[i]-mydata$y.resid[i+1])) * mydata$st_sl[i+1]
                mydata$dx.fitted[i+1]=mydata$x.fitted[i]-mydata$x.fitted[i+1]
                if (mydata$half_loop[i+1] != mydata$half_loop[i]) {
                        mydata$dy.resid[i+1]<-0
                        mydata$dx.fitted[i+1]<-0
                }
        }
        mydata
}



#####################################################################################################
#####################################################################################################


Smoothing <- function(mydata){
        mydata$ySmooth <- 0
        k <-mydata$y
        
        for (i in 3:(length(mydata$y)-1)) {
                m = c(k[i-2], k[i-1], k[i], k[i+1], k[i+2])
                k[i] <- mean(m)
        }
        
        mydata$ySmooth<-k
        mydata
        
}


#####################################################################################################
#####################################################################################################


InsertResidual<-function(mydata = mydata){
        mydata$x.fitted<-0
        mydata$y.resid<-0
        mydata$reg.slope<-0
        mydata$reg.intercept<-0
        for (s in unique(mydata$system)){
                mys<-mydata[mydata$system == s, ]
                for (l in unique(mys$half_loop)){
                        mysl <- mys[mys$half_loop == l, ]
                        
                        fit<-lm(y ~ x, data = mysl)
                        df <- augment(fit)
                       
                        mydata[mydata$half_loop == l,]$x.fitted <- as.numeric(df$.fitted)
                        mydata[mydata$half_loop == l,]$y.resid <- as.numeric(df$.resid)
                        mydata[mydata$half_loop == l,]$reg.slope <-  coef(fit)[2]
                        mydata[mydata$half_loop == l,]$reg.intercept<-  coef(fit)[1]
                }
        }
        mydata
}



#####################################################################################################
#####################################################################################################


PlotColFit<-function(mydata = mydata, myname = "sys.time",  p = 0, c= 3, w = 15, h = 18, xfitted=1, res = 200){
        plots<-list()
        n=0
        if (myname == "sys.time") myname <- paste0(Sys.time(),".png")
        
        for (s in unique(mydata$system)){
                mys<-mydata[mydata$system == s, ]
                f<-mys$f_on_samp
                v<-mys$line_rate
                d<-mys$date
                t<-mys$tip
                tr<-mys$trace_d
                for (l in unique(mys$half_loop)){
                        mysl <- mys[mys$half_loop == l, ]
                        if (xfitted==1){
                                g1 <- ggplot(data =mysl, aes(x = x.fitted, y = y.resid)) + 
                                ggtitle(paste0(s, ", ",f, ", ", v, " Hz", ", trace/re-trace: ", tr, ", loop: ", l, ", ", t, ", ", d))
                                g1 <- g1 + geom_point(size = 0.05)
                                g1 <- g1 + geom_line(size = 0.05, colour = "black")
                        }
                        if (xfitted == 0) {
                                g1 <- ggplot(data = mysl, aes(x = x, y = y.resid)) +
                                        ggtitle(paste0(s, ", ",f, ", ", v, " Hz", ", trace/re-trace: ", tr, ", loop: ", l, ", ", t, ", ", d))
                                g1 <- g1 + geom_point(size = 0.05)
                                g1 <- g1 + geom_line(size = 0.05, colour = "black")
                        }
                                
                        n=n+1
                        plots[[n]] <- g1
                }
        }
        if (p == 1) png(myname, units="in", width=w, height=h, res=res)
        multiplot(plotlist = plots, cols = c)
}


#####################################################################################################
#####################################################################################################


PlotCol_RegLine<-function(mydata = mydata, myname = "sys.time", p = 0,  c= 3, w = 15, h = 18, res = 200){
        plots<-list()
        n=0
        if (myname == "sys.time") myname <- paste0(Sys.time(),".png")
        
        for (s in unique(mydata$system)){
                mys<-mydata[mydata$system == s, ]
                f<-mys$f_on_samp
                v<-mys$line_rate
                d<-mys$date
                t<-mys$tip
                tr<-mys$trace_d
                for (l in unique(mys$half_loop)){
                        mysl <- mys[mys$half_loop == l, ]
                        
                        # curve + regression line
                        g1 <- ggplot(mysl, aes(x=x,y=y, group = half_loop))
                        g1 = g1 + geom_line(size = 0.05, colour = "black")
                        g1 = g1 + geom_point(size = 0.05, colour = "lightblue")
                        g1 = g1 + geom_smooth(method='lm',formula=y~x)
                        
                        # Segments and peaks, y-dy, dy = y+1 - y
                        PeakSegm<-mysl
                        PeakSegm<-PeakSegm[order(-abs(PeakSegm$dy)), ]
                        PeakSegm<-PeakSegm[1:5, ]
                        g1 = g1 + geom_segment(data = PeakSegm, aes(x = x, y = y, xend = x+dx, yend = y-dy), size = 0.8, colour = "red")
                        g1 = g1 + geom_segment(data = PeakSegm, aes(x = xmin, y = ymin, xend = xmax, yend = ymax), size = 0.8, linetype = 2, colour = "blue")
                        g1 = g1 + ggtitle(paste0(s, ", ",f, ", ", v, " Hz", ", trace/re-trace: ", tr, ", loop: ", l, ", ", t, ", ", d))
                        n=n+1
                        plots[[n]] <- g1
                }
        }
        if (p == 1) png(myname, units="in", width=w, height=h, res=res)
        multiplot(plotlist = plots, cols = c)
}



#####################################################################################################
#####################################################################################################

PlotCol_RegLine_Fit<-function(mydata = mydata, myname = "sys.time", p = 0,  c= 3, w = 15, h = 18, res = 200){
        plots<-list()
        n=0
        if (myname == "sys.time") myname <- paste0(Sys.time(),".png")
        
        for (s in unique(mydata$system)){
                mys<-mydata[mydata$system == s, ]
                f<-mys$f_on_samp
                v<-mys$line_rate
                d<-mys$date
                t<-mys$tip
                tr<-mys$trace_d
                for (l in unique(mys$half_loop)){
                        mysl <- mys[mys$half_loop == l, ]
                        
                        # curve + regression line
                        g1 <- ggplot(mysl, aes(x=x.fitted,y=y.resid, group = half_loop))
                        g1 = g1 + geom_line(size = 0.05, colour = "black")
                        g1 = g1 + geom_point(size = 0.05, colour = "lightblue")
                        #g1 = g1 + geom_smooth(method='lm',formula=y.resid~x.fitted)
                        
                        # Segments and peaks, y-dy, dy = y+1 - y
                        PeakSegm<-mysl
                        PeakSegm<-PeakSegm[order(-abs(PeakSegm$dy.resid)), ]
                        PeakSegm<-PeakSegm[1:5, ]
                        g1 = g1 + geom_segment(data = PeakSegm, aes(x = x.fitted, y = y.resid, xend = x.fitted+dx.fitted, yend = y.resid-dy.resid), size = 0.8, colour = "red")
                        g1 = g1 + geom_segment(data = PeakSegm, aes(x = x.fitted.min, y = y.resid.min, xend = x.fitted.max, yend = y.resid.max), size = 0.8, linetype = 2, colour = "blue")
                        g1 = g1 + ggtitle(paste0(s, ", ",f, ", ", v, " Hz", ", trace/re-trace: ", tr, ", loop: ", l, ", ", t, ", ", d))
                        n=n+1
                        plots[[n]] <- g1
                }
        }
        if (p == 1) png(myname, units="in", width=w, height=h, res=res)
        multiplot(plotlist = plots, cols = c)
}


#####################################################################################################
                                #   Boxplots   #
#####################################################################################################

PlotBoxSystemForce_AllinOne<-function(mydata = mydata, myname = "sys.time", p=0,  c= 3, w = 15, h = 18, res = 200, colS = colSet){
        if (myname == "sys.time") myname <- paste0(Sys.time(),".png")
        g1 <- ggplot(mydata, aes(system, reg.slope))
        g1 <- g1 + geom_boxplot(aes(fill = system), outlier.colour = "red") + geom_point(size = 1, aes(colour = f_on_samp))
        g1 <- g1 + ggtitle(paste0("Regression Slopes"))
        g1 <- g1 + geom_abline(slope = 0, intercept = 0, linetype = 2, colour = 'red')
        
        if (p==1) png(myname, units="in", width=w, height=h, res=res)
        g1 <- g1 + scale_fill_brewer(palette = colS)
        g1
}


PlotBoxSystemForce<-function(mydata = mydata, myname = "sys.time", p=0,  c= 3, w = 15, h = 18, res = 200, colScale = colSc){
        plots<-list()
        n=0
        if (myname == "sys.time") myname <- paste0(Sys.time(),".png")
        
        for (s in unique(mydata$system)){
                mys<-mydata[mydata$system == s, ]
                v<-mys$line_rate
                d<-mys$date
                t<-mys$tip
                tr<-mys$trace_d
                for (l in unique(mys$f_on_samp)){
                        mysl <- mys[mys$f_on_samp == l, ]
                        
                        g1 <- ggplot(mysl, aes(system, reg.slope, colour = system))
                        g1 <- g1 + ggtitle(paste0("Regression Slopes", ", ",l, ", ", v, " Hz", ", trace/re-trace: ", tr,  ", ", t, ", ", d))
                        g1 <- g1 + geom_boxplot(aes(), outlier.colour = "red") + geom_point(size = 1, colour = "red")
                        n=n+1
                        plots[[n]] <- g1 + colScale
                }
        }
        if (p==1) png(myname, units="in", width=w, height=h, res=res)
        multiplot(plotlist = plots, cols =3)
}


#####################################################################################################
#####################################################################################################


PlotColLm<-function(mydata = mydata, myname = "sys.time", p=0,  c= 3, w = 15, h = 18, res = 200){
        plots<-list()
        n=0
        if (myname == "sys.time") myname <- paste0(Sys.time(),".png")
        
        for (s in unique(mydata$system)){
                mys<-mydata[mydata$system == s, ]
                f<-mys$f_on_samp
                v<-mys$line_rate
                d<-mys$date
                t<-mys$tip
                tr<-mys$trace_d
                for (l in unique(mys$half_loop)){
                        mysl <- mys[mys$half_loop == l, ]
                        
                        fit<-lm(y ~ x, data = mysl)
                        df <- augment(fit)
                        
                        g1 <- ggplot(df, aes(x = .fitted, y = .resid)) + #xlim(0.2, 1.2)
                                ggtitle(paste0(s, ", ",f, ", ", v, " Hz", ", trace/re-trace: ", tr, ", loop: ", l, ", ", t, ", ", d))
                        g1 <- g1 + geom_point(size = 0.05)
                        g1 <- g1 + geom_line(size = 0.05, colour = "black")
                        
                        n=n+1
                        plots[[n]] <- g1
                }
        }
        if (p==1) png(myname, units="in", width=w, height=h, res=res)
        multiplot(plotlist = plots, cols =c)
}


#####################################################################################################
                        # Slope Plots #
#####################################################################################################

SlopeHistogram<-function(mydata = mydata, low = 0, up = 0, set_x_marker = 0, logscale = FALSE, bi = 500, pr = 0, w =5, h = 3, res = 200, myname = "sys.time", colS = colSet){
        if (low == 0 & up == 0 ) {
                low<-min(mydata$slope)
                up<-max(mydata$slope)}
        h1 <- ggplot()
        h1 <- h1 + geom_vline(xintercept = set_x_marker)
        if (logscale == FALSE){
                h1 <- h1 + geom_histogram(data = mydata[(mydata$slope>low & mydata$slope<up), ], aes(slope, fill = system), bins = bi)
        }
        if (logscale == TRUE){
                min_plus <- as.logical(low<0 & up>=0)
                min_min <- as.logical(low<0 & up<0)
                plus_plus <- as.logical(low>0 & up>0)
                
                if(min_min == TRUE) 
                        h1 <- h1 + geom_histogram(data = mydata[(mydata$slope>low & mydata$slope<up), ], aes(-log(abs(slope)), fill = system), bins = bi)
                if(min_plus == TRUE)
                        h1 <- h1 + geom_histogram(data = mydata[(mydata$slope>low & mydata$slope<0), ], aes(-log(abs(slope)), fill = system), bins = bi)
                        h1 <- h1 + geom_histogram(data = mydata[(mydata$slope>0 & mydata$slope<up), ], aes(log(slope), fill = system), bins = bi)
                if(plus_plus == TRUE)
                        h1 <- h1 + geom_histogram(data = mydata[(mydata$slope>low & mydata$slope<up), ], aes(log(slope), fill = system), bins = bi)
                
        }
        if (myname == "sys.time") myname <- paste0(Sys.time(),".png")
        if (pr == 1) png(myname, units="in", width=w, height=h, res=res)
        h1 <- h1 + scale_fill_brewer(palette = colS)
        h1
}

#####################################################################################################
#####################################################################################################

PlotSlope_1<-function(mydata = mydata, logscale = FALSE,  myname = "sys.time", p = 0, c= 3, w = 10, h = 6, res = 200, colScale = colSc){
        plots<-list()
        n=0
        if (myname == "sys.time") myname <- paste0(Sys.time(),".png")
        for (fsamp in unique(mydata$f_on_samp)){
                mys<-mydata[mydata$f_on_samp == fsamp, ]
                v<-mys$line_rate
                d<-mys$date
                t<-mys$tip
                tr<-mys$trace_d
                for (s in unique(mys$system)){
                        mysl <- mys[mys$system == s, ]
                        g1 <- ggplot()
                        if (logscale == TRUE){
                                g1 <- g1 + geom_point(data = mysl[mysl$slope<0,], aes(x=-log(abs(dy)) ,y=-log(abs(slope)), colour = system), size = 0.5)
                                g1 <- g1 + geom_point(data = mysl[mysl$slope>0,], aes(x=log(dy) ,y=log(slope), colour = system), size = 0.5)
                                g1 <- g1 + colScale
                        }
                        if (logscale == FALSE){
                                g1 <- g1 + geom_point(data = mysl, aes(x=dy,y=slope, colour = system), size = 0.5)
                                g1 <- g1 + colScale
                        }
                        g1 <- g1 + ggtitle(paste0(s, ", ",fsamp, ", ", v, " Hz", ", trace/re-trace: ", tr, ", ", t, ", ", d))
                        n=n+1
                        plots[[n]] <- g1
                }
        }
        if (p == 1) png(myname, units="in", width=w, height=h, res=res)
        multiplot(plotlist = plots, cols =c)
}


PlotSlope_2<-function(mydata = mydata, logscale = FALSE,  myname = "sys.time", p = 0, c= 2, w = 15, h = 18, res = 200, colScale = colSc){
        plots<-list()
        n=0
        if (myname == "sys.time") myname <- paste0(Sys.time(),".png")
        for (fsamp in unique(mydata$f_on_samp)){
                mys<-mydata[mydata$f_on_samp == fsamp, ]
                v<-mys$line_rate
                d<-mys$date
                t<-mys$tip
                tr<-mys$trace_d
                for (s in unique(mys$system)){
                        mysl <- mys[mys$system == s, ]
                        g1 <- ggplot()
                        if (logscale == TRUE){
                                g1 <- g1 + geom_point(data = mysl[mysl$dy<0,], aes(x=-log(abs(dy)),y=slope, colour = system), size = 0.5)
                                g1 <- g1 + geom_point(data = mysl[mysl$dy>0,], aes(x= log(dy) ,y=slope, colour = system), size = 0.5)
                                g1 <- g1 + colScale
                        }
                        if (logscale == FALSE){
                                g1 <- g1 + geom_point(data = mysl, aes(x=dy,y=slope, colour = system), size = 0.5)
                                g1 <- g1 + colScale
                        }
                        g1 <- g1 + ggtitle(paste0(s, ", ",fsamp, ", ", v, " Hz", ", trace/re-trace: ", tr, ", ", t, ", ", d))
                        n=n+1
                        plots[[n]] <- g1
                }
        }
        if (p == 1) png(myname, units="in", width=w, height=h, res=res)
        multiplot(plotlist = plots, cols =c)
}

#####################################################################################################
#####################################################################################################

PlotSlope<-function(mydata = mydata, myname = "sys.time",  p = 0, c= 3, w = 15, h = 18, res = 200){
        plots<-list()
        n=0
        if (myname == "sys.time") myname <- paste0(Sys.time(),".png")
        
        for (s in unique(mydata$system)){
                mys<-mydata[mydata$system == s, ]
                f<-mys$f_on_samp
                v<-mys$line_rate
                d<-mys$date
                t<-mys$tip
                tr<-mys$trace_d
                for (l in unique(mys$half_loop)){
                        mysl <- mys[mys$half_loop == l, ]
                        
                        g1 <- ggplot(mysl, aes(x=dy,y=slope, group = half_loop, colour = system))
                        g1 <- g1 + geom_point(size = 0.2)
                        g1 <- g1 + ggtitle(paste0(s, ", ",f, ", ", v, " Hz", ", trace/re-trace: ", tr, ", loop: ", l, ", ", t, ", ", d))
                        
                        n=n+1
                        plots[[n]] <- g1
                }
        }
        if (p == 1) png(myname, units="in", width=w, height=h, res=res)
        multiplot(plotlist = plots, cols =c)
}

#####################################################################################################
#####################################################################################################


dyHistogram<-function(mydata, low = 0, up = 0, c = 2, set_x_marker = 0,  same_scale = FALSE, logscale = FALSE, bi = 500, pr = 0, w =5, h = 3, res = 200, myname = "sys.time",  colS = colSet){
        plots<-list()
        n=0
        mindy = 0
        maxdy = 0
        mincount = 0
        maxcount = 0
        if (low == 0 & up == 0 ) {
                low<-min(mydata$slope)
                up<-max(mydata$slope)}
        
        for (i in unique(mydata$f_on_samp)){
                myfdata<-mydata[mydata$f_on_samp == i & mydata$dy!=0,]
                h1 <- ggplot()
                h1 <- h1 + ggtitle(paste0("Load = ",i))
                #h1 <- h1 + geom_vline(xintercept = set_x_marker)
                if (logscale == FALSE){
                        h1 <- h1 + geom_histogram(data = myfdata[(myfdata$slope>low & myfdata$slope<up), ], aes(dy, fill = system), bins = bi)
                }
                if (logscale == TRUE){
                        min_plus <- as.logical(low<0 & up>=0)
                        min_min <- as.logical(low<0 & up<=0)
                        plus_plus <- as.logical(low>0 & up>=0)
                        
                        if(min_min == TRUE) 
                                h1 <- h1 + geom_histogram(data = myfdata[(myfdata$slope>low & myfdata$slope<up), ], aes(-log(abs(dy)), colour = system), bins = bi)
                        if(min_plus == TRUE)
                                h1 <- h1 + geom_histogram(data = myfdata[(myfdata$slope>low & myfdata$slope<=0), ], aes(-log(abs(dy)), colour = system), bins = bi)
                                h1 <- h1 + geom_histogram(data = myfdata[(myfdata$slope>=0 & myfdata$slope<up), ], aes(log(dy), colour = system), bins = bi)
                        if(plus_plus == TRUE)
                                h1 <- h1 + geom_histogram(data = myfdata[(myfdata$slope>low & myfdata$slope<up), ], aes(log(dy), colour = system), bins = bi)
                }
                # scales 
                if (same_scale == TRUE){
                        min_ct<-as.numeric(min(hist(myfdata[myfdata$slope>low & myfdata$slope<up, ]$dy, plot = FALSE, breaks = bi)$counts))
                        max_ct<-as.numeric(max(hist(myfdata[myfdata$slope>low & myfdata$slope<up, ]$dy, plot = FALSE, breaks = bi)$counts))
                        if (mindy > min(myfdata[myfdata$slope>low & myfdata$slope<up, ]$dy)) mindy <- min(myfdata$dy)
                        if (maxdy < max(myfdata[myfdata$slope>low & myfdata$slope<up, ]$dy)) maxdy <- max(myfdata$dy)
                        if (mincount > min_ct) mincount <- min_ct
                        if (maxcount < max_ct) maxcount <- max_ct
                        #print(paste("min_ct = ", min_ct, "max_ct = ", max_ct))
                        #print(paste("counts histogram = ", hist(myfdata[myfdata$slope>low & myfdata$slope<up & myfdata$dy != 0, ]$dy, breaks = bi)$counts))
                }
                n=n+1
                plots[[n]] <- h1 + scale_fill_brewer(palette = colS)
        }
        
        if (same_scale == TRUE){
                for (i in 1:n){
                        plots[[i]] <- plots[[i]] + ylim(c(mincount, maxcount)) + xlim(c(mindy, maxdy))
                }
        }
        
        if (myname == "sys.time") myname <- paste0(Sys.time(),".png")
        if (pr == 1) png(myname, units="in", width=w, height=h, res=res)
        multiplot(plotlist = plots, cols = c)
}

#####################################################################################################
#####################################################################################################


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
        library(grid)
        
        # Make a list from the ... arguments and plotlist
        plots <- c(list(...), plotlist)
        
        numPlots = length(plots)
        
        # If layout is NULL, then use 'cols' to determine layout
        if (is.null(layout)) {
                # Make the panel
                # ncol: Number of columns of plots
                # nrow: Number of rows needed, calculated from # of cols
                layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                                 ncol = cols, nrow = ceiling(numPlots/cols))
        }
        
        if (numPlots==1) {
                print(plots[[1]])
                
        } else {
                # Set up the page
                grid.newpage()
                pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
                
                # Make each plot, in the correct location
                for (i in 1:numPlots) {
                        # Get the i,j matrix positions of the regions that contain this subplot
                        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
                        
                        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                                        layout.pos.col = matchidx$col))
                }
        }
}
