
# require(arcdiagram)
# require(plotrix)

source("scripts/code_alleles/arcplot.r")
source("scripts/code_alleles/myarc.R")

print("<A5 script>")

  filename=paste0(MIOutFolder,tag,"MI.csv")
  MIPSexist<- file.exists(filename)

  if(MIPSexist){

        myData= read.table(filename, header=TRUE, sep=',', stringsAsFactors=FALSE)
        # filename=paste0("MIPS_output/",tag,"out.pdf")
        filename=paste0(MIOutFolder,tag,"_cutoff_",miPlotThreshold,"out.png")

        print(paste("<file out>", filename))
        print("<plotting results>")
        # pdf(filename,w=9,h=6,bg="white",title=tag)
        png(filename,w=900,h=600,bg="white",title=tag)

                # cutoff= quantile(myData$MInorm,probs=miPlotThreshold)
                cutoff= miPlotThreshold

                mmain=paste("arcs= MI>",cutoff)
                xxlab="Genome positions"
                yylab="Scaled mutual information index"

                xxlim<- c(min(myData$siteA,myData$siteB)-20, max(myData$siteA,myData$siteB)+20)

                yylim<- c(0,1)

                plot(myData$siteA, myData$MInorm, pch=21, col="black", xlab=xxlab, ylab=yylab, main=mmain, ylim=yylim, xlim=xxlim)

                points(myData$siteA, myData$MInorm, pch=21, col="black")
                points(myData$siteB, myData$MInorm, pch=21, col="black")
                points(myData$siteA, myData$MInorm, pch=20, col="skyblue")
                points(myData$siteB, myData$MInorm, pch=20, col="skyblue")

                # wantedIndexes<- which(myData$MInorm > cutoff)
                # if(length(wantedIndexes)==0){
                #     print("note, the cutoff for arcs is too high, no points selected")
                # }else{
                #     arcs= cbind(myData$siteA[wantedIndexes], myData$siteB[wantedIndexes])
                #     myArcs(edgelist=arcs, maxAbsHeight=max(yylim), col.nodes="grey11", bg.nodes= "grey66", pch=21, col.arcs = "orange")
                # }

                box(lwd=1.5)

        a=dev.off()
  }
