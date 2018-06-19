##Install the packages below if you have not already
##install.packages(c("ggplot2","plyr","shinythemes"))

library(ggplot2)
library(plyr)
library(shinythemes)

##User must specify full path for the directory that contains Winnow Output Test Files
##Change the directory in quotes below to the location where your files are stored
##This is the only step required to run the app
##⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇︎⬇

dir="~/Stapleton_Lab/Projects/Visualization/Demonstrate/Winnow Results"

##⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎⬆︎

setwd(dir)
##Creates a list of all files in the folder with a type of "txt"
files <- Sys.glob("*.txt")
##Reads all of the observations in to a list of each of the files
lstF <- lapply(files, function(x) read.table(x, header=TRUE))
##Merges all the files into one large dataset
dt <- rbindlist(lstF, fill=TRUE)
##Adds a column for the filenames
dt[,filenm:=rep(files, sapply(lstF,nrow))]
##Creates a new data set with just True and False Positives with filenames and number of recurring for each file
TPFP <- ddply(dt, .(tp, fp, filenm), summarise, count=length(filenm))

##Determines maximum value of all true positives
tpmax<-function(list){
    tps<-list()
    for (i in 1:length(list)){
        tps[[i]]<-list[[i]]$tp
    }
    y<-unlist(lapply(tps, max))
    return(max(y))
}
ta<-tpmax(lstF)
#Determines minimum value of all true positives
tpmin<-function(list){
    tps<-list()
    for (i in 1:length(list)){
        tps[[i]]<-list[[i]]$tp
    }
    y<-unlist(lapply(tps, min))
    return(min(y))
}
tb<-tpmin(lstF)
#Determines median value of all true positives
tpmed<-function(list){
    tps<-list()
    for (i in 1:length(list)){
        tps[[i]]<-list[[i]]$tp
    }
    y<-unlist(lapply(tps, median))
    return(median(y))
}
tc<-tpmed(lstF)
#Determines maximum value of all false positives
fpmax<-function(list){
    fps<-list()
    for (i in 1:length(list)){
        fps[[i]]<-list[[i]]$fp
    }
    y<-unlist(lapply(fps, max))
    return(max(y))
}
fa<-fpmax(lstF)
#Determines minimum value of all false positives
fpmin<-function(list){
    fps<-list()
    for (i in 1:length(list)){
        fps[[i]]<-list[[i]]$fp
    }
    y<-unlist(lapply(fps, min))
    return(min(y))
}
fb<-fpmin(lstF)
#Determines median value of all false positives
fpmed<-function(list){
    fps<-list()
    for (i in 1:length(list)){
        fps[[i]]<-list[[i]]$fp
        y<-unlist(lapply(fps, median))
        return(median(y))
    }
}
fc<-fpmed(lstF)


shinyApp(
  
  ui = fluidPage(
      ##there are several themes that can be used from the shinytheme package if you want a different look
      theme = shinytheme("superhero"),

      titlePanel(
          fluidRow(
              column(3, strong("Winnow Output")), 
              column(2, img(height = 150, width = 150, src = "spirit_black.jpg")),
              column(2, img(height = 150, width = 150, src = "nsf.jpg")),
              column(5, img(height = 150, width = 450, src = "cyverse.jpg"))
          )
      ),
      
    sidebarLayout(
      sidebarPanel(
        fluidRow(
          column(width=12,
                 div(class = "option_group",
                     radioButtons("plot_type", "Plot type",
                     c("Scatter", "Linear Rates", "AUC by MAE","True-Positives by False-Positives"), inline = FALSE),
                     
                     conditionalPanel("input.plot_type === 'Scatter'",
                                      selectInput("pvar1", "x-Var",
                                                  c("tp" = "tp",
                                                    "fp" = "fp"
                                                    ),
                                                  selected = "tp"
                                      ),
                                      selectInput("pvar2", "y-Var",
                                                  c("tp" = "tp",
                                                    "fp" = "fp"
                                                  ),
                                                  selected = "fp"
                                      )
                     ),
                     conditionalPanel("input.plot_type === 'Linear Rates'",
                                      selectInput("lvar1", "x-Var",
                                                  c("tpr" = "tpr",
                                                    "fpr" = "fpr"
                                                    ),
                                                  selected = "tpr"
                                                  ),
                                      selectInput("lvar2", "y-Var",
                                                  c("tpr" = "tpr",
                                                    "fpr" = "fpr"
                                                    ),
                                                  selected = "fpr"
                                                  )
                                      ),
                     conditionalPanel("input.plot_type === 'AUC by MAE'",
                                      sliderInput("auc.min", "AUC axis minimum", min = 0, max = 2, value = .5,step=0.05),
                                      sliderInput("auc.max", "AUC axis maximum", min = 0, max = 2, value = 1,step=0.05),
                                      sliderInput("mae.min", "MAE axis minimum", min = 0, max = 2, value = 0,step=0.05),
                                      sliderInput("mae.max", "MAE axis maximum", min = 0, max = 2, value = 1,step=0.05)
                     )
                  )#ends div
                 )#ends column
        )#ends fluidRow
      )#ends sidebarPanel
          ,
      mainPanel(
        tabsetPanel(
          tabPanel("Plot",plotOutput("plot1")),
          tabPanel("Table",
                   
                   fluidRow(
                     column(4,
                            h4("means"),
                            tableOutput("meanresults")
                     ),
                     column(4,
                            h4("sums"),
                            tableOutput("sumresults")
                     ),
                     column(4,
                            h4("means"),
                            tableOutput("meanresults2")
                     )
                   )#end fluidRow
          )#ends tablPale
        )#ends tabsetPanel  
      )#ends mainPanel
     )#ends sidebarLayout
    ),#ends ui = fluidPage

  
  
  server = function(input, output) {
    
   
      
    output$plot1 = renderPlot({
      if (input$plot_type == "Scatter") {
        ggplot(dt,aes_string(x=input$pvar1,y=input$pvar2))+geom_point(color="firebrick")
      } else if (input$plot_type == "Linear Rates") {
        ggplot(dt,aes_string(x=input$lvar1,y=input$lvar2))+geom_line(size = 1, alpha = 1 )+
          labs(title= "Comparison of Rates")
      }
        else if (input$plot_type == "AUC by MAE"){
            plot(lstF[[1]]$mae, lstF[[1]]$auc, main="Plot of AUC by MAE", xlab="Mean Absolute Error (MAE)", ylab="Area under R-O Curve (AUC)", 
                 pch=21, bg="black", xlim=c(input$mae.min, input$mae.max), ylim=c(input$auc.min,input$auc.max))
            plotcol<-c("black")
           
             if (length(lstF) > 1){
                #Create overlapping data plots to compare potentially by GWAS tool
                #assuming that the length of the Winnow files is at least 2
                for (i in 2:length(lstF)){
                    points(lstF[[i]]$mae, lstF[[i]]$auc, main="Plot of AUC by MAE", xlab="Mean Absolute Error (MAE)", ylab="Area under R-O Curve (AUC)",
                           pch=21, bg=rainbow(i+1)[i], xlim=c(input$mae.min, input$mae.max), ylim=c(input$auc.min, input$auc.max))
                    plotcol[i]<-rainbow(i+1)[i]
                    
                }
             }
        }
        else if (input$plot_type == "True-Positives by False-Positives"){
        p <- ggplot(TPFP, aes(x=fp, y=tp),environment=environment())
   
        p2 <- p +
            geom_rect(data=dt[1,], aes(xmin=fc, xmax=fa, ymin=tc, ymax=ta),
                      alpha=0.2, fill="blue", linetype=0) +
            geom_rect(data=dt[1,], aes(xmin=fb, xmax=fc, ymin=tc, ymax=ta),
                      alpha=0.2,fill="green", linetype=0) +
            geom_rect(data=dt[1,], aes(xmin=fb, xmax=fc, ymin=tb, ymax=tc),
                      alpha=0.2, fill="blue", linetype=0) +
            geom_rect(data=dt[1,], aes(xmin=fc, xmax=fa, ymin=tb, ymax=tc),
                      alpha=0.2, fill="gray", linetype=0) +
            theme(panel.background=element_rect(fill='white', colour='black')) +
            theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
            geom_point(aes(colour=filenm, size=count)) +
            scale_size_continuous(range=c(2,8)) +
            xlab("False Positives") +
            ylab("True Positives") +
            ggtitle("False Positives by True Positves") +
            xlim(0, fa) + ylim(0, ta)
        print(p2)
        }
        
    })
    
    output$meanresults = renderTable({apply(dt[,1:5],2,mean)},rownames=TRUE)
    
    output$sumresults = renderTable(apply(dt[,6:9],2,sum),rownames=TRUE)
    
    output$meanresults2 = renderTable(apply(dt[,10:16],2,mean),rownames=TRUE)
    
    
    
  }#ends server
  
  
  )#ends shinyApp


