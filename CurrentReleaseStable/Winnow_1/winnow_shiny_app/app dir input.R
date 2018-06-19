#install.packages(c("shinythemes", "shinyFiles"))
library(ggplot2)
library(Cairo)  
library(DT)
library(plyr)
library(shinythemes)
library(shinyFiles)
#source("var.R")
##shinyFilesExample()
##User must specify full path for the directory that contains Winnow Output Test Files
# dir="~/Stapleton_Lab/Projects/Visualization/Demonstrate/Winnow Results"



shinyApp(
  
  ui = fluidPage(
      
      theme = shinytheme("superhero"),

      titlePanel(
          fluidRow(
              column(9, strong("Winnow Output")), 
              column(3, img(height = 150, width = 150, src = "spirit_black.jpg"))
          )
      ),
      
    sidebarLayout(
      sidebarPanel(
        fluidRow(
          column(width=12,
                 shinyDirButton('directory', 'Folder select', 'Please select a folder'),
                 h3("Directory Selected"),
                 verbatimTextOutput("dirpath"),
                 
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
      volumes <- getVolumes()   
      folderInput1 <- reactive({
          shinyDirChoose(input, 'directory', roots = volumes, session = session, 
                         restrictions = system.file(package = 'base'))
          return(parseDirPath(volumes, input$directory))
      })
      
      output$directorypath = renderPrint({  folderInput1()  })
      
      files1 <- reactive({
          list.files(path = folderInput1(), pattern = "*.txt", full.names = T)
      })
      
      nFiles1 <- reactive({ length(files1() ) })
      
     
      # readFiles <- function(dir) {
      #     setwd(input$dir)
      #     files <- Sys.glob("*.txt")
      #     listOfFiles <- lapply(files, function(x) read.table(x, header=TRUE))
      #     return(listOfFiles)
      # }
      # myfiles <-readFiles(input$dir)
      filenames <- reactive ({file_path_sans_ext(files1)})
      
      myfiles <- reactive({lapply(files1, function(x) read.table(x, header=TRUE))})
      
      tpmax<-function(list){
          tps<-list()
          for (i in 1:length(list)){
              tps[[i]]<-list[[i]]$tp
          }
          y<-unlist(lapply(tps, max))
          return(max(y))
      }
      ta<-reactive({tpmax(myfiles)})
      #Determines minimum value of all true positives
      tpmin<-function(list){
          tps<-list()
          for (i in 1:length(list)){
              tps[[i]]<-list[[i]]$tp
          }
          y<-unlist(lapply(tps, min))
          return(min(y))
      }
      tb<-reactive({tpmin(myfiles)})
      #Determines median value of all true positives
      tpmed<-function(list){
          tps<-list()
          for (i in 1:length(list)){
              tps[[i]]<-list[[i]]$tp
          }
          y<-unlist(lapply(tps, median))
          return(median(y))
      }
      tc<-reactive({tpmed(myfiles)})
      #Determines maximum value of all false positives
      fpmax<-function(list){
          fps<-list()
          for (i in 1:length(list)){
              fps[[i]]<-list[[i]]$fp
          }
          y<-unlist(lapply(fps, max))
          return(max(y))
      }
      fa<-reactive ({fpmax(myfiles)})
      #Determines minimum value of all false positives
      fpmin<-function(list){
          fps<-list()
          for (i in 1:length(list)){
              fps[[i]]<-list[[i]]$fp
          }
          y<-unlist(lapply(fps, min))
          return(min(y))
      }
      fb<-reactive({fpmin(myfiles)})
      #Determines median value of all false positives
      fpmed<-function(list){
          fps<-list()
          for (i in 1:length(list)){
              fps[[i]]<-list[[i]]$fp
              y<-unlist(lapply(fps, median))
              return(median(y))
          }
      }
      fc<-reactive({fpmed(myfiles)})
      
      


      all.data<<-reactive ({do.call("rbind", myfiles)})

      TPFP <- reactive ({ddply(all.data, .(tp, fp, filenames), summarize, count=length(filenames))})

      
    output$plot1 = renderPlot({
      if (input$plot_type == "Scatter") {
        ggplot(all.data(),aes_string(x=input$pvar1,y=input$pvar2))+geom_point(color="firebrick")
      } else if (input$plot_type == "Linear Rates") {
        ggplot(all.data,aes_string(x=input$lvar1,y=input$lvar2))+geom_line(size = 1, alpha = 1 )+
          labs(title= "Comparison of Rates")
      }
        else if (input$plot_type == "AUC by MAE"){
            plot(myfiles[[1]]$mae, myfiles[[1]]$auc, main="Plot of AUC by MAE", xlab="Mean Absolute Error (MAE)", ylab="Area under R-O Curve (AUC)", 
                 pch=21, bg="black", xlim=c(input$mae.min, input$mae.max), ylim=c(input$auc.min,input$auc.max))
            plotcol<-c("black")
           
             if (length(myfiles) > 1){
                #Create overlapping data plots to compare potentially by GWAS tool
                #assuming that the length of the Winnow files is at least 2
                for (i in 2:length(myfiles)){
                    points(myfiles[[i]]$mae, myfiles[[i]]$auc, main="Plot of AUC by MAE", xlab="Mean Absolute Error (MAE)", ylab="Area under R-O Curve (AUC)",
                           pch=21, bg=rainbow(i+1)[i], xlim=c(input$mae.min, input$mae.max), ylim=c(input$auc.min, input$auc.max))
                    plotcol[i]<-rainbow(i+1)[i]
                    
                }
             }
        }
        else if (input$plot_type == "True-Positives by False-Positives"){
        p <- ggplot(TPFP(), aes_string(x=fp, y=tp),environment=environment())
   
        p2 <- p +
            geom_rect(data=all.data[1,], aes(xmin=fc, xmax=fa, ymin=tc, ymax=ta),
                      alpha=0.2, fill="blue", linetype=0) +
            geom_rect(data=all.data[1,], aes(xmin=fb, xmax=fc, ymin=tc, ymax=ta),
                      alpha=0.2,fill="green", linetype=0) +
            geom_rect(data=all.data[1,], aes(xmin=fb, xmax=fc, ymin=tb, ymax=tc),
                      alpha=0.2, fill="blue", linetype=0) +
            geom_rect(data=all.data[1,], aes(xmin=fc, xmax=fa, ymin=tb, ymax=tc),
                      alpha=0.2, fill="gray", linetype=0) +
            theme(panel.background=element_rect(fill='white', colour='black')) +
            theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
            geom_point(aes(colour=filenames, size=count)) +
            scale_size_continuous(range=c(2,8)) +
            xlab("False Positives") +
            ylab("True Positives") +
            ggtitle("False Positives by True Positves") +
            xlim(0, fa) + ylim(0, ta)
        print(p2)
        }
        
    })
    
    output$meanresults = renderTable({apply(all.data()[,1:5],2,mean)},rownames=TRUE)
    
    output$sumresults = renderTable(apply(all.data()[,6:9],2,sum),rownames=TRUE)
    
    output$meanresults2 = renderTable(apply(all.data()[,10:16],2,mean),rownames=TRUE)
    
    
    
  }#ends server
  
  
  )#ends shinyApp


