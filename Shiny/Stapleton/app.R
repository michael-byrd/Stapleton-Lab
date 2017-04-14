library(shiny)
# Define UI for application that draws a histogram
ui <- fluidPage(

   # Application title
   titlePanel("Expectation and Variance"),

   # Sidebar with a slider input for number of bins
   sidebarLayout(
      sidebarPanel(
         #For the radio buttons, color and black/white are offered. color is noted "rb" black/white is"bw"
         radioButtons(inputId="clr",label="Color Choices", selected = "bw",choiceNames = c("Black and White","Color"),choiceValues = c("bw","rb"))
        )
      ,
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("plt",width="100%",height= "400px")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
      
  #Render the Plot
   output$plt <- renderPlot({
     #Trying to recreate the Difference between mean and variance shiny plot
     library(ggplot2)
     library(reshape2)

     #generating 100 random variables for 1dt data set
     df=data.frame(
       ord=c(1:100),
       GroupA=rnorm(100, mean=3, sd=(3/6)),
       GroupB=rnorm(100, mean=3, sd=(3/6))
     )
     #alternate data frame
     dfa=data.frame(
       ord=c(151:250),
       GroupA=rnorm(100, mean=1.2*.6*3, sd=1.2*.6*(3/6)),
       GroupB=rnorm(100, mean=.8*.6*3, sd=.8*.6*(3/6))
     )
     #melting data
     md<- melt(df[,c(1,2:3)], id="ord")
     mda<- melt(dfa[,c(1,2:3)], id="ord")
     #set up for the X axis labels
     breaks.major <- c(-10,50,200,260)
     breaks.minor <- c(0,100,110,140,150,250)
     labels.minor <- c("","Control", "Stress","")
     lims <-c(-10,260)
     #contructing charts
     #check the radio buttons
     colr<- ifelse(input$clr=="rb",1,0)
     if(colr==1){
     ggplot(xlab="",ylab="") +
       geom_point(data=df,aes(x=ord,y=GroupA,fill="Group A"),shape=21,color="black",stroke=.2*(2+sqrt(3)),alpha=.8,size=2+sqrt(3))+
       geom_point(data=df,aes(x=ord,y=GroupB,fill="Group B"),shape=22,color="black",stroke=.2*(2+sqrt(3)),alpha=.8,size=2+sqrt(3))+
       geom_point(data=dfa,aes(x=ord,y=GroupA,fill="Group A"),shape=21,color="black",stroke=.2*(2+sqrt(3)),alpha=.8,size=2+sqrt(3))+
       geom_point(data=dfa,aes(x=ord,y=GroupB,fill="Group B"),shape=22,color="black",stroke=.2*(2+sqrt(3)),alpha=.8,size=2+sqrt(3))+
       scale_fill_manual(name="", values=c("Group A"="darkblue", "Group B"="red"), guide="legend")+
       geom_segment(aes(x=-10, xend=110, y=3, yend=3))+
       geom_segment(aes(x=135, xend=260, y=.6*3, yend=.6*3))+
       labs(title="Control vs. Stress Environments", face= "bold")+
       guides(fill= guide_legend(override.aes=list(shape=c(21,22))))+
       scale_x_continuous(expand=c(0,0), limit=lims, minor_breaks=breaks.minor, breaks=breaks.major, labels = labels.minor)+
       theme(axis.title.x=element_blank(),
             axis.text.x=element_text(size=18, face="bold", color="black"),
             axis.ticks.x=element_blank(),
             axis.title.y=element_blank(),
             panel.background = element_rect(fill = "white"),
             plot.title = element_text(size = rel(2),hjust=.5),
             legend.background = element_rect(fill= "white"),
             legend.key = element_blank(),
             legend.text= element_text(size=14)
            )
     }
     else{
         #contructing charts
         ggplot(xlab="",ylab="") +
           geom_point(data=df,aes(x=ord,y=GroupA,shape="Group A"),fill="white",color="black",stroke=.2*(2+sqrt(3)),alpha=.8,size=2+sqrt(3))+
           geom_point(data=df,aes(x=ord,y=GroupB,shape="Group B"),fill="white",color="black",stroke=.2*(2+sqrt(3)),alpha=.8,size=2+sqrt(3))+
           geom_point(data=dfa,aes(x=ord,y=GroupA,shape="Group A"),fill="white",color="black",stroke=.2*(2+sqrt(3)),alpha=.8,size=2+sqrt(3))+
           geom_point(data=dfa,aes(x=ord,y=GroupB,shape="Group B"),fill="white",color="black",stroke=.2*(2+sqrt(3)),alpha=.8,size=2+sqrt(3))+
           geom_segment(aes(x=-10, xend=110, y=3, yend=3))+
           geom_segment(aes(x=135, xend=260, y=.6*3, yend=.6*3))+
           labs(title="Control vs. Stress Environments", face= "bold")+
           scale_shape_manual(name='',values=c("Group A"=21,"Group B"=22))+
           guides(shape= guide_legend(override.aes=list(shape=c(21,22))))+
           scale_x_continuous(expand=c(0,0), limit=lims, minor_breaks=breaks.minor, breaks=breaks.major, labels = labels.minor)+
           theme(axis.title.x=element_blank(),
                 axis.text.x=element_text(size=18, face="bold", color="black"),
                 axis.ticks.x=element_blank(),
                 axis.title.y=element_blank(),
                 panel.background = element_rect(fill = "white"),
                 plot.title = element_text(size = rel(2),hjust=.5),
                 legend.background = element_rect(fill= "white"),
                 legend.key = element_blank(),
                 legend.text= element_text(size=14)
           )
       }
     }
   )
   } 
# Run the application
shinyApp(ui = ui, server = server)