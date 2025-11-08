#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyWidgets)
library(ggplot2)

# Define UI for application that draws a plot
ui <- fluidPage(
    
    # Application title
    titlePanel("Better Gini Impact Calculator"),
    
    # Sidebar with a slider input
    sidebarLayout(
        sidebarPanel(
            tags$head(
                # Note the wrapping of the string in HTML()
                tags$style(HTML("
      @import url('https://fonts.googleapis.com/css2?family=Yusei+Magic&display=swap');
      body {
        # background-color: gray;
        # color: black;
      }
      h3 {
        # font-family: 'Yusei Magic', sans-serif;
        font-size: 20px;
      }
      h4 {
        font-size: 18px;
      }

      .shiny-input-container {
       # color: #00ff00;
      }"))
            ),
            selectInput("METHOD", "Curve function", 
                        choices = c("MidNormal", "MidFractal", "BiFractal", "BiNormal"),
                        selected = "MidNormal"),
            setSliderColor(c("dodgerblue2", "olivedrab", "lightcoral", "firebrick", "gold", "lightblue"),c(1,2,3,4, 5,6)),
            sliderInput("GINI1",
                        "GINI 1:",
                        min = 0,
                        max = 1,
                        value = 0.45),
            sliderInput("GINI2",
                        "GINI 2:",
                        min = 0,
                        max = 1,
                        value = 0.65),
            sliderInput("B",
                        "Population bad rate:",
                        min = 0,
                        max = 100,
                        value = 10),
            sliderInput("a0",
                        "Initial approval rate:",
                        min = 0,
                        max = 100,
                        value = 60),
            conditionalPanel(
                condition = "input.METHOD == 'BiFractal'",
                sliderInput("BETA",
                            "Beta:",
                            min = 0,
                            max = 1,
                            value = 0.5)
            ),
            conditionalPanel(
                condition = "input.METHOD == 'BiNormal'",
                sliderInput("SHAPE",
                            "Shape:",
                            min = 0,
                            max = 1,
                            value = 0.5,
                            step = 0.000001)
                
            )
            
            
        ),
        
        mainPanel(
            plotOutput("curvePlot", height = "600", width = 600)
        ),
        
    ), br(), br(), br(),
    
    
    wellPanel(
        htmlOutput("text1")
        # verbatimTextOutput("verb1"),
        # textOutput("text2"),
        # verbatimTextOutput("verb2"),
        # textOutput("text3"),
        # verbatimTextOutput("verb3"),
        # textOutput("text4"),
        # verbatimTextOutput("verb4"),
        # textOutput("text5"),
        # verbatimTextOutput("verb5")
    )
    # ,
    # submitButton("Create a plot!")
)

# FUNCTIONS:
FuncMidNormal<-function(x,g){pnorm(qnorm((g+1)/2)*sqrt(2)+qnorm(x))}
FuncMidFractal<-function(x,g){0.5*(1-(1-x)^((1+g)/(1-g)))+0.5*x^((1-g)/(1+g))}

FuncBiFractal<-function(x,g, beta){beta*(1-(1-x)^((1+g)/(1-g)))+(1-beta)*x^((1-g)/(1+g))}
FuncBiNormal<-function(x,g,shape){pnorm(qnorm((g+1)/2)*sqrt(1+shape^2)+shape*qnorm(x))}

server <- function(input, output) {
    
    # G1_react <- reactive({get(input$GINI1, 'package:datasets')})
    GINI1<-reactive({
        validate(
            need(input$GINI1 < 1 && input$GINI1 > 0, "Please select GINI1 value between 0 and 1")
        )
        input$GINI1})
    GINI2<-reactive({
        validate(
            need(input$GINI2 < 1 && input$GINI2 > 0, "Please select GINI2 value between 0 and 1")
        )
        input$GINI2})
    B <- reactive({
        validate(
            need(input$B < 100, "Please select B value smaller than 100 %")
        )
        input$B/100})
    a0 <- reactive({
        validate(
            need(input$a0 > 0, "Please select a0 value greater than 0 %")
        )
        input$a0/100})
    
    beta <- reactive({
        validate(
            need(input$BETA > 0 && input$BETA < 1, "Please select BETA value between 0 and 1")
        )
        input$BETA})
    
    shape <- reactive({
        validate(
            need(input$SHAPE > 0, "Please select shape value between 0 and 1")
        )
        -((log(1-input$SHAPE)) / log(2) )  })
    
    
    method <- reactive({input$METHOD})
    
    # MID NORMAL
    y0<-function(x){FuncMidNormal(x,GINI1())}
    y1<-function(x){FuncMidNormal(x,GINI2())}
    
    phi1<-function(x){(1-B())*(1-x)+B()*(1-y0(x))-a0()}
    x0<-reactive({as.numeric(uniroot(phi1,lower=0,upper=1,tol = .Machine$double.eps))[1]})
    
    phi2<-function(x){(1-B())*(1-x)+B()*(1-y1(x))-a0()}
    x1<-reactive({as.numeric(uniroot(phi2,lower=0,upper=1,tol = .Machine$double.eps))[1]})
    
    b0 <- reactive({B()*(1-y0(x0()))/a0()})
    b1 <- reactive({B()*(1-y1(x1()))/a0()})
    
    phi3<-function(x){B()*(1-y1(x))/((1-B())*(1-x)+B()*(1-y1(x)))-b0()}
    x1prime <- reactive({as.numeric(uniroot(phi3,lower=0,upper=.99999999999,tol = .Machine$double.eps))[1]})
    a1 <- reactive({((1-B())*(1-x1prime())+B()*(1-y1(x1prime())))})
    
    
    # MID FRACTAL
    y0_mf<-function(x){FuncMidFractal(x,GINI1())}
    y1_mf<-function(x){FuncMidFractal(x,GINI2())}
    
    phi1_mf<-function(x){(1-B())*(1-x)+B()*(1-y0_mf(x))-a0()}
    x0_mf<-reactive({as.numeric(uniroot(phi1_mf,lower=0,upper=1,tol = .Machine$double.eps))[1]})
    
    phi2_mf<-function(x){(1-B())*(1-x)+B()*(1-y1_mf(x))-a0()}
    x1_mf<-reactive({as.numeric(uniroot(phi2_mf,lower=0,upper=1,tol = .Machine$double.eps))[1]})
    
    b0_mf <- reactive({B()*(1-y0_mf(x0_mf()))/a0()})
    b1_mf <- reactive({B()*(1-y1_mf(x1_mf()))/a0()})
    
    phi3_mf<-function(x){B()*(1-y1_mf(x))/((1-B())*(1-x)+B()*(1-y1_mf(x)))-b0_mf()}
    x1prime_mf <- reactive({as.numeric(uniroot(phi3_mf,lower=0,upper=.99999999999,tol = .Machine$double.eps))[1]})
    a1_mf <- reactive({((1-B())*(1-x1prime_mf())+B()*(1-y1_mf(x1prime_mf())))})
    
    # BI FRACTAL
    y0_bf<-function(x){FuncBiFractal(x,GINI1(), beta())}
    y1_bf<-function(x){FuncBiFractal(x,GINI2(), beta())}
    
    phi1_bf<-function(x){(1-B())*(1-x)+B()*(1-y0_bf(x))-a0()}
    x0_bf<-reactive({as.numeric(uniroot(phi1_bf,lower=0,upper=1,tol = .Machine$double.eps))[1]})
    
    phi2_bf<-function(x){(1-B())*(1-x)+B()*(1-y1_bf(x))-a0()}
    x1_bf<-reactive({as.numeric(uniroot(phi2_bf,lower=0,upper=1,tol = .Machine$double.eps))[1]})
    
    b0_bf <- reactive({B()*(1-y0_bf(x0_bf()))/a0()})
    b1_bf <- reactive({B()*(1-y1_bf(x1_bf()))/a0()})
    
    phi3_bf<-function(x){B()*(1-y1_bf(x))/((1-B())*(1-x)+B()*(1-y1_bf(x)))-b0_bf()}
    x1prime_bf <- reactive({as.numeric(uniroot(phi3_bf,lower=0,upper=.99999999999,tol = .Machine$double.eps))[1]})
    a1_bf <- reactive({((1-B())*(1-x1prime_bf())+B()*(1-y1_bf(x1prime_bf())))})
    
    # BI NORMAL
    y0_bn<-function(x){FuncBiNormal(x,GINI1(), shape())}
    y1_bn<-function(x){FuncBiNormal(x,GINI2(), shape())}
    
    phi1_bn<-function(x){(1-B())*(1-x)+B()*(1-y0_bn(x))-a0()}
    x0_bn<-reactive({as.numeric(uniroot(phi1_bn,lower=0,upper=1,tol = .Machine$double.eps))[1]})
    
    phi2_bn<-function(x){(1-B())*(1-x)+B()*(1-y1_bn(x))-a0()}
    x1_bn<-reactive({as.numeric(uniroot(phi2_bn,lower=0,upper=1,tol = .Machine$double.eps))[1]})
    
    b0_bn <- reactive({B()*(1-y0_bn(x0_bn()))/a0()})
    b1_bn <- reactive({B()*(1-y1_bn(x1_bn()))/a0()})
    
    phi3_bn<-function(x){B()*(1-y1_bn(x))/((1-B())*(1-x)+B()*(1-y1_bn(x)))-b0_bn()}
    x1prime_bn <- reactive({as.numeric(uniroot(phi3_bn,lower=0,upper=.99999999999,tol = .Machine$double.eps))[1]})
    a1_bn <- reactive({((1-B())*(1-x1prime_bn())+B()*(1-y1_bn(x1prime_bn())))})
    
    
    output$curvePlot <- renderPlot({
        if(method() == "MidNormal"){
            xcord <- c(x0(), x1(), x1prime())
            ycord <- c(y0(x0()), y1(x1()), y1(x1prime()))
        }
        else if(method() == "MidFractal"){
            xcord <- c(x0_mf(), x1_mf(), x1prime_mf())
            ycord <- c(y0_mf(x0_mf()), y1_mf(x1_mf()), y1_mf(x1prime_mf()))
            y0 <- y0_mf
            y1 <- y1_mf
        } 
        else if(method() == "BiFractal") {
            xcord <- c(x0_bf(), x1_bf(), x1prime_bf())
            ycord <- c(y0_bf(x0_bf()), y1_bf(x1_bf()), y1_bf(x1prime_bf()))
            y0 <- y0_bf
            y1 <- y1_bf
        }
        else if(method() == "BiNormal") {
            xcord <- c(x0_bn(), x1_bn(), x1prime_bn())
            ycord <- c(y0_bn(x0_bn()), y1_bn(x1_bn()), y1_bn(x1prime_bn()))
            y0 <- y0_bn
            y1 <- y1_bn
        }
        point_id <- c("a", "b", "c")
        df <- data.frame(xcord, ycord, point_id)
        
        
        ggplot(data.frame(x=c(0,1)), aes(x)) +
            stat_function(fun = y0, geom = "line", aes(colour = "y0"), lwd = 1.1) +
            stat_function(fun = y1, geom = "line", aes(colour = "y1"), linetype= "dashed", lwd = 1.1) +
            scale_colour_manual("ROC Curve", values = c("dodgerblue2", "olivedrab"), breaks = waiver(), labels=c("scoring 1", "scoring 2")) +
            geom_point(data = df[1,], aes(xcord, ycord, shape = factor(point_id)),  colour = "black", fill = "dodgerblue2", size = 4, stroke = 1.5) +
            geom_point(data = df[c(2,3),], aes(xcord, ycord, shape = factor(point_id)), colour = "black", fill = "olivedrab", size = 4, stroke = 1.5) +
            scale_shape_manual("Point",values = c('a' = 21, 'b' = 22, 'c' = 24), labels = c("Currently", "Bad rate reduction scenario", "Approval rate improvement scenario")) +
            xlab("cumulative good proportion") + ylab("cumulative bad proportion") + ggtitle("") + 
            theme(legend.position="bottom", legend.box = "vertical", axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold")) +
            guides(color = guide_legend(order = 1))
        # code responsible for the square shape of the plot
    }
    # ,height=reactive(ifelse(!is.null(input$innerWidth),input$innerWidth*2/5,0))
    
    )
    
    
    
    
    
    output$text1 <- renderText({
        if(method() == "MidNormal"){
            a1 <- a1()
            b0 <- b0()
            b1 <- b1()
        }
        else if(method() == "MidFractal"){
            a1 <- a1_mf()
            b0 <- b0_mf()
            b1 <- b1_mf()
        }
        else if(method() == "BiFractal") {
            a1 <- a1_bf()
            b0 <- b0_bf()
            b1 <- b1_bf()
        }
        else if(method() == "BiNormal") {
            a1 <- a1_bn()
            b0 <- b0_bn()
            b1 <- b1_bn()
        }
        HTML(paste("<h3>With the new scoring you can achieve a bad rate reduction of:</h3><br/><h4><b>", paste(round((b1/b0-1)*100,4), " %", sep=''), "</b></h4></br>"),
             paste("<h3>... or approval rate improvement of:</h3><br/><h4><b>", paste(round((a1/a0()-1)*100,4), " %", sep=''), "</b><h4></br>"),
             paste("<h3>Portfolio bad rate with scoring 1:</h3><br/><h4><b>", paste(round(b0*100,4), " %", sep=''), "</b><h4></br>"),
             paste("<h3>Portfolio bad rate with scoring 2:</h3><br/><h4><b>", paste(round(b1*100,4), " %", sep=''), "</b><h4></br>"),
             paste("<h3>Approval rate with scoring 2:</h3><br/><h4><b>", paste(round(a1*100,4)," %", sep=''), "</b><h4></br>")
             
        )
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
