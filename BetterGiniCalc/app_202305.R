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
library(ggrepel)
library(numDeriv)
library(data.table)

# Define UI for application that draws a plot
ui <- fluidPage(
    
    # Application title
    titlePanel("Better Gini Impact Calculator"),
      sidebarLayout(
    # Sidebar with a slider input    
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
            setSliderColor(c("dodgerblue2", "olivedrab", "lightcoral", "firebrick", "gold", "gold", "lightblue", "lightblue"),c(1,2,3,4, 5,6, 7, 8)),
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
                sliderInput("BETA1",
                            "Beta1:",
                            min = 0,
                            max = 1,
                            value = 0.5)
            ),
            conditionalPanel(
              condition = "input.METHOD == 'BiFractal'",
              sliderInput("BETA2",
                          "Beta2:",
                          min = 0,
                          max = 1,
                          value = 0.5)
            ),
            conditionalPanel(
                condition = "input.METHOD == 'BiNormal'",
                sliderInput("SHAPE1",
                            "Shape1:",
                            min = 0.7,
                            max = 1.4,
                            value = 1.0,
                            step = 0.001)

            ),
            conditionalPanel(
              condition = "input.METHOD == 'BiNormal'",
              sliderInput("SHAPE2",
                          "Shape2:",
                          min = 0.7,
                          max = 1.4,
                          value = 1.0,
                          step = 0.001)
              
            )
            

        ),

        mainPanel(
          tabsetPanel(
            tabPanel("Plot", plotOutput("curvePlot", height = "600", width = 600)),
            tabPanel("Tables", 
                     tableOutput("table1"),
                     tableOutput("table2"),
                     tableOutput("table3")
                     )
            
          )
        ),

    ) , br(), br(), br(),
    
    
    # wellPanel(
    #     tableOutput("table1"),
    #     tableOutput("table2"),
    #     tableOutput("table3")
    #     #htmlOutput("text1")
    #     # verbatimTextOutput("verb1"),
    #     # textOutput("text2")
    #     # verbatimTextOutput("verb2"),
    #     # textOutput("text3"),
    #     # verbatimTextOutput("verb3"),
    #     # textOutput("text4"),
    #     # verbatimTextOutput("verb4"),
    #     # textOutput("text5"),
    #     # verbatimTextOutput("verb5")
    # )
    # ,
    # submitButton("Create a plot!")
)

# FUNCTIONS:
FuncMidNormal<-function(x,g){pnorm(qnorm((g+1)/2)*sqrt(2)+qnorm(x))}
FuncMidFractal<-function(x,g){0.5*(1-(1-x)^((1+g)/(1-g)))+0.5*x^((1-g)/(1+g))}

FuncBiFractal<-function(x,g, beta){beta*(1-(1-x)^((1+g)/(1-g)))+(1-beta)*x^((1-g)/(1+g))}
FuncBiNormal<-function(x,g,shape){pnorm(qnorm((g+1)/2)*sqrt(1+shape^2)+shape*qnorm(x))}

GiniP<-function(f, x){2*(integrate(f,x,1)$value-(1-x)*f(x))/((1-x)*(1-f(x)))-1
}

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
    
    beta1 <- reactive({
        validate(
            need(input$BETA1 >= 0 && input$BETA1 <= 1, "Please select BETA1 value between 0 and 1")
        )
        input$BETA1})
    
    beta2 <- reactive({
      validate(
        need(input$BETA2 >= 0 && input$BETA2 <= 1, "Please select BETA2 value between 0 and 1")
      )
      input$BETA2})
    
    shape1 <- reactive({
        validate(
            need(input$SHAPE1 >=0.7 && input$SHAPE1<=1.4, "Please select shape1 value closer to 1")
        )
        input$SHAPE1  })
    
    shape2 <- reactive({
      validate(
        need(input$SHAPE2 >=0.7 && input$SHAPE2<=1.4, "Please select shape2 value closer to 1")
      )
      input$SHAPE2  })
    
    method <- reactive({input$METHOD})
    
    # MID NORMAL
    y0<-function(x){FuncMidNormal(x,GINI1())}
    y1<-function(x){FuncMidNormal(x,GINI2())}
    
    phi0<-function(x){(1-B())*(1-x)+B()*(1-y0(x))-a0()}
    x0<-reactive({as.numeric(uniroot(phi0,lower=0,upper=1,tol = .Machine$double.eps)[1])})
    b0 <- reactive({B()*(1-y0(x0()))/a0()})
    
    deriv0 <-  reactive({numDeriv::grad(y0, x0())})
    mbr0 <- reactive({(1+(1-B())/B()/deriv0())^(-1)})
    ir0 <-  reactive({mbr0()/(1-mbr0())})
    profit0 <- reactive({a0()*(ir0()*(1-b0())-b0())})
    ginip0 <- reactive({GiniP(y0,x0())})
    
    
    a1 <- reactive({a0()})
    phi1<-function(x){(1-B())*(1-x)+B()*(1-y1(x))-a1()}
    x1<-reactive({as.numeric(uniroot(phi1,lower=0,upper=1,tol = .Machine$double.eps)[1])})
    b1 <- reactive({B()*(1-y1(x1()))/a0()})
    
    deriv1 <- reactive({numDeriv::grad(y1, x1())})
    mbr1 <- reactive({(1+(1-B())/B()/deriv1())^(-1)})
    ir1 <- reactive({ir0()})
    profit1 <- reactive({a1()*(ir1()*(1-b1())-b1())})
    ginip1 <- reactive({GiniP(y1,x1())})
    a1_change <- reactive({a1()/a0()-1})
    b1_change <- reactive({b1()/b0()-1})
    profit1_change <- reactive({profit1()/profit0()-1})
    
    b2 <- reactive({b0()})
    phi2<-function(x){B()*(1-y1(x))/((1-B())*(1-x)+B()*(1-y1(x)))-b2()}
    x2 <- reactive({as.numeric(uniroot(phi2,lower=0.001,upper=.999,tol = .Machine$double.eps)[1])})
    a2 <- reactive({((1-B())*(1-x2())+B()*(1-y1(x2())))})
    
    deriv2 <- reactive({numDeriv::grad(y1, x2())})
    mbr2 <- reactive({(1+(1-B())/B()/deriv2())^(-1)})
    ir2 <- reactive({ir0()})
    profit2 <- reactive({a2()*(ir2()*(1-b2())-b2())})
    ginip2 <- reactive({GiniP(y1, x2())})
    
    a2_change <- reactive({a2()/a0()-1})
    b2_change <- reactive({b2()/b0()-1})
    profit2_change <- reactive({profit2()/profit0()-1})
    
   
    ## NEW SCENARIO
    phi3<-function(x){numDeriv::grad(y1, x)-deriv0()}
    x3 <- reactive({as.numeric(uniroot(phi3,lower=0+1/1000,upper=1-1/1000,tol = .Machine$double.eps)[1])})
    a3 <- reactive({((1-B())*(1-x3())+B()*(1-y1(x3())))})
    b3 <- reactive({B()*(1-y1(x3()))/a3()})
    
    
    deriv3 <- reactive({numDeriv::grad(y1, x3())})
    mbr3 <- reactive({(1+(1-B())/B()/deriv3())^(-1)})
    ir3 <- reactive({ir0()})
    profit3 <- reactive({a3()*(ir3()*(1-b3())-b3())})
    ginip3 <- reactive({GiniP(y1,x3())})
    
    a3_change <- reactive({a3()/a0()-1})
    b3_change <- reactive({b3()/b0()-1})
    profit3_change <- reactive({profit3()/profit0()-1})
    
    # MID FRACTAL
    y0_mf<-function(x){FuncMidFractal(x,GINI1())}
    y1_mf<-function(x){FuncMidFractal(x,GINI2())}
    
    # phi1_mf<-function(x){(1-B())*(1-x)+B()*(1-y0_mf(x))-a0()}
    # x0_mf<-reactive({as.numeric(uniroot(phi1_mf,lower=0,upper=1,tol = .Machine$double.eps))[1]})
    # 
    # phi2_mf<-function(x){(1-B())*(1-x)+B()*(1-y1_mf(x))-a0()}
    # x1_mf<-reactive({as.numeric(uniroot(phi2_mf,lower=0,upper=1,tol = .Machine$double.eps))[1]})
    # 
    # b0_mf <- reactive({B()*(1-y0_mf(x0_mf()))/a0()})
    # b1_mf <- reactive({B()*(1-y1_mf(x1_mf()))/a0()})
    # 
    # phi3_mf<-function(x){B()*(1-y1_mf(x))/((1-B())*(1-x)+B()*(1-y1_mf(x)))-b0_mf()}
    # x1prime_mf <- reactive({as.numeric(uniroot(phi3_mf,lower=0,upper=.99999999999,tol = .Machine$double.eps))[1]})
    # a1_mf <- reactive({((1-B())*(1-x1prime_mf())+B()*(1-y1_mf(x1prime_mf())))})
    # 
    # phi_yn_x3_mf<-function(x){numDeriv::grad(y1, x)-deriv0()}
    # x3_mf <- reactive({as.numeric(uniroot(phi_yn_x3_mf,lower=0+1/1000,upper=1-1/1000,tol = .Machine$double.eps))[1]})
    # a3_mf <- reactive({((1-B())*(1-x3())+B()*(1-y1(x3())))})
    # 
    ########
    phi0_mf<-function(x){(1-B())*(1-x)+B()*(1-y0_mf(x))-a0()}
    x0_mf<-reactive({as.numeric(uniroot(phi0_mf,lower=0,upper=1,tol = .Machine$double.eps)[1])})
    b0_mf <- reactive({B()*(1-y0_mf(x0_mf()))/a0()})
    
    deriv0_mf <-  reactive({numDeriv::grad(y0_mf, x0_mf())})
    mbr0_mf <- reactive({(1+(1-B())/B()/deriv0_mf())^(-1)})
    ir0_mf <-  reactive({mbr0_mf()/(1-mbr0_mf())})
    profit0_mf <- reactive({a0()*(ir0_mf()*(1-b0_mf())-b0_mf())})
    ginip0_mf <- reactive({GiniP(y0_mf,x0_mf())})
    
    
    a1_mf <- reactive({a0()})
    phi1_mf<-function(x){(1-B())*(1-x)+B()*(1-y1_mf(x))-a1_mf()}
    x1_mf<-reactive({as.numeric(uniroot(phi1_mf,lower=0,upper=1,tol = .Machine$double.eps)[1])})
    b1_mf <- reactive({B()*(1-y1_mf(x1_mf()))/a0()})
    
    deriv1_mf <- reactive({numDeriv::grad(y1_mf, x1_mf())})
    mbr1_mf <- reactive({(1+(1-B())/B()/deriv1_mf())^(-1)})
    ir1_mf <- reactive({ir0_mf()})
    profit1_mf <- reactive({a1_mf()*(ir1_mf()*(1-b1_mf())-b1_mf())})
    ginip1_mf <- reactive({GiniP(y1_mf,x1_mf())})
    a1_change_mf <- reactive({a1_mf()/a0()-1})
    b1_change_mf <- reactive({b1_mf()/b0_mf()-1})
    profit1_change_mf <- reactive({profit1_mf()/profit0_mf()-1})
    
    b2_mf <- reactive({b0_mf()})
    phi2_mf <- function(x){B()*(1-y1_mf(x))/((1-B())*(1-x)+B()*(1-y1_mf(x)))-b2_mf()}
    x2_mf <- reactive({as.numeric(uniroot(phi2_mf,lower=0.001,upper=.999,tol = .Machine$double.eps)[1])})
    a2_mf <- reactive({((1-B())*(1-x2_mf())+B()*(1-y1_mf(x2_mf())))})
    
    deriv2_mf <- reactive({numDeriv::grad(y1_mf, x2_mf())})
    mbr2_mf <- reactive({(1+(1-B())/B()/deriv2_mf())^(-1)})
    ir2_mf <- reactive({ir0_mf()})
    profit2_mf <- reactive({a2_mf()*(ir2_mf()*(1-b2_mf())-b2_mf())})
    ginip2_mf <- reactive({GiniP(y1_mf, x2_mf())})
    
    a2_change_mf <- reactive({a2_mf()/a0()-1})
    b2_change_mf <- reactive({b2_mf()/b0_mf()-1})
    profit2_change_mf <- reactive({profit2_mf()/profit0_mf()-1})
    
   
    phi3_mf <- function(x){numDeriv::grad(y1_mf, x)-deriv0_mf()}
    x3_mf <- reactive({as.numeric(uniroot(phi3_mf,lower=0+1/1000,upper=1-1/1000,tol = .Machine$double.eps)[1])})
    a3_mf <- reactive({((1-B())*(1-x3_mf())+B()*(1-y1_mf(x3_mf())))})
    b3_mf <- reactive({B()*(1-y1_mf(x3_mf()))/a3_mf()})
    
    
    deriv3_mf <- reactive({numDeriv::grad(y1_mf, x3_mf())})
    mbr3_mf <- reactive({(1+(1-B())/B()/deriv3_mf())^(-1)})
    ir3_mf <- reactive({ir0_mf()})
    profit3_mf <- reactive({a3_mf()*(ir3_mf()*(1-b3_mf())-b3_mf())})
    ginip3_mf <- reactive({GiniP(y1_mf,x3_mf())})
    
    a3_change_mf <- reactive({a3_mf()/a0()-1})
    b3_change_mf <- reactive({b3_mf()/b0_mf()-1})
    profit3_change_mf <- reactive({profit3_mf()/profit0_mf()-1})
    ########
    
    # BI FRACTAL
    y0_bf<-function(x){FuncBiFractal(x,GINI1(), beta1())}
    y1_bf<-function(x){FuncBiFractal(x,GINI2(), beta2())}
    # 
    # phi1_bf<-function(x){(1-B())*(1-x)+B()*(1-y0_bf(x))-a0()}
    # x0_bf<-reactive({as.numeric(uniroot(phi1_bf,lower=0,upper=1,tol = .Machine$double.eps))[1]})
    # 
    # phi2_bf<-function(x){(1-B())*(1-x)+B()*(1-y1_bf(x))-a0()}
    # x1_bf<-reactive({as.numeric(uniroot(phi2_bf,lower=0,upper=1,tol = .Machine$double.eps))[1]})
    # 
    # b0_bf <- reactive({B()*(1-y0_bf(x0_bf()))/a0()})
    # b1_bf <- reactive({B()*(1-y1_bf(x1_bf()))/a0()})
    # 
    # phi3_bf<-function(x){B()*(1-y1_bf(x))/((1-B())*(1-x)+B()*(1-y1_bf(x)))-b0_bf()}
    # x1prime_bf <- reactive({as.numeric(uniroot(phi3_bf,lower=0,upper=.99999999999,tol = .Machine$double.eps))[1]})
    # a1_bf <- reactive({((1-B())*(1-x1prime_bf())+B()*(1-y1_bf(x1prime_bf())))})
    # 
    # phi_yn_x3_bf<-function(x){numDeriv::grad(y1, x)-deriv0()}
    # x3_bf <- reactive({as.numeric(uniroot(phi_yn_x3_bf,lower=0+1/1000,upper=1-1/1000,tol = .Machine$double.eps))[1]})
    # a3_bf <- reactive({((1-B())*(1-x3())+B()*(1-y1(x3())))})
    # 
    ########
    phi0_bf<-function(x){(1-B())*(1-x)+B()*(1-y0_bf(x))-a0()}
    x0_bf<-reactive({as.numeric(uniroot(phi0_bf,lower=0.001,upper=0.999,tol = .Machine$double.eps)[1])})
    b0_bf <- reactive({B()*(1-y0_bf(x0_bf()))/a0()})
    
    deriv0_bf <-  reactive({numDeriv::grad(y0_bf, x0_bf())})
    mbr0_bf <- reactive({(1+(1-B())/B()/deriv0_bf())^(-1)})
    ir0_bf <-  reactive({mbr0_bf()/(1-mbr0_bf())})
    profit0_bf <- reactive({a0()*(ir0_bf()*(1-b0_bf())-b0_bf())})
    ginip0_bf <- reactive({GiniP(y0_bf,x0_bf())})
    
    
    a1_bf <- reactive({a0()})
    phi1_bf<-function(x){(1-B())*(1-x)+B()*(1-y1_bf(x))-a1_bf()}
    x1_bf<-reactive({as.numeric(uniroot(phi1_bf,lower=0,upper=1,tol = .Machine$double.eps)[1])})
    b1_bf <- reactive({B()*(1-y1_bf(x1_bf()))/a0()})
    
    deriv1_bf <- reactive({numDeriv::grad(y1_bf, x1_bf())})
    mbr1_bf <- reactive({(1+(1-B())/B()/deriv1_bf())^(-1)})
    ir1_bf <- reactive({ir0_bf()})
    profit1_bf <- reactive({a1_bf()*(ir1_bf()*(1-b1_bf())-b1_bf())})
    ginip1_bf <- reactive({GiniP(y1_bf,x1_bf())})
    a1_change_bf <- reactive({a1_bf()/a0()-1})
    b1_change_bf <- reactive({b1_bf()/b0_bf()-1})
    profit1_change_bf <- reactive({profit1_bf()/profit0_bf()-1})
    
    b2_bf <- reactive({b0_bf()})
    phi2_bf <- function(x){B()*(1-y1_bf(x))/((1-B())*(1-x)+B()*(1-y1_bf(x)))-b2_bf()}
    x2_bf <- reactive({as.numeric(uniroot(phi2_bf,lower=0.001,upper=.999,tol = .Machine$double.eps)[1])})
    a2_bf <- reactive({((1-B())*(1-x2_bf())+B()*(1-y1_bf(x2_bf())))})
    
    deriv2_bf <- reactive({numDeriv::grad(y1_bf, x2_bf())})
    mbr2_bf <- reactive({(1+(1-B())/B()/deriv2_bf())^(-1)})
    ir2_bf <- reactive({ir0_bf()})
    profit2_bf <- reactive({a2_bf()*(ir2_bf()*(1-b2_bf())-b2_bf())})
    ginip2_bf <- reactive({GiniP(y1_bf, x2_bf())})
    
    a2_change_bf <- reactive({a2_bf()/a0()-1})
    b2_change_bf <- reactive({b2_bf()/b0_bf()-1})
    profit2_change_bf <- reactive({profit2_bf()/profit0_bf()-1})
    
    
    phi3_bf <- function(x){numDeriv::grad(y1_bf, x)-deriv0_bf()}
    x3_bf <- reactive({as.numeric(uniroot(phi3_bf,lower=0+1/1000,upper=1-1/1000,tol = .Machine$double.eps)[1])})
    a3_bf <- reactive({((1-B())*(1-x3_bf())+B()*(1-y1_bf(x3_bf())))})
    b3_bf <- reactive({B()*(1-y1_bf(x3_bf()))/a3_bf()})
    
    
    deriv3_bf <- reactive({numDeriv::grad(y1_bf, x3_bf())})
    mbr3_bf <- reactive({(1+(1-B())/B()/deriv3_bf())^(-1)})
    ir3_bf <- reactive({ir0_bf()})
    profit3_bf <- reactive({a3_bf()*(ir3_bf()*(1-b3_bf())-b3_bf())})
    ginip3_bf <- reactive({GiniP(y1_bf,x3_bf())})
    
    a3_change_bf <- reactive({a3_bf()/a0()-1})
    b3_change_bf <- reactive({b3_bf()/b0_bf()-1})
    profit3_change_bf <- reactive({profit3_bf()/profit0_bf()-1})
    ########
    # BI NORMAL
    y0_bn<-function(x){FuncBiNormal(x,GINI1(), shape1())}
    y1_bn<-function(x){FuncBiNormal(x,GINI2(), shape2())}
    
    # phi1_bn<-function(x){(1-B())*(1-x)+B()*(1-y0_bn(x))-a0()}
    # x0_bn<-reactive({as.numeric(uniroot(phi1_bn,lower=0,upper=1,tol = .Machine$double.eps))[1]})
    # 
    # phi2_bn<-function(x){(1-B())*(1-x)+B()*(1-y1_bn(x))-a0()}
    # x1_bn<-reactive({as.numeric(uniroot(phi2_bn,lower=0,upper=1,tol = .Machine$double.eps))[1]})
    # 
    # b0_bn <- reactive({B()*(1-y0_bn(x0_bn()))/a0()})
    # b1_bn <- reactive({B()*(1-y1_bn(x1_bn()))/a0()})
    # 
    # phi3_bn<-function(x){B()*(1-y1_bn(x))/((1-B())*(1-x)+B()*(1-y1_bn(x)))-b0_bn()}
    # x1prime_bn <- reactive({as.numeric(uniroot(phi3_bn,lower=0,upper=.99999999999,tol = .Machine$double.eps))[1]})
    # a1_bn <- reactive({((1-B())*(1-x1prime_bn())+B()*(1-y1_bn(x1prime_bn())))})
    # 
    # phi_yn_x3_bn<-function(x){numDeriv::grad(y1, x)-deriv0()}
    # x3_bn <- reactive({as.numeric(uniroot(phi_yn_x3_bn,lower=0+1/1000,upper=1-1/1000,tol = .Machine$double.eps))[1]})
    # a3_bn <- reactive({((1-B())*(1-x3())+B()*(1-y1(x3())))})
    ######## bn
    phi0_bn <- function(x){(1-B())*(1-x)+B()*(1-y0_bn(x))-a0()}
    x0_bn<-reactive({as.numeric(uniroot(phi0_bn,lower=0,upper=1,tol = .Machine$double.eps)[1])})
    b0_bn <- reactive({B()*(1-y0_bn(x0_bn()))/a0()})
    
    deriv0_bn <-  reactive({numDeriv::grad(y0_bn, x0_bn())})
    mbr0_bn <- reactive({(1+(1-B())/B()/deriv0_bn())^(-1)})
    ir0_bn <-  reactive({mbr0_bn()/(1-mbr0_bn())})
    profit0_bn <- reactive({a0()*(ir0_bn()*(1-b0_bn())-b0_bn())})
    ginip0_bn <- reactive({GiniP(y0_bn,x0_bn())})
    
    
    a1_bn <- reactive({a0()})
    phi1_bn<-function(x){(1-B())*(1-x)+B()*(1-y1_bn(x))-a1_bn()}
    x1_bn<-reactive({as.numeric(uniroot(phi1_bn,lower=0,upper=1,tol = .Machine$double.eps)[1])})
    b1_bn <- reactive({B()*(1-y1_bn(x1_bn()))/a0()})
    
    deriv1_bn <- reactive({numDeriv::grad(y1_bn, x1_bn())})
    mbr1_bn <- reactive({(1+(1-B())/B()/deriv1_bn())^(-1)})
    ir1_bn <- reactive({ir0_bn()})
    profit1_bn <- reactive({a1_bn()*(ir1_bn()*(1-b1_bn())-b1_bn())})
    ginip1_bn <- reactive({GiniP(y1_bn,x1_bn())})
    a1_change_bn <- reactive({a1_bn()/a0()-1})
    b1_change_bn <- reactive({b1_bn()/b0_bn()-1})
    profit1_change_bn <- reactive({profit1_bn()/profit0_bn()-1})
    
    b2_bn <- reactive({b0_bn()})
    phi2_bn <- function(x){B()*(1-y1_bn(x))/((1-B())*(1-x)+B()*(1-y1_bn(x)))-b2_bn()}
    x2_bn <- reactive({as.numeric(uniroot(phi2_bn,lower=0.001,upper=.999,tol = .Machine$double.eps))[1]})
    a2_bn <- reactive({((1-B())*(1-x2_bn())+B()*(1-y1_bn(x2_bn())))})
    
    deriv2_bn <- reactive({numDeriv::grad(y1_bn, x2_bn())})
    mbr2_bn <- reactive({(1+(1-B())/B()/deriv2_bn())^(-1)})
    ir2_bn <- reactive({ir0_bn()})
    profit2_bn <- reactive({a2_bn()*(ir2_bn()*(1-b2_bn())-b2_bn())})
    ginip2_bn <- reactive({GiniP(y1_bn, x2_bn())})
    
    a2_change_bn <- reactive({a2_bn()/a0()-1})
    b2_change_bn <- reactive({b2_bn()/b0_bn()-1})
    profit2_change_bn <- reactive({profit2_bn()/profit0_bn()-1})
    
    
    phi3_bn <- function(x){numDeriv::grad(y1_bn, x)-deriv0_bn()}
    x3_bn <- reactive({as.numeric(uniroot(phi3_bn,lower=0+1/1000,upper=1-1/1000,tol = .Machine$double.eps)[1])})
    a3_bn <- reactive({((1-B())*(1-x3_bn())+B()*(1-y1_bn(x3_bn())))})
    b3_bn <- reactive({B()*(1-y1_bn(x3_bn()))/a3_bn()})
    
    
    deriv3_bn <- reactive({numDeriv::grad(y1_bn, x3_bn())})
    mbr3_bn <- reactive({(1+(1-B())/B()/deriv3_bn())^(-1)})
    ir3_bn <- reactive({ir0_bn()})
    profit3_bn <- reactive({a3_bn()*(ir3_bn()*(1-b3_bn())-b3_bn())})
    ginip3_bn <- reactive({GiniP(y1_bn,x3_bn())})
    
    a3_change_bn <- reactive({a3_bn()/a0()-1})
    b3_change_bn <- reactive({b3_bn()/b0_bn()-1})
    profit3_change_bn <- reactive({profit3_bn()/profit0_bn()-1})
    ########
    
    output$curvePlot <- renderPlot({
        if(method() == "MidNormal"){
            xcord <- c(x0(), x1(), x2(), x3())
            ycord <- c(y0(x0()), y1(x1()), y1(x2()), y1(x3()))
            b1_change <- b1_change()
            a2_change <- a2_change()
            profit3_change <- profit3_change()
        }
        else if(method() == "MidFractal"){
            xcord <- c(x0_mf(), x1_mf(), x2_mf(), x3_mf())
            ycord <- c(y0_mf(x0_mf()), y1_mf(x1_mf()), y1_mf(x2_mf()), y1_mf(x3_mf()))
            y0 <- y0_mf
            y1 <- y1_mf
            b1_change <- b1_change_mf()
            a2_change <- a2_change_mf()
            profit3_change <- profit3_change_mf()
        } 
        else if(method() == "BiFractal") {
            xcord <- c(x0_bf(), x1_bf(), x2_bf(), x3_bf())
            ycord <- c(y0_bf(x0_bf()), y1_bf(x1_bf()), y1_bf(x2_bf()), y1_bf(x3_bf()))
            y0 <- y0_bf
            y1 <- y1_bf
            b1_change <- b1_change_bf()
            a2_change <- a2_change_bf()
            profit3_change <- profit3_change_bf()
        }
        else if(method() == "BiNormal") {
            xcord <- c(x0_bn(), x1_bn(), x2_bn(), x3_bn())
            ycord <- c(y0_bn(x0_bn()), y1_bn(x1_bn()), y1_bn(x2_bn()), y1_bn(x3_bn()))
            y0 <- y0_bn
            y1 <- y1_bn
            b1_change <- b1_change_bn()
            a2_change <- a2_change_bn()
            profit3_change <- profit3_change_bn()
        }
        point_id <- c("a", "b", "c", "d")
        df <- data.frame(xcord, ycord, point_id)
        
        
        ggplot(data.frame(x=c(0,1)), aes(x)) +
            stat_function(fun = y0, geom = "line", aes(colour = "y0"), lwd = 1.1) +
            stat_function(fun = y1, geom = "line", aes(colour = "y1"), linetype= "dashed", lwd = 1.1) +
            scale_colour_manual("ROC Curve", values = c("dodgerblue2", "olivedrab"), breaks = waiver(), labels=c("scoring 1", "scoring 2")) +
            geom_point(data = df[1,], aes(xcord, ycord, shape = factor(point_id)),  colour = "black", fill = "dodgerblue2", size = 4, stroke = 1.5) +
            geom_point(data = df[c(2,3, 4),], aes(xcord, ycord, shape = factor(point_id)), colour = "black", fill = "olivedrab", size = 4, stroke = 1.5) +
            scale_shape_manual("Point",values = c('a' = 21, 'b' = 22, 'c' = 24, 'd' = 23), labels = c("Currently", "Bad rate reduction scenario", "Approval rate improvement scenario", "Profit increase scenario")) +
            xlab("cumulative good proportion") + ylab("cumulative bad proportion") + ggtitle("") + 
            theme(legend.position="bottom", legend.box = "vertical", axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold")) +
            guides(color = guide_legend(order = 1)) + 
            geom_text(x=xcord[2], y=ycord[2], label = paste(round(100*b1_change, 3), "%", sep = " "), hjust = 1, vjust = -1, aes(size = 2), show.legend = FALSE ) +
            geom_text(x=xcord[3], y=ycord[3], label = paste(round(100*a2_change, 3), "%", sep = " "), hjust = 1, vjust = -1, aes(size = 2), show.legend = FALSE ) +
            geom_text(x=xcord[4], y=ycord[4], label = paste(round(100*profit3_change, 3), "%", sep = " "), hjust = 1, vjust = -1, aes(size = 2), show.legend = FALSE )
        # # code responsible for the square shape of the plot
    }
    # ,height=reactive(ifelse(!is.null(input$innerWidth),input$innerWidth*2/5,0))
    
    )
    
    
    output$table1 <- renderTable(
      if(method() == "MidNormal"){
      setDT(
      data.frame(no=1:3, scenario=c('keep approval', 'keep bad rate', 'keep marginal badrate'), approval_change=c(paste(round(100*a1_change(), 3), "%", sep = " "), paste(round(100*a2_change(), 3), "%", sep = " "), paste(round(100*a3_change(), 3), "%", sep = " ")),
                 badrate_change=c(paste(round(100*b1_change(), 3),"%", sep= " "), paste(round(100*b2_change(), 3),"%", sep= " "), paste(round(100*b3_change(), 3),"%", sep= " ")),
                 profit_change=c(paste(round(100*profit1_change(), 3), "%", sep = " "), paste(round(100*profit2_change(), 3), "%", sep = " "), paste(round(100*profit3_change(), 3), "%", sep = " "))
                 )
    
      ) 
      }
      else if(method() == "MidFractal") {
        setDT(
          data.frame(no=1:3, scenario=c('keep approval', 'keep bad rate', 'keep marginal badrate'), approval_change=c(paste(round(100*a1_change_mf(), 3), "%", sep = " "), paste(round(100*a2_change_mf(), 3), "%", sep = " "), paste(round(100*a3_change_mf(), 3), "%", sep = " ")),
                     badrate_change=c(paste(round(100*b1_change_mf(), 3), "%", sep = " "), paste(round(100*b2_change_mf(), 3), "%", sep = " "), paste(round(100*b3_change_mf(), 3), "%", sep = " ")),
                     profit_change=c(paste(round(100*profit1_change_mf(), 3), "%", sep = " "), paste(round(100*profit2_change_mf(), 3), "%", sep = " "), paste(round(100*profit3_change_mf(), 3), "%", sep = " "))
          )
          
        )
      }
      else if(method() == "BiFractal") {
        setDT(
          data.frame(no=1:3, scenario=c('keep approval', 'keep bad rate', 'keep marginal badrate'), approval_change=c(paste(round(100*a1_change_bf(), 3), "%", sep = " "), paste(round(100*a2_change_bf(), 3), "%", sep = " "), paste(round(100*a3_change_bf(), 3), "%", sep = " ")),
                     badrate_change=c(paste(round(100*b1_change_bf(), 3), "%", sep = " "), paste(round(100*b2_change_bf(), 3), "%", sep = " "), paste(round(100*b3_change_bf(), 3), "%", sep = " ")),
                     profit_change=c(paste(round(100*profit1_change_bf(), 3), "%", sep = " "), paste(round(100*profit2_change_bf(), 3), "%", sep = " "), paste(round(100*profit3_change_bf(), 3), "%", sep = " "))
          )
          
        )
      }
      else if(method() == "BiNormal") {
        setDT(
          data.frame(no=1:3, scenario=c('keep approval', 'keep bad rate', 'keep marginal bad rate'), approval_change=c(paste(round(100*a1_change_bn(), 3), "%", sep = " "), paste(round(100*a2_change_bn(), 3), "%", sep = " "), paste(round(100*a3_change_bn(), 3), "%", sep = " ")),
                     badrate_change=c(paste(round(100*b1_change_bn(), 3), "%", sep = " "), paste(round(100*b2_change_bn(), 3), "%", sep = " "), paste(round(100*b3_change_bn(), 3), "%", sep = " ")),
                     profit_change=c(paste(round(100*profit1_change_bn(), 3), "%", sep = " "), paste(round(100*profit2_change_bn(), 3), "%", sep = " "), paste(round(100*profit3_change_bn(), 3), "%", sep = " "))
          )
          
        )
      }
    )
    
    
    
    
    output$table2 <- renderTable(
      if(method() == "MidNormal"){
        setDT(
          data.frame(scenario=1:3, existing_approval=rep(paste(round(100*a0(), 3), "%", sep = " "),3), 
                     new_approval= c(paste(round(100*a1(), 3), "%", sep = " "), paste(round(100*a2(), 3), "%", sep = " "), paste(round(100*a3(), 3), "%", sep = " ")),
                     existing_badrate= rep(paste(round(100*b0(), 3), "%", sep = " "),3),
                     new_badrate= c(paste(round(100*b1(), 3), "%", sep = " "), paste(round(100*b2(), 3), "%", sep = " "), paste(round(100*b3(), 3), "%", sep = " ")),
                     existing_profit= 100*rep(profit0(), 3),
                     new_profit= 100*c(profit1(), profit2(), profit3())
          )
          
        )
      }
      else if(method() == "MidFractal") {
        setDT(
          data.frame(scenario=1:3, existing_approval=rep(paste(round(100*a0(), 3), "%", sep = " "),3), 
                     new_approval=c(paste(round(100*a1_mf(), 3), "%", sep = " "), paste(round(100*a2_mf(), 3), "%", sep = " "), paste(round(100*a3_mf(), 3), "%", sep = " ")),
                     existing_badrate=rep(paste(round(100*b0_mf(), 3), "%", sep = " "),3),
                     new_badrate=c(paste(round(100*b1_mf(), 3), "%", sep = " "),paste(round(100*b2_mf(), 3), "%", sep = " "),paste(round(100*b3_mf(), 3), "%", sep = " ")),
                     existing_profit=100*rep(profit0_mf(), 3),
                     new_profit=100*c(profit1_mf(),profit2_mf(),profit3_mf())
          )
        )
      }
      else if(method() == "BiFractal") {
        setDT(
        data.frame(scenario=1:3, existing_approval=rep(paste(round(100*a0(), 3), "%", sep = " "),3), 
                   new_approval=c(paste(round(100*a1_bf(), 3), "%", sep = " "), paste(round(100*a2_bf(), 3), "%", sep = " "), paste(round(100*a3_bf(), 3), "%", sep = " ")),
                   existing_badrate=rep(paste(round(100*b0_bf(), 3), "%", sep = " "),3),
                   new_badrate=c(paste(round(100*b1_bf(), 3), "%", sep = " "),paste(round(100*b2_bf(), 3), "%", sep = " "),paste(round(100*b3_bf(), 3), "%", sep = " ")),
                   existing_profit=100*rep(profit0_bf(), 3),
                   new_profit=100*c(profit1_bf(), profit2_bf(), profit3_bf())
          )
        )
      }
      else if(method() == "BiNormal") {
        setDT(
        data.frame(scenario=1:3, existing_approval=rep(paste(round(100*a0(), 3), "%", sep = " "),3), 
                   new_approval=c(paste(round(100*a1_bn(), 3), "%", sep = " "), paste(round(100*a2_bn(), 3), "%", sep = " "), paste(round(100*a3_bn(), 3), "%", sep = " ")),
                   existing_badrate=rep(paste(round(100*b0_bn(), 3), "%", sep = " "),3),
                   new_badrate=c(paste(round(100*b1_bn(), 3), "%", sep = " "),paste(round(100*b2_bn(), 3), "%", sep = " "),paste(round(100*b3_bn(), 3), "%", sep = " ")),
                   existing_profit=100*rep(profit0_bn(), 3),
                   new_profit=100*c(profit1_bn(), profit2_bn(), profit3_bn())
        )    
        )
      }
     )
    
    output$table3 <- renderTable(
      if(method() == "MidNormal"){
        setDT(
          data.frame(scenario=1:3, 
                     existing_portfolio_gini = rep(ginip0(),3),
                     new_portfolio_gini = c(ginip1(), ginip2(), ginip3()),
                     existing_marginal_bad_rate = rep(paste(round(100*mbr0(), 3), "%", sep = " "),3)
                     ,new_marginal_bad_rate = c(paste(round(100*mbr1(), 3), "%", sep = " "), paste(round(100*mbr2(), 3), "%", sep = " "), paste(round(100*mbr3(), 3), "%", sep = " "))
          )
          
        )
      }
      else if(method() == "MidFractal") {
        setDT(
          data.frame(scenario=1:3, 
                     existing_portfolio_gini = rep(ginip0_mf(),3),
                     new_portfolio_gini = c(ginip1_mf(), ginip2_mf(), ginip3_mf()),
                     existing_marginal_bad_rate = rep(paste(round(100*mbr0_mf(), 3), "%", sep = " "),3), 
                     new_marginal_bad_rate = c(paste(round(100*mbr1_mf(), 3), "%", sep = " "), paste(round(100*mbr2_mf(), 3), "%", sep = " "), paste(round(100*mbr3_mf(), 3), "%", sep = " "))
          )
        )
      }
      else if(method() == "BiFractal") {
        setDT(
          data.frame(scenario=1:3, 
                     existing_portfolio_gini = rep(ginip0_bf(),3),
                     new_portfolio_gini = c(ginip1_bf(), ginip2_bf(), ginip3_bf()),
                     existing_marginal_bad_rate = rep(paste(round(100*mbr0_bf(), 3), "%", sep = " "),3), 
                     new_marginal_bad_rate = c(paste(round(100*mbr1_bf(), 3), "%", sep = " "), paste(round(100*mbr2_bf(), 3), "%", sep = " "), paste(round(100*mbr3_bf(), 3), "%", sep = " "))
          )
        )
      }
      else if(method() == "BiNormal") {
        setDT(
          data.frame(scenario=1:3, 
                     existing_portfolio_gini = rep(ginip0_bn(),3),
                     new_portfolio_gini = c(ginip1_bn(), ginip2_bn(), ginip3_bn()),
                     existing_marginal_bad_rate = rep(paste(round(100*mbr0_bn(), 3), "%", sep = " "),3), 
                     new_marginal_bad_rate = c(paste(round(100*mbr1_bn(), 3), "%", sep = " "), paste(round(100*mbr2_bn(), 3), "%", sep = " "), paste(round(100*mbr3_bn(), 3), "%", sep = " "))
          ) 
        )
      }
    )
    
    output$text1 <- renderText({

        # a1_change <- 100*a1_change()
        # a2_change <- 100*a2_change()
        # a3_change <- 100*a3_change()
        # b1_change <- 100*b1_change() 
        # b2_change <- 100*b2_change()
        # b3_change <- 100*b3_change()
        # profit1_change <- 100*profit1_change()
        # profit2_change <- 100*profit2_change()
        # profit3_change <- 100*profit3_change()
        # 
        # HTML(paste("<h3>With the new scoring you can achieve a bad rate reduction of:</h3><br/><h4><b>", paste(round((b1/b0-1)*100,4), " %", sep=''), "</b></h4></br>"),
        #      paste("<h3>... or approval rate improvement of:</h3><br/><h4><b>", paste(round((a1/a0()-1)*100,4), " %", sep=''), "</b><h4></br>"),
        #      paste("<h3>Portfolio bad rate with scoring 1:</h3><br/><h4><b>", paste(round(b0*100,4), " %", sep=''), "</b><h4></br>"),
        #      paste("<h3>Portfolio bad rate with scoring 2:</h3><br/><h4><b>", paste(round(b1*100,4), " %", sep=''), "</b><h4></br>"),
        #      paste("<h3>Approval rate with scoring 2:</h3><br/><h4><b>", paste(round(a1*100,4)," %", sep=''), "</b><h4></br>"),
        # 
        # 
        # )
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

