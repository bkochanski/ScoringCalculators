library(shiny)
library(shinyWidgets)
library(ggplot2)
library(ggrepel)
library(numDeriv)
library(data.table)

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
                        step = 0.001,
                        value = 0.45),
            sliderInput("GINI2",
                        "GINI 2:",
                        min = 0,
                        max = 1,
                        step = 0.001,
                        value = 0.65),
            sliderInput("B",
                        "Population bad rate:",
                        min = 0,
                        max = 100,
                        step = 0.1, 
                        value = 10),
            sliderInput("a0",
                        "Initial approval rate:",
                        min = 0,
                        max = 100,
                        step = 0.1, 
                        value = 60),
            conditionalPanel(
                condition = "input.METHOD == 'BiFractal'",
                sliderInput("BETA1",
                            "Shape 1:",
                            min = 0,
                            max = 1,
                            value = 0.5)
            ),
            conditionalPanel(
              condition = "input.METHOD == 'BiFractal'",
              checkboxInput("SAMEBETA",
                          "Same shape",
                          value = TRUE)
            ),
            conditionalPanel(
              condition = "input.METHOD == 'BiFractal' && input.SAMEBETA==0",
              sliderInput("BETA2",
                          "Shape 2:",
                          min = 0,
                          max = 1,
                          value = 0.5)
            ),
            conditionalPanel(
                condition = "input.METHOD == 'BiNormal'",
                sliderInput("SHAPE1",
                            "Shape 1:",
                            min = 0.7,
                            max = 1.4,
                            value = 1.0,
                            step = 0.001)
            ),
            conditionalPanel(
              condition = "input.METHOD == 'BiNormal'",
              checkboxInput("SAMESHAPE",
                            "Same shape",
                            value = TRUE)
            ),            conditionalPanel(
              condition = "input.METHOD == 'BiNormal' && input.SAMESHAPE==0",
              sliderInput("SHAPE2",
                          "Shape 2:",
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
                     ),
            # tabPanel("Instructions", 
            #          mainPanel(
            #            h2("Curve function")
            #            , p("Choose a curve function: Binormal and Bifractal currently available. Midnormal is a binormal with shape parameter set to 1. Midfractal is bifractal with beta parameter set to 0.5")
            #            , a("Bifractal", href="https://www.google.com/")
            #            , 
            #          )
            #          )
            
          )
        ),

    ) , br(), br(), br(),
)

# FUNCTIONS:
FuncMidNormal<-function(x,g){pnorm(qnorm((g+1)/2)*sqrt(2)+qnorm(x))}
FuncMidFractal<-function(x,g){0.5*(1-(1-x)^((1+g)/(1-g)))+0.5*x^((1-g)/(1+g))}

FuncBiFractal<-function(x,g, beta){beta*(1-(1-x)^((1+g)/(1-g)))+(1-beta)*x^((1-g)/(1+g))}
FuncBiNormal<-function(x,g,shape){pnorm(qnorm((g+1)/2)*sqrt(1+shape^2)+shape*qnorm(x))}

# Gini above the cutoff (truncated Gini) calculation
# f - ROC curve model (e.g. binormal with set Gini & shape parameters)
# x - cut-off point representation on the ROC curve (x coordinate: cumulative good proportion below the cutoff)
GiniP<-function(f, x){2*(integrate(f,x,1)$value-(1-x)*f(x))/((1-x)*(1-f(x)))-1}

server <- function(input, output) {
    
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
            need(input$B < 100, "Please select B value smaller than 100%")
        )
        input$B/100})
    a0 <- reactive({
        validate(
            need(input$a0 > 0, "Please select a0 value greater than 0%")
        )
        input$a0/100})
    
    beta1 <- reactive({
        validate(
            need(input$BETA1 >= 0 && input$BETA1 <= 1, "Please select shape 1 value between 0 and 1")
        )
        input$BETA1})
    
    samebeta <- reactive(input$SAMEBETA)
    
    beta2 <- reactive({
      validate(
        need(input$BETA2 >= 0 && input$BETA2 <= 1, "Please select shape 2 value between 0 and 1")
      )
      input$BETA2})

    shape1 <- reactive({
        validate(
            need(input$SHAPE1 >=0.7 && input$SHAPE1<=1.4, "Please select shape 1 value closer to 1")
        )
        input$SHAPE1  })
 
    sameshape <- reactive(input$SAMESHAPE)
       
    shape2 <- reactive({
      validate(
        need(input$SHAPE2 >=0.7 && input$SHAPE2<=1.4, "Please select shape 2 value closer to 1")
      )
      input$SHAPE2  })
    
    method <- reactive({input$METHOD})
    
    ############# 
    # MID NORMAL
    #############
    
    # setting midnormal functions with GINI1 and GINI2 parameters
    y0<-function(x){FuncMidNormal(x,GINI1())}
    y1<-function(x){FuncMidNormal(x,GINI2())}
    
    # searching for representation of the cutoff point 
    phi0<-function(x){(1-B())*(1-x)+B()*(1-y0(x))-a0()}
    x0<-reactive({as.numeric(uniroot(phi0,lower=0,upper=1,tol = .Machine$double.eps)[1])})

    # determining portfolio bad rate for scorecard 1
    b0 <- reactive({B()*(1-y0(x0()))/a0()})
    
    # determining marginal bad rate at the cutoff for scorecard 1 
    deriv0 <-  reactive({numDeriv::grad(y0, x0())})
    mbr0 <- reactive({(1+(1-B())/B()/deriv0())^(-1)})
    
    # determining interest rate assuming profit = 0 at marginal bad rate 
    ir0 <-  reactive({mbr0()/(1-mbr0())})
    
    # calculation of the portfolio profit
    profit0 <- reactive({a0()*(ir0()*(1-b0())-b0())})
    
    # Gini for scorecard 1 on accepted portfolio
    ginip0 <- reactive({GiniP(y0,x0())})
    
    ## 1. Bad rate reduction scenario (keep approval)
    
    # same approval
    a1 <- reactive({a0()})
    
    # finding the cutoff representation on ROC curve 2 for scenario 1
    phi1<-function(x){(1-B())*(1-x)+B()*(1-y1(x))-a1()}
    x1<-reactive({as.numeric(uniroot(phi1,lower=0,upper=1,tol = .Machine$double.eps)[1])})
    
    # new portfolio bad rate
    b1 <- reactive({B()*(1-y1(x1()))/a0()})
    
    # new marginal bad rate at the cutoff
    deriv1 <- reactive({numDeriv::grad(y1, x1())})
    mbr1 <- reactive({(1+(1-B())/B()/deriv1())^(-1)})
    
    # same interest rate
    ir1 <- reactive({ir0()})
    
    # profit for scorecard 2
    profit1 <- reactive({a1()*(ir1()*(1-b1())-b1())})
    
    # calculating Gini on approved portfolio for scorecard 2 
    ginip1 <- reactive({GiniP(y1,x1())})
    
    # approval change, portfolio bad rate change, profit change
    a1_change <- reactive({a1()/a0()-1})
    b1_change <- reactive({b1()/b0()-1})
    profit1_change <- reactive({profit1()/profit0()-1})
    
    ## 2. Approval rate increase scenario (keep bad rate)

    # same portfolio bad rate
    b2 <- reactive({b0()})
    
    # finding the cutoff representation on ROC curve 2 for scenario 2
    phi2<-function(x){B()*(1-y1(x))/((1-B())*(1-x)+B()*(1-y1(x)))-b2()}
    x2 <- reactive({as.numeric(uniroot(phi2,lower=0.001,upper=.999,tol = .Machine$double.eps)[1])})
    
    # new approval rate
    a2 <- reactive({((1-B())*(1-x2())+B()*(1-y1(x2())))})
    
    # as above
    deriv2 <- reactive({numDeriv::grad(y1, x2())})
    mbr2 <- reactive({(1+(1-B())/B()/deriv2())^(-1)})
    ir2 <- reactive({ir0()})
    profit2 <- reactive({a2()*(ir2()*(1-b2())-b2())})
    ginip2 <- reactive({GiniP(y1, x2())})
    
    a2_change <- reactive({a2()/a0()-1})
    b2_change <- reactive({b2()/b0()-1})
    profit2_change <- reactive({profit2()/profit0()-1})
    
    ## 3. Profit increase scenario (keep marginal bad rate)
    
    # same marginal bad rate - finding the cutoff representation on ROC curve 2 for scenario 3 
    phi3<-function(x){numDeriv::grad(y1, x)-deriv0()}
    x3 <- reactive({as.numeric(uniroot(phi3,lower=0+1/1000,upper=1-1/1000,tol = .Machine$double.eps)[1])})
    
    # new approval rate and portfolio bad rate in scenario 3
    a3 <- reactive({((1-B())*(1-x3())+B()*(1-y1(x3())))})
    b3 <- reactive({B()*(1-y1(x3()))/a3()})
    
    # as above
    deriv3 <- reactive({numDeriv::grad(y1, x3())})
    mbr3 <- reactive({(1+(1-B())/B()/deriv3())^(-1)})
    ir3 <- reactive({ir0()})
    profit3 <- reactive({a3()*(ir3()*(1-b3())-b3())})
    ginip3 <- reactive({GiniP(y1,x3())})
    
    a3_change <- reactive({a3()/a0()-1})
    b3_change <- reactive({b3()/b0()-1})
    profit3_change <- reactive({profit3()/profit0()-1})

    ############# 
    # MIDFRACTAL
    #############

    y0_mf<-function(x){FuncMidFractal(x,GINI1())}
    y1_mf<-function(x){FuncMidFractal(x,GINI2())}
    
    ########
    phi0_mf<-function(x){(1-B())*(1-x)+B()*(1-y0_mf(x))-a0()}
    x0_mf<-reactive({as.numeric(uniroot(phi0_mf,lower=0,upper=1,tol = .Machine$double.eps)[1])})
    b0_mf <- reactive({B()*(1-y0_mf(x0_mf()))/a0()})
    
    deriv0_mf <-  reactive({numDeriv::grad(y0_mf, x0_mf())})
    mbr0_mf <- reactive({(1+(1-B())/B()/deriv0_mf())^(-1)})
    ir0_mf <-  reactive({mbr0_mf()/(1-mbr0_mf())})
    profit0_mf <- reactive({a0()*(ir0_mf()*(1-b0_mf())-b0_mf())})
    ginip0_mf <- reactive({GiniP(y0_mf,x0_mf())})
    
    ## 1. Bad rate reduction scenario (keep approval)
    
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

    ## 2. Approval rate increase scenario (keep bad rate)

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

    ## 3. Profit increase scenario (keep marginal bad rate)
    
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
    # BIFRACTAL
    ########
    
    y0_bf<-function(x){FuncBiFractal(x,GINI1(), beta1())}
    y1_bf<-function(x){FuncBiFractal(x,GINI2(), ifelse(samebeta(), beta1(), beta2()))}
    
    #####
    
    phi0_bf<-function(x){(1-B())*(1-x)+B()*(1-y0_bf(x))-a0()}
    x0_bf<-reactive({as.numeric(uniroot(phi0_bf,lower=0.001,upper=0.999,tol = .Machine$double.eps)[1])})
    b0_bf <- reactive({B()*(1-y0_bf(x0_bf()))/a0()})
    
    deriv0_bf <-  reactive({numDeriv::grad(y0_bf, x0_bf())})
    mbr0_bf <- reactive({(1+(1-B())/B()/deriv0_bf())^(-1)})
    ir0_bf <-  reactive({mbr0_bf()/(1-mbr0_bf())})
    profit0_bf <- reactive({a0()*(ir0_bf()*(1-b0_bf())-b0_bf())})
    ginip0_bf <- reactive({GiniP(y0_bf,x0_bf())})

    ## 1. Bad rate reduction scenario (keep approval)
    
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
    
    ## 2. Approval rate increase scenario (keep bad rate)
    
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
    
    ## 3. Profit increase scenario (keep marginal bad rate)
    
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
    # BINORMAL
    ########
    
    y0_bn<-function(x){FuncBiNormal(x,GINI1(), shape1())}
    y1_bn<-function(x){FuncBiNormal(x,GINI2(), ifelse(sameshape(), shape1(), shape2()))}
    
    ######## bn

    phi0_bn <- function(x){(1-B())*(1-x)+B()*(1-y0_bn(x))-a0()}
    x0_bn<-reactive({as.numeric(uniroot(phi0_bn,lower=0,upper=1,tol = .Machine$double.eps)[1])})
    b0_bn <- reactive({B()*(1-y0_bn(x0_bn()))/a0()})
    
    deriv0_bn <-  reactive({numDeriv::grad(y0_bn, x0_bn())})
    mbr0_bn <- reactive({(1+(1-B())/B()/deriv0_bn())^(-1)})
    ir0_bn <-  reactive({mbr0_bn()/(1-mbr0_bn())})
    profit0_bn <- reactive({a0()*(ir0_bn()*(1-b0_bn())-b0_bn())})
    ginip0_bn <- reactive({GiniP(y0_bn,x0_bn())})
    
    ## 1. Bad rate reduction scenario (keep approval)
    
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

    ## 2. Approval rate increase scenario (keep bad rate)
    
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
    
    ## 3. Profit increase scenario (keep marginal bad rate)

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
        df$labels <- c('', paste("[1] Bad rate change:\n", ifelse(b1_change>0,"+",""), round(100*b1_change, 3), "%", sep = ""),
                       paste("[2] Approval change:\n", ifelse(a2_change>0,"+",""), round(100*a2_change, 3), "%", sep = ""),
                       paste("[3] Profit change:\n", ifelse(profit3_change>0,"+",""), round(100*profit3_change, 3), "%", sep = ""))
        
        
        ggplot(data.frame(x=c(0,1)), aes(x)) +
            stat_function(fun = y0, geom = "line", aes(colour = "y0"), lwd = 1.1) +
            stat_function(fun = y1, geom = "line", aes(colour = "y1"), linetype= "dashed", lwd = 1.1) +
            scale_colour_manual("ROC Curve", values = c("dodgerblue2", "olivedrab"), breaks = waiver(), labels=c("Scorecard 1", "Scorecard 2")) +
            geom_point(data = df[1,], aes(xcord, ycord, shape = factor(point_id)),  colour = "black", fill = "dodgerblue2", size = 4, stroke = 1.5) +
            geom_point(data = df[c(2,3, 4),], aes(xcord, ycord, shape = factor(point_id)), colour = "black", fill = "olivedrab", size = 4, stroke = 1.5) +
            scale_shape_manual("Cut-off\npoint",
                               values = c('a' = 21, 'b' = 22, 'c' = 24, 'd' = 23),
                               labels = c("Currently", "[1] Bad rate reduction\n scenario (keep approval)", "[2] Approval rate\n improvement scenario\n (keep bad rate)", "[3] Profit increase\n scenario (keep marginal\n bad rate)")) +
            xlab("cumulative good proportion") + ylab("cumulative bad proportion") + ggtitle("") + 
            theme(legend.position="bottom", 
                  legend.box = "vertical", 
                  axis.text=element_text(size=14), 
                  axis.title=element_text(size=16,face="bold"),
                  legend.text=element_text(size=11)) +
            guides(color = guide_legend(order = 1)) +
            geom_text_repel(data = df[2:4,],
                mapping = aes(x=xcord, y=ycord, label = labels),
                size=5,
                min.segment.length = 0,
                hjust = 1, vjust = -1, show.legend = FALSE)
    }
    )
    
    f3<-function(x){format(round(100*x,3), nsmall=3)}
    fperc<-function(x){paste(f3(x),"%", sep="")}
    
    
    output$table1 <- renderTable(
      if(method() == "MidNormal"){
      setDT(
      data.frame(no=1:3, scenario=c('keep approval (bad rate reduction)', 'keep bad rate (approval rate improvement)', 'keep marginal bad rate (profit increase)'), 
                 approval_change=c(fperc(a1_change()), fperc(a2_change()), fperc(a3_change())),
                 bad_rate_change=c(fperc(b1_change()), fperc(b2_change()), fperc(b3_change())),
                 profit_change=c(fperc(profit1_change()), fperc(profit2_change()), fperc(profit3_change()))
                 )
    
      ) 
      }
      else if(method() == "MidFractal") {
        setDT(
          data.frame(no=1:3, 
                     scenario=c('keep approval (bad rate reduction)', 'keep bad rate (approval rate improvement)', 'keep marginal bad rate (profit increase)'), 
                     approval_change=c(fperc(a1_change_mf()), fperc(a2_change_mf()), fperc(a3_change_mf())),
                     bad_rate_change=c(fperc(b1_change_mf()), fperc(b2_change_mf()), fperc(b3_change_mf())),
                     profit_change=c(fperc(profit1_change_mf()), fperc(profit2_change_mf()), fperc(profit3_change_mf()))
          )
          
        )
      }
      else if(method() == "BiFractal") {
        setDT(
          data.frame(no=1:3, scenario=c('keep approval (bad rate reduction)', 'keep bad rate (approval rate improvement)', 'keep marginal bad rate (profit increase)'), 
                     approval_change=c(fperc(a1_change_bf()), fperc(a2_change_bf()), fperc(a3_change_bf())),
                     bad_rate_change=c(fperc(b1_change_bf()), fperc(b2_change_bf()), fperc(b3_change_bf())),
                     profit_change=c(fperc(profit1_change_bf()), fperc(profit2_change_bf()), fperc(profit3_change_bf()))
          )
          
        )
      }
      else if(method() == "BiNormal") {
        setDT(
          data.frame(no=1:3, scenario=c('keep approval (bad rate reduction)', 'keep bad rate (approval rate improvement)', 'keep marginal bad rate (profit increase)'), 
                     approval_change=c(fperc(a1_change_bn()), fperc(a2_change_bn()), fperc(a3_change_bn())),
                     bad_rate_change=c(fperc(b1_change_bn()), fperc(b2_change_bn()), fperc(b3_change_bn())),
                     profit_change=c(fperc(profit1_change_bn()), fperc(profit2_change_bn()), fperc(profit3_change_bn()))
          )
          
        )
      }
      , align='llrrr'
    )
    
    
    
    
    output$table2 <- renderTable(
      if(method() == "MidNormal"){
        setDT(
          data.frame(scenario=1:3, 
                     existing_approval = rep(fperc(a0()),3), 
                     new_approval = c(fperc(a1()), fperc(a2()), fperc(a3())),
                     existing_bad_rate = rep(fperc(b0()),3),
                     new_bad_rate = c(fperc(b1()), fperc(b2()), fperc(b3())),
                     existing_profit = rep(f3(profit0()), 3),
                     new_profit = c(f3(profit1()), f3(profit2()), f3(profit3()))
          )
          
        )
      }
      else if(method() == "MidFractal") {
        setDT(
          data.frame(scenario=1:3, 
                     existing_approval = rep(fperc(a0()),3), 
                     new_approval = c(fperc(a1_mf()), fperc(a2_mf()), fperc(a3_mf())),
                     existing_bad_rate = rep(fperc(b0_mf()),3),
                     new_bad_rate = c(fperc(b1_mf()), fperc(b2_mf()), fperc(b3_mf())),
                     existing_profit = rep(f3(profit0_mf()), 3),
                     new_profit = c(f3(profit1_mf()), f3(profit2_mf()), f3(profit3_mf()))
          )
        )
      }
      else if(method() == "BiFractal") {
        setDT(
        data.frame(scenario=1:3, 
                     existing_approval = rep(fperc(a0()),3), 
                     new_approval = c(fperc(a1_bf()), fperc(a2_bf()), fperc(a3_bf())),
                     existing_bad_rate = rep(fperc(b0_bf()),3),
                     new_bad_rate = c(fperc(b1_bf()), fperc(b2_bf()), fperc(b3_bf())),
                     existing_profit = rep(f3(profit0_bf()), 3),
                     new_profit = c(f3(profit1_bf()), f3(profit2_bf()), f3(profit3_bf()))
          )
        )
      }
      else if(method() == "BiNormal") {
        setDT(
        data.frame(scenario=1:3, 
                   existing_approval = rep(fperc(a0()),3), 
                   new_approval = c(fperc(a1_bn()), fperc(a2_bn()), fperc(a3_bn())),
                   existing_bad_rate = rep(fperc(b0_bn()),3),
                   new_bad_rate = c(fperc(b1_bn()), fperc(b2_bn()), fperc(b3_bn())),
                   existing_profit = rep(f3(profit0_bn()), 3),
                   new_profit = c(f3(profit1_bn()), f3(profit2_bn()), f3(profit3_bn()))
        )    
        )
      }, align='rrrrrrr'
     )
    
    output$table3 <- renderTable(
      if(method() == "MidNormal"){
        setDT(
          data.frame(scenario=1:3, 
                     existing_portfolio_gini = rep(f3(ginip0()/100),3),
                     new_portfolio_gini = c(f3(ginip1()/100), f3(ginip2()/100), f3(ginip3()/100)),
                     existing_marginal_bad_rate = rep(fperc(mbr0()),3),
                     new_marginal_bad_rate = c(fperc(mbr1()), fperc(mbr2()), fperc(mbr3()))
          )
          
        )
      }
      else if(method() == "MidFractal") {
        setDT(
          data.frame(scenario=1:3, 
                     existing_portfolio_gini = rep(f3(ginip0_mf()/100),3),
                     new_portfolio_gini = c(f3(ginip1_mf()/100), f3(ginip2_mf()/100), f3(ginip3_mf()/100)),
                     existing_marginal_bad_rate = rep(fperc(mbr0_mf()),3),
                     new_marginal_bad_rate = c(fperc(mbr1_mf()), fperc(mbr2_mf()), fperc(mbr3_mf()))
          )
        )
      }
      
      else if(method() == "BiFractal") {
        setDT(
          data.frame(scenario=1:3, 
                     existing_portfolio_gini = rep(f3(ginip0_bf()/100),3),
                     new_portfolio_gini = c(f3(ginip1_bf()/100), f3(ginip2_bf()/100), f3(ginip3_bf()/100)),
                     existing_marginal_bad_rate = rep(fperc(mbr0_bf()),3),
                     new_marginal_bad_rate = c(fperc(mbr1_bf()), fperc(mbr2_bf()), fperc(mbr3_bf()))
          )
        )
      }
      else if(method() == "BiNormal") {
        setDT(
          data.frame(scenario=1:3, 
                     existing_portfolio_gini = rep(f3(ginip0_bn()/100),3),
                     new_portfolio_gini = c(f3(ginip1_bn()/100), f3(ginip2_bn()/100), f3(ginip3_bn()/100)),
                     existing_marginal_bad_rate = rep(fperc(mbr0_bn()),3),
                     new_marginal_bad_rate = c(fperc(mbr1_bn()), fperc(mbr2_bn()), fperc(mbr3_bn()))
          ) 
        )
      }, align='rrrrr'
    )
    
    output$text1 <- renderText({
        # additional output
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

