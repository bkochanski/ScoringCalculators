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

    # setting functions with GINI1/SHAPE1 and GINI2/SHAPE2 parameters
    if(method() == "MidNormal"){
      reactive({y0<-function(x){FuncMidNormal(x,GINI1())}})
      reactive({y1<-function(x){FuncMidNormal(x,GINI2())}})
    }
    # else if(method() == "MidFractal"){
    #   y0<-function(x){FuncMidFractal(x,GINI1())}
    #   y1<-function(x){FuncMidFractal(x,GINI2())}
    # }
    # else if(method() == "BiNormal"){
    #   y0<-function(x){FuncBiNormal(x,GINI1(), shape1())}
    #   y1<-function(x){FuncBiNormal(x,GINI2(), ifelse(sameshape(), shape1(), shape2()))}
    # }
    # else if(method() == "BiFractal"){
    #   y0<-function(x){FuncBiFractal(x,GINI1(), beta1())}
    #   y1<-function(x){FuncBiFractal(x,GINI2(), ifelse(samebeta(), beta1(), beta2()))}
    # }
    
    # searching for representation of the cutoff point 
    phi0<-function(x){(1-B())*(1-x)+B()*(1-y0()(x))-a0()}
    x0<-reactive({as.numeric(uniroot(phi0,lower=0,upper=1,tol = .Machine$double.eps)[1])})

    # determining portfolio bad rate for scorecard 1
    b0 <- reactive({B()*(1-y0()(x0()))/a0()})
    
    # determining marginal bad rate at the cutoff for scorecard 1 
    deriv0 <-  reactive({numDeriv::grad(y0(), x0())})
    mbr0 <- reactive({(1+(1-B())/B()/deriv0())^(-1)})
    
    # determining interest rate assuming profit = 0 at marginal bad rate 
    ir0 <-  reactive({mbr0()/(1-mbr0())})
    
    # calculation of the portfolio profit
    profit0 <- reactive({a0()*(ir0()*(1-b0())-b0())})
    
    # Gini for scorecard 1 on accepted portfolio
    ginip0 <- reactive({GiniP(y0(),x0())})
    
    ## 1. Bad rate reduction scenario (keep approval)
    
    # same approval
    a1 <- reactive({a0()})
    
    # finding the cutoff representation on ROC curve 2 for scenario 1
    phi1<-function(x){(1-B())*(1-x)+B()*(1-y1()(x))-a1()}
    x1<-reactive({as.numeric(uniroot(phi1,lower=0,upper=1,tol = .Machine$double.eps)[1])})
    
    # new portfolio bad rate
    b1 <- reactive({B()*(1-y1()(x1()))/a0()})
    
    # new marginal bad rate at the cutoff for scenario 1
    deriv1 <- reactive({numDeriv::grad(y1(), x1())})
    mbr1 <- reactive({(1+(1-B())/B()/deriv1())^(-1)})
    
    # same interest rate
    ir1 <- reactive({ir0()})
    
    # profit for scorecard 2
    profit1 <- reactive({a1()*(ir1()*(1-b1())-b1())})
    
    # calculating Gini on approved portfolio for scorecard 2 
    ginip1 <- reactive({GiniP(y1(),x1())})
    
    # approval change, portfolio bad rate change, profit change
    a1_change <- reactive({a1()/a0()-1})
    b1_change <- reactive({b1()/b0()-1})
    profit1_change <- reactive({profit1()/profit0()-1})
    
    ## 2. Approval rate increase scenario (keep bad rate)

    # same portfolio bad rate
    b2 <- reactive({b0()})
    
    # finding the cutoff representation on ROC curve 2 for scenario 2
    phi2<-function(x){B()*(1-y1()(x))/((1-B())*(1-x)+B()*(1-y1()(x)))-b2()}
    x2 <- reactive({as.numeric(uniroot(phi2,lower=0.001,upper=.999,tol = .Machine$double.eps)[1])})
    
    # new approval rate
    a2 <- reactive({((1-B())*(1-x2())+B()*(1-y1()(x2())))})
    
    # as above
    deriv2 <- reactive({numDeriv::grad(y1(), x2())})
    mbr2 <- reactive({(1+(1-B())/B()/deriv2())^(-1)})
    ir2 <- reactive({ir0()})
    profit2 <- reactive({a2()*(ir2()*(1-b2())-b2())})
    ginip2 <- reactive({GiniP(y1(), x2())})
    
    a2_change <- reactive({a2()/a0()-1})
    b2_change <- reactive({b2()/b0()-1})
    profit2_change <- reactive({profit2()/profit0()-1})
    
    ## 3. Profit increase scenario (keep marginal bad rate)
    
    # same marginal bad rate - finding the cutoff representation on ROC curve 2 for scenario 3 
    phi3<-function(x){numDeriv::grad(y1(), x)-deriv0()}
    x3 <- reactive({as.numeric(uniroot(phi3,lower=0+1/1000,upper=1-1/1000,tol = .Machine$double.eps)[1])})
    
    # new approval rate and portfolio bad rate in scenario 3
    a3 <- reactive({((1-B())*(1-x3())+B()*(1-y1()(x3())))})
    b3 <- reactive({B()*(1-y1()(x3()))/a3()})
    
    # as above
    deriv3 <- reactive({numDeriv::grad(y1(), x3())})
    mbr3 <- reactive({(1+(1-B())/B()/deriv3())^(-1)})
    ir3 <- reactive({ir0()})
    profit3 <- reactive({a3()*(ir3()*(1-b3())-b3())})
    ginip3 <- reactive({GiniP(y1(),x3())})
    
    a3_change <- reactive({a3()/a0()-1})
    b3_change <- reactive({b3()/b0()-1})
    profit3_change <- reactive({profit3()/profit0()-1})

    ########
    
    output$curvePlot <- renderPlot({
            xcord <- c(x0(), x1(), x2(), x3())
            ycord <- c(y0()(x0()), y1()(x1()), y1()(x2()), y1()(x3()))
            b1_change <- b1_change()
            a2_change <- a2_change()
            profit3_change <- profit3_change()
        point_id <- c("a", "b", "c", "d")
        df <- data.frame(xcord, ycord, point_id)
        df$labels <- c('', paste("[1] Bad rate change:\n", ifelse(b1_change>0,"+",""), round(100*b1_change, 3), "%", sep = ""),
                       paste("[2] Approval change:\n", ifelse(a2_change>0,"+",""), round(100*a2_change, 3), "%", sep = ""),
                       paste("[3] Profit change:\n", ifelse(profit3_change>0,"+",""), round(100*profit3_change, 3), "%", sep = ""))
        
        
        ggplot(data.frame(x=c(0,1)), aes(x)) +
            stat_function(fun = y0(), geom = "line", aes(colour = "y0"), lwd = 1.1) +
            stat_function(fun = y1(), geom = "line", aes(colour = "y1"), linetype= "dashed", lwd = 1.1) +
            scale_colour_manual("ROC Curve", values = c("dodgerblue2", "olivedrab"), breaks = waiver(), labels=c("Scoring 1", "Scoring 2")) +
            geom_point(data = df[1,], aes(xcord, ycord, shape = factor(point_id)),  colour = "black", fill = "dodgerblue2", size = 4, stroke = 1.5) +
            geom_point(data = df[c(2,3, 4),], aes(xcord, ycord, shape = factor(point_id)), colour = "black", fill = "olivedrab", size = 4, stroke = 1.5) +
            scale_shape_manual("Point",
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
    
    #formats
    f3<-function(x){format(round(100*x,3), nsmall=3)}
    fperc<-function(x){paste(f3(x),"%", sep="")}

    output$table1 <- renderTable(
      setDT(
      data.frame(no=1:3, scenario=c('keep approval (bad rate reduction)', 'keep bad rate (approval rate improvement)', 'keep marginal bad rate (profit increase)'), 
                 approval_change=c(fperc(a1_change()), fperc(a2_change()), fperc(a3_change())),
                 bad_rate_change=c(fperc(b1_change()), fperc(b2_change()), fperc(b3_change())),
                 profit_change=c(fperc(profit1_change()), fperc(profit2_change()), fperc(profit3_change()))
                 )
      ) 
      , align='llrrr'
    )

    output$table2 <- renderTable(
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
      , align='rrrrrrr'
     )
    
    output$table3 <- renderTable(
        setDT(
          data.frame(scenario=1:3, 
                     existing_portfolio_gini = rep(f3(ginip0()/100),3),
                     new_portfolio_gini = c(f3(ginip1()/100), f3(ginip2()/100), f3(ginip3()/100)),
                     existing_marginal_bad_rate = rep(fperc(mbr0()),3),
                     new_marginal_bad_rate = c(fperc(mbr1()), fperc(mbr2()), fperc(mbr3()))
          )
        )
      , align='rrrrr'
    )
    
    output$text1 <- renderText({
        # additional output
    })
}

# Run the application 
shinyApp(ui = ui, server = server)