library(shiny)
library(shinydashboard)

ui <- dashboardPage(
    
    # HEADER
    dashboardHeader(title="Combined Scorecards Calculator", titleWidth = 580),
    
    # SIDEBAR
    dashboardSidebar(
        numericInput("g1",
                     "GINI 1:",
                     min = 0,
                     max = 1,
                     value = .6,
                     step = .01),
        numericInput("g2",
                     "GINI 2:",
                     min = 0,
                     max = 1,
                     value = .6,
                     step = .01),
        numericInput("corr",
                     "Correlation:",
                     min = 0,
                     max = 1,
                     value = .6,
                     step = .01),
        numericInput("defrate",
                     "Default rate:",
                     min = 0,
                     max = 1,
                     value = .1,
                     step = .01)
        
    ),
    
    # BODY
    dashboardBody(
        fluidRow(
            valueBoxOutput("box_new_gini", width = 12)            
            
        ),
        
        fluidRow(
            infoBoxOutput("box_weight1", width = 6),
            infoBoxOutput("box_weight2", width = 6)
        ),
        br(),br(),br(),
        
        fluidRow(
            valueBoxOutput("box_new_rho", width = 12)
        ),
        
        fluidRow(
            infoBoxOutput("box_rho1", width = 6),
            infoBoxOutput("box_rho2", width = 6)
        ),
        br(),br(),br(),
        

        fluidRow(
            valueBoxOutput("box_a_opt", width = 12)
        )
    )
)

# CALCULATOR LOGIC
gini_combine_calculator<-function(g1, g2, corr, defaultrate){
    ginic<-function(bc, gc){
        #function for gini when we have cumulative goods and bads vectors
        sum((gc[2:(length(gc))]-gc[1:(length(gc)-1)])*
                (bc[2:(length(bc))]+bc[1:(length(bc)-1)]))-1
    } 
    
    gini_from_r<-function(rho=.5, defrate=.1){
        F1<-function(s, d, rho){integrate(function(x){pnorm((qnorm(d)-rho*x)/sqrt(1-rho^2))*dnorm(x)}, lower=-Inf, upper=s)$value}
        F2<-function(s,d,rho){pnorm((-qnorm(d)+rho*s)/sqrt(1-rho^2))*dnorm(s)}
        F1_<-Vectorize(function(x){F1(x, defrate, rho)})
        F2_<-Vectorize(function(x){F2(x, defrate, rho)})
        2*integrate(function(x){F1_(x)*F2_(x)}, 
                    lower=-Inf, upper=Inf, subdivisions=200)$value/defrate/(1-defrate)-1
    }
    
    #rho_s1
    phi_s1<-function(x){gini_from_r(rho=x, defrate=defaultrate)-g1}
    rho_s1<-uniroot(phi_s1,lower=0,upper=1,tol = .Machine$double.eps)$root
    #rho_s2
    phi_s2<-function(x){gini_from_r(rho=x, defrate=defaultrate)-g2}
    rho_s2<-uniroot(phi_s2,lower=0,upper=1,tol = .Machine$double.eps)$root
    
    (a_opt<-(corr*rho_s2-rho_s1)/(corr*rho_s1-rho_s2))
    corr_opt<-(a_opt*rho_s1+rho_s2)/sqrt(a_opt^2+2*a_opt*corr+1)
    g_result0<-gini_from_r(corr_opt, defaultrate)
    g_result1<-if(a_opt<0 | a_opt>1000) {NaN} else {g_result0}
    return(c(new_gini=g_result1, 
             #new_gini_2=g_result0, 
             a_opt=a_opt, score_1_weight=a_opt/(1+a_opt), score_2_weight=1/(1+a_opt), 
             rho1=rho_s1, rho2=rho_s2, new_corr=corr_opt))
}


# Define server logic required to draw a histogram
server <- function(input, output) {
    # inout values cinerted into reactive values
    G1 <- reactive({
        input$g1
    })
    G2 <- reactive({
        input$g2
    })
    CORR <- reactive({
        input$corr
    })
    DR <- reactive({
        input$defrate
    })
    
    # vector of 7 values returned by the calculator function
    resultVector <- reactive({gini_combine_calculator(G1(),G2(),CORR(),DR())})
    
    # NOT USED outputs - can by displayed anytime
    # main output:
    output$new_gini <- renderText({
        paste(resultVector()[1], sep='')
    })
    output$weight1 <- renderText({
        paste(resultVector()[3], sep='')
    })
    output$weight2 <- renderText({
        paste(resultVector()[4], sep='')
    })
    
    # auxillary output:
    output$a_opt <- renderText({
        paste(resultVector()[2], sep='')
    })
    output$rho1 <- renderText({
        paste(resultVector()[5], sep='')
    })
    output$rho2 <- renderText({
        paste(resultVector()[6], sep='')
    })
    output$rho_new <- renderText({
        paste(resultVector()[7], sep='')
    })
    
    # USED outputs
    # BOXES
    output$box_new_gini <- renderValueBox({
        valueBox("GINI of combined models", tags$p(paste(resultVector()[1]), style = 'font-size: 200%; font-weight: bold;'), icon=icon("compress-arrows-alt"), color = "light-blue")
    })
    
    output$box_weight1 <- renderInfoBox({
        infoBox("Scoring 1 weight", paste(resultVector()[3]), icon=icon("balance-scale-left"), color = "aqua")
    })
    
    output$box_weight2 <- renderInfoBox({
        infoBox("Scoring 2 weight", paste(resultVector()[4]), icon=icon("balance-scale-right"), color = "aqua")
    })
    
    output$box_new_rho <- renderValueBox({
        valueBox("New rho", tags$p(paste(resultVector()[7]), style = 'font-size: 200%; font-weight: bold;'), icon=icon("chart-line"), color = "olive")
    })
    
    output$box_rho1 <- renderInfoBox({
        infoBox("Rho 1", paste(resultVector()[5]), icon=icon("chart-bar"), color = "green")
    })
    
    output$box_rho2 <- renderInfoBox({
        infoBox("Rho 2", paste(resultVector()[6]), icon=icon("chart-bar"), color = "green")
    })
    
    output$box_a_opt <- renderValueBox({
        valueBox("A opt", tags$p(paste(resultVector()[2]), style = 'font-size: 200%; font-weight: bold;'), icon=icon("cloudscale"), color = "yellow")
    })
    
    
}

# Run the application
shinyApp(ui = ui, server = server)

