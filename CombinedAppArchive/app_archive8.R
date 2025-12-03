library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(ggrepel)

ui <- dashboardPage(
  
  # HEADER
  dashboardHeader(title="Combined Scorecards Calculator", titleWidth = 580),
  
  # SIDEBAR
  dashboardSidebar(
    width = 350,
    h4("Input Parameters", style = "text-align: center; padding: 10px; font-size: 22px; font-weight: bold;"),
    
    div(style = "padding: 0 15px; margin-bottom: 30px;",
        tags$style(HTML("
                .irs-bar {height: 14px;}
                .irs-line {height: 14px;}
                .irs-slider {width: 26px; height: 26px; top: 19px;}
                .irs-from, .irs-to, .irs-single {font-size: 16px;}
                .irs-grid-text {font-size: 14px;}
            ")),
        sliderInput("g1",
                    "GINI 1:",
                    min = 0,
                    max = 1,
                    value = .4,
                    step = .01,
                    width = "65%"),
        numericInput("g1_num",
                     NULL,
                     min = 0,
                     max = 1,
                     value = .4,
                     step = .01,
                     width = "200px")
    ),
    
    div(style = "padding: 0 15px; margin-bottom: 30px;",
        sliderInput("g2",
                    "GINI 2:",
                    min = 0,
                    max = 1,
                    value = .3,
                    step = .01,
                    width = "65%"),
        numericInput("g2_num",
                     NULL,
                     min = 0,
                     max = 1,
                     value = .3,
                     step = .01,
                     width = "200px")
    ),
    
    div(style = "padding: 0 15px; margin-bottom: 30px;",
        sliderInput("corr",
                    "Correlation:",
                    min = 0,
                    max = 1,
                    value = .15,
                    step = .01,
                    width = "65%"),
        numericInput("corr_num",
                     NULL,
                     min = 0,
                     max = 1,
                     value = .15,
                     step = .01,
                     width = "200px")
    ),
    
    div(style = "padding: 0 15px; margin-bottom: 30px;",
        sliderInput("defrate",
                    "Default rate:",
                    min = 0,
                    max = 1,
                    value = .1,
                    step = .01,
                    width = "65%"),
        numericInput("defrate_num",
                     NULL,
                     min = 0,
                     max = 1,
                     value = .1,
                     step = .01,
                     width = "200px")
    ),
    
    tags$style(HTML("
            .sidebar label {font-size: 19px; font-weight: 700;}
            .sidebar input[type='number'] {font-size: 19px; height: 48px;}
        "))
  ),
  
  # BODY
  dashboardBody(
    tags$style(HTML("
            .small-box {height: 65px;}
            .small-box .icon-large {font-size: 40px;}
            .small-box h3 {font-size: 22px; margin: 3px 0;}
            .small-box p {font-size: 11px; margin: 0;}
            .info-box {height: 60px; min-height: 60px;}
            .info-box-icon {height: 60px; line-height: 60px; font-size: 24px;}
            .info-box-content {padding: 5px 8px;}
            .info-box-text {font-size: 11px;}
            .info-box-number {font-size: 18px; font-weight: bold;}
        ")),
    
    fluidRow(
      box(
        title = "GINI Combination Chart",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        plotOutput("gini_chart", height = "700px")
      )
    ),
    
    fluidRow(
      valueBoxOutput("box_new_gini", width = 12)            
    ),
    
    fluidRow(
      infoBoxOutput("box_weight1", width = 6),
      infoBoxOutput("box_weight2", width = 6)
    ),
    
    fluidRow(
      valueBoxOutput("box_new_rho", width = 12)
    ),
    
    fluidRow(
      infoBoxOutput("box_rho1", width = 6),
      infoBoxOutput("box_rho2", width = 6)
    ),
    
    fluidRow(
      valueBoxOutput("box_a_opt", width = 12)
    )
  )
)

# CALCULATOR LOGIC
gini_combine_calculator <- function(g1, g2, corr, defaultrate){
  ginic <- function(bc, gc){
    sum((gc[2:(length(gc))]-gc[1:(length(gc)-1)])*
          (bc[2:(length(bc))]+bc[1:(length(bc)-1)]))-1
  } 
  
  gini_crd <- function(rho=0.5, defrate=0.1, gran=10000) {
    drates_i <- pnorm(qnorm(defrate), rho*(qnorm(1:(gran-1)/gran)),sqrt(1-rho^2))
    drates_2i <- (c(1, drates_i)+c(drates_i, 0))/2
    return(ginic(cumsum((drates_2i)/sum(drates_2i)), cumsum((1-drates_2i)/sum(1-drates_2i))))
  }
  
  # rho_s1
  phi_s1 <- function(x){gini_crd(rho=x, defrate=defaultrate, gran=100000)-g1}
  rho_s1 <- uniroot(phi_s1,lower=0,upper=1,tol = .Machine$double.eps)$root
  
  # rho_s2
  phi_s2 <- function(x){gini_crd(rho=x, defrate=defaultrate, gran=100000)-g2}
  rho_s2 <- uniroot(phi_s2,lower=0,upper=1,tol = .Machine$double.eps)$root
  
  a_opt <- (corr*rho_s2-rho_s1)/(corr*rho_s1-rho_s2)
  corr_opt <- (a_opt*rho_s1+rho_s2)/sqrt(a_opt^2+2*a_opt*corr+1)
  g_result0 <- gini_crd(corr_opt, defaultrate, gran=100000)
  g_result1 <- if(a_opt<0 | a_opt>1000) {NaN} else {g_result0}
  
  return(c(new_gini=g_result1, 
           a_opt=a_opt, 
           score_1_weight=a_opt/(1+a_opt), 
           score_2_weight=1/(1+a_opt), 
           rho1=rho_s1, 
           rho2=rho_s2, 
           new_corr=corr_opt,
           gini1=g1,
           gini2=g2))
}

# Define server logic
server <- function(input, output, session) {
  # Track which input was last changed
  lastChanged <- reactiveValues(
    g1 = "slider",
    g2 = "slider", 
    corr = "slider",
    defrate = "slider"
  )
  
  # Input values converted into reactive values
  G1 <- reactive({
    if (lastChanged$g1 == "number") input$g1_num else input$g1
  })
  
  G2 <- reactive({
    if (lastChanged$g2 == "number") input$g2_num else input$g2
  })
  
  CORR <- reactive({
    if (lastChanged$corr == "number") input$corr_num else input$corr
  })
  
  DR <- reactive({
    if (lastChanged$defrate == "number") input$defrate_num else input$defrate
  })
  
  # Synchronize sliders and numeric inputs
  observeEvent(input$g1, {
    lastChanged$g1 <- "slider"
    updateNumericInput(session, "g1_num", value = input$g1)
  }, ignoreInit = TRUE)
  
  observeEvent(input$g1_num, {
    lastChanged$g1 <- "number"
    updateSliderInput(session, "g1", value = input$g1_num)
  }, ignoreInit = TRUE)
  
  observeEvent(input$g2, {
    lastChanged$g2 <- "slider"
    updateNumericInput(session, "g2_num", value = input$g2)
  }, ignoreInit = TRUE)
  
  observeEvent(input$g2_num, {
    lastChanged$g2 <- "number"
    updateSliderInput(session, "g2", value = input$g2_num)
  }, ignoreInit = TRUE)
  
  observeEvent(input$corr, {
    lastChanged$corr <- "slider"
    updateNumericInput(session, "corr_num", value = input$corr)
  }, ignoreInit = TRUE)
  
  observeEvent(input$corr_num, {
    lastChanged$corr <- "number"
    updateSliderInput(session, "corr", value = input$corr_num)
  }, ignoreInit = TRUE)
  
  observeEvent(input$defrate, {
    lastChanged$defrate <- "slider"
    updateNumericInput(session, "defrate_num", value = input$defrate)
  }, ignoreInit = TRUE)
  
  observeEvent(input$defrate_num, {
    lastChanged$defrate <- "number"
    updateSliderInput(session, "defrate", value = input$defrate_num)
  }, ignoreInit = TRUE)
  
  # Vector of values returned by the calculator function
  resultVector <- reactive({
    gini_combine_calculator(G1(), G2(), CORR(), DR())
  })
  
  # Chart output
  output$gini_chart <- renderPlot({
    g <- resultVector()
    corr <- CORR()
    defrate <- DR()
    
    # Helper function
    w2a <- function(w){w/(1-w)}
    ws <- 1:199/200
    
    # Gini from rho function
    gini_crd <- function(rho=0.5, defrate=0.1, gran=10000) {
      ginic <- function(bc, gc){
        sum((gc[2:(length(gc))]-gc[1:(length(gc)-1)])*
              (bc[2:(length(bc))]+bc[1:(length(bc)-1)]))-1
      }
      drates_i <- pnorm(qnorm(defrate), rho*(qnorm(1:(gran-1)/gran)),sqrt(1-rho^2))
      drates_2i <- (c(1, drates_i)+c(drates_i, 0))/2
      return(ginic(cumsum((drates_2i)/sum(drates_2i)), cumsum((1-drates_2i)/sum(1-drates_2i))))
    }
    
    gini_from_aop <- function(a, corr_val){
      rho_combined <- (a * g["rho1"] + g["rho2"]) / sqrt(a^2 + 2 * a * corr_val + 1)
      gini_crd(rho_combined, defrate, gran=10000)
    }
    gini_from_aop <- Vectorize(gini_from_aop)
    
    # Create data frame for the curve
    gdf <- data.frame(
      ws = ws, 
      corr1 = gini_from_aop(w2a(ws), corr_val = corr)
    )
    
    # Create data frame for optimal point
    gresdf <- data.frame(
      score_1_weight = g["score_1_weight"],
      new_gini = g["new_gini"],
      y30 = as.numeric(gdf[which.min(abs(gdf$ws - g["score_2_weight"])), 2])
    )
    
    # Create the plot
    p <- ggplot() + 
      theme_bw(base_size = 20) + 
      xlab("Normalized weight of scorecard 1") + 
      ylab("Gini of the combined scorecard") + 
      scale_x_continuous(breaks=0:5/5) +
      theme(
        axis.title = element_text(size = 24, face = "bold"),
        axis.text = element_text(size = 18),
        plot.margin = margin(25, 25, 25, 25)
      ) +
      geom_line(data=gdf, aes(x=ws, y=corr1), color="steelblue", size=2) +
      geom_point(data=gresdf, aes(x=score_1_weight, y=new_gini), 
                 color="red", size=7) +
      geom_text_repel(data=gresdf, 
                      aes(x=score_1_weight, y=new_gini, 
                          label=round(new_gini, 3)),
                      color="red", fontface="bold", size=8,
                      box.padding = 1.5) +
      geom_label(data=gresdf, 
                 aes(x=score_1_weight, y=y30, 
                     label=paste0("corr = ", corr)),
                 nudge_y = -0.01, size=6.5) +
      geom_text_repel(data=data.frame(x=-0.05, y=g["gini2"], 
                                      la=paste0("Gini of scorecard 2 = ", round(g["gini2"], 3))),
                      aes(x=x, y=y, label=la), size=6.5,
                      box.padding = 1.5) +
      geom_text_repel(data=data.frame(x=1.05, y=g["gini1"], 
                                      la=paste0("Gini of scorecard 1 = ", round(g["gini1"], 3))),
                      aes(x=x, y=y, label=la), size=6.5,
                      box.padding = 1.5)
    
    print(p)
  })
  
  # BOXES
  output$box_new_gini <- renderValueBox({
    valueBox("GINI of combined models", 
             tags$p(round(resultVector()[1], 4), style = 'font-size: 120%; font-weight: bold;'), 
             icon=icon("compress-arrows-alt"), 
             color = "light-blue")
  })
  
  output$box_weight1 <- renderInfoBox({
    infoBox("Scoring 1 weight", 
            round(resultVector()[3], 4), 
            icon=icon("balance-scale-left"), 
            color = "aqua")
  })
  
  output$box_weight2 <- renderInfoBox({
    infoBox("Scoring 2 weight", 
            round(resultVector()[4], 4), 
            icon=icon("balance-scale-right"), 
            color = "aqua")
  })
  
  output$box_new_rho <- renderValueBox({
    valueBox("New rho", 
             tags$p(round(resultVector()[7], 4), style = 'font-size: 120%; font-weight: bold;'), 
             icon=icon("chart-line"), 
             color = "olive")
  })
  
  output$box_rho1 <- renderInfoBox({
    infoBox("Rho 1", 
            round(resultVector()[5], 4), 
            icon=icon("chart-bar"), 
            color = "green")
  })
  
  output$box_rho2 <- renderInfoBox({
    infoBox("Rho 2", 
            round(resultVector()[6], 4), 
            icon=icon("chart-bar"), 
            color = "green")
  })
  
  output$box_a_opt <- renderValueBox({
    valueBox("A opt", 
             tags$p(round(resultVector()[2], 4), style = 'font-size: 120%; font-weight: bold;'), 
             icon=icon("cloudscale"), 
             color = "yellow")
  })
}

# Run the application
shinyApp(ui = ui, server = server)