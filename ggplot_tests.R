library(shiny)
library(shinyWidgets)
library(ggplot2)
library(ggrepel)
library(numDeriv)
library(data.table)

# FUNCTIONS:
FuncMidNormal<-function(x,g){pnorm(qnorm((g+1)/2)*sqrt(2)+qnorm(x))}
FuncMidFractal<-function(x,g){0.5*(1-(1-x)^((1+g)/(1-g)))+0.5*x^((1-g)/(1+g))}

FuncBiFractal<-function(x,g, beta){beta*(1-(1-x)^((1+g)/(1-g)))+(1-beta)*x^((1-g)/(1+g))}
FuncBiNormal<-function(x,g,shape){pnorm(qnorm((g+1)/2)*sqrt(1+shape^2)+shape*qnorm(x))}

GiniP<-function(f, x){2*(integrate(f,x,1)$value-(1-x)*f(x))/((1-x)*(1-f(x)))-1}

  GINI1<-.45
  GINI2<-.65
  B <- .10
  a0 <- .60
    
  # MID NORMAL
  y0<-function(x){FuncMidNormal(x,GINI1)}
  y1<-function(x){FuncMidNormal(x,GINI2)}
  
  phi0<-function(x){(1-B)*(1-x)+B*(1-y0(x))-a0}
  x0<-as.numeric(uniroot(phi0,lower=0,upper=1,tol = .Machine$double.eps)[1])
  b0 <- B*(1-y0(x0))/a0
  
  deriv0 <-  numDeriv::grad(y0, x0)
  mbr0 <- (1+(1-B)/B/deriv0)^(-1)
  ir0 <-  mbr0/(1-mbr0)
  profit0 <- a0*(ir0*(1-b0)-b0)
  ginip0 <- GiniP(y0,x0)
  
  
  a1 <- a0
  phi1<-function(x){(1-B)*(1-x)+B*(1-y1(x))-a1}
  x1<-as.numeric(uniroot(phi1,lower=0,upper=1,tol = .Machine$double.eps)[1])
  b1 <- B*(1-y1(x1))/a0
  
  deriv1 <- numDeriv::grad(y1, x1)
  mbr1 <- (1+(1-B)/B/deriv1)^(-1)
  ir1 <- ir0
  profit1 <- a1*(ir1*(1-b1)-b1)
  ginip1 <- GiniP(y1,x1)
  a1_change <- a1/a0-1
  b1_change <- b1/b0-1
  profit1_change <- profit1/profit0-1
  
  b2 <- b0
  phi2<-function(x){B*(1-y1(x))/((1-B)*(1-x)+B*(1-y1(x)))-b2}
  x2 <- as.numeric(uniroot(phi2,lower=0.001,upper=.999,tol = .Machine$double.eps)[1])
  a2 <- ((1-B)*(1-x2)+B*(1-y1(x2)))
  
  deriv2 <- numDeriv::grad(y1, x2)
  mbr2 <- (1+(1-B)/B/deriv2)^(-1)
  ir2 <- ir0
  profit2 <- a2*(ir2*(1-b2)-b2)
  ginip2 <- GiniP(y1, x2)
  
  a2_change <- a2/a0-1
  b2_change <- b2/b0-1
  profit2_change <- profit2/profit0-1
  
  
  ## NEW SCENARIO
  phi3<-function(x){numDeriv::grad(y1, x)-deriv0}
  x3 <- as.numeric(uniroot(phi3,lower=0+1/1000,upper=1-1/1000,tol = .Machine$double.eps)[1])
  a3 <- ((1-B)*(1-x3)+B*(1-y1(x3)))
  b3 <- B*(1-y1(x3))/a3
  
  
  deriv3 <- numDeriv::grad(y1, x3)
  mbr3 <- (1+(1-B)/B/deriv3)^(-1)
  ir3 <- ir0
  profit3 <- a3*(ir3*(1-b3)-b3)
  ginip3 <- GiniP(y1,x3)
  
  a3_change <- a3/a0-1
  b3_change <- b3/b0-1
  profit3_change <- profit3/profit0-1
  

      xcord <- c(x0, x1, x2, x3)
      ycord <- c(y0(x0), y1(x1), y1(x2), y1(x3))
      b1_change <- b1_change
      a2_change <- a2_change
      profit3_change <- profit3_change
    point_id <- c("a", "b", "c", "d")
    df <- data.frame(xcord, ycord, point_id)
    df$labels <- c('', paste("[1] Bad rate change:\n", ifelse(b1_change>0,"+",""), round(100*b1_change, 3), "%", sep = ""),
                  paste("[2] Approval change:\n", ifelse(a2_change>0,"+",""), round(100*a2_change, 3), "%", sep = ""),
                  paste("[3] Profit change:\n", ifelse(profit3_change>0,"+",""), round(100*profit3_change, 3), "%", sep = ""))
    
    ggplot(data.frame(x=c(0,1)), aes(x)) +
      stat_function(fun = y0, geom = "line", aes(colour = "y0"), lwd = 1.1) +
      stat_function(fun = y1, geom = "line", aes(colour = "y1"), linetype= "dashed", lwd = 1.1) +
      scale_colour_manual("ROC Curve", values = c("dodgerblue2", "olivedrab"), breaks = waiver(), labels=c("Scoring 1", "Scoring 2")) +
      geom_point(data = df[1,], aes(xcord, ycord, shape = factor(point_id)),  colour = "black", fill = "dodgerblue2", size = 4, stroke = 1.5) +
      geom_point(data = df[c(2,3, 4),], aes(xcord, ycord, shape = factor(point_id)), colour = "black", fill = "olivedrab", size = 4, stroke = 1.5) +
      #            geom_point(data = df[1:4], aes(xcord, ycord, shape = factor(point_id)),  colour = "black", fill = c("dodgerblue2", "olivedrab", "olivedrab", "olivedrab"), size = 4, stroke = 1.5) +
      scale_shape_manual("Point",
                         values = c('a' = 21, 'b' = 22, 'c' = 24, 'd' = 23), 
                         labels = c("Currently", "[1] Bad rate reduction\n scenario (keep approval)", "[2] Approval rate improvement\n scenario (keep bad rate)", "[3] Profit increase\n scenario (keep marginal\n bad rate)")) +
      xlab("cumulative good proportion") + ylab("cumulative bad proportion") + ggtitle("") + 
      theme(legend.position="bottom", legend.box = "vertical", axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold")) +
      guides(color = guide_legend(order = 1)) +
      geom_text_repel(data = df[2:4,],
                mapping = aes(x=xcord, y=ycord, label = labels),
                size=5,
                min.segment.length = 0,
                hjust = 1, vjust = -1, show.legend = FALSE)
    
    
#      geom_text_repel(x=xcord[2], y=ycord[2], label = paste("[1] Bad rate change:\n", round(100*b1_change, 3), "%", sep = " "), hjust = 1, vjust = -1, aes(size = 2), show.legend = FALSE, label.padding=0 ) +
#      geom_text(x=xcord[3], y=ycord[3], label = paste("[2] Approval change:\n", round(100*a2_change, 3), "%", sep = " "), hjust = 1, vjust = -1, aes(size = 2), show.legend = FALSE, label.padding=0 ) +
#      geom_text(x=xcord[4], y=ycord[4], label = paste("[3] Profit change:\n", round(100*profit3_change, 3), "%", sep = " "), hjust = 1, vjust = -1, aes(size = 2), show.legend = FALSE, label.padding=0 )
    # # code responsible for the square shape of the plot
  }
  # ,height=reactive(ifelse(!is.null(input$innerWidth),input$innerWidth*2/5,0))
  
  )
  