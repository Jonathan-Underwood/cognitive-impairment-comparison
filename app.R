library(shiny)
require(ggplot2)
require(MASS)
require(cowplot)


#### user interface logic ####
ui <- fluidPage(
  
  # Application title
  titlePanel("Compare classification techniques"),
  
  # Sidebar with a slider input 
  sidebarLayout(
    sidebarPanel(
      h5("Jonathan Underwood 2016"),
      selectInput("select", label = h3("Select type of display"), 
                  choices = list("Histograms", "PC plots"), 
                  selected = "Histograms"),
      sliderInput("mean",
                  "Population mean",
                  min = 0,
                  max = 100,
                  value = 50),
      sliderInput("sd",
                  "Population standard deviation",
                  min = 0,
                  max = 50,
                  value = 10),
      sliderInput("imp",
                  "Percentage impaired (%)",
                  min = 0,
                  max = 50,
                  value = 10),
      sliderInput("imp_mean",
                  "Impaired mean",
                  min = 0,
                  max = 100,
                  value = 30),
      sliderInput("imp_sd",
                  "Impaired standard deviation",
                  min = 0,
                  max = 50,
                  value = 10),
      sliderInput("alpha",
                  "Mahalanobis distance alpha (%)",
                  min = 1,
                  max = 20,
                  value = 5)
    ),
    
    # Show a plot of the generated distributions
    mainPanel(
      plotOutput("plot"),
      hr(),
      tableOutput("table"),
      p(strong("Abbreviations:"), "PPV = positive predictive value; NPV = negative predictive value"),
      hr(),
      p("This models a simulated population (n=10,000) with with user defined parameters to test the assumptions of various definitions of cognitive impairment. It models testing cognitive function across six domains, which follow a multivariate normal distribution with covariance informed by the ",
        a("POPPY study.",
        href = "http://doi.org/10.1186/s12879-016-1970-8")),
      br(),
      p("A user defined population of 'impaired' individuals can also be added. The ", 
        a("Frascati",
          href = "http://doi.org/10.1212/01.WNL.0000287431.88658.8b"), 
        " and ", 
        a("global deficit score", 
          href = "http://doi.org/10.1080/13803390490510031"), 
        " criteria are then applied to the data assuming a normative mean T-score of 50 and standard deviation of 10. Additonally, a novel way of defining impairment using the multivariate statistic, the Mahalanobis distance, is also applied, again assuming a normative mean T-score of 50 and standard deviation of 10.")
    )
  )
)


#### server logic ####

server <- (function(input, output) {
  # POPPY cognitive correlation matrix
  syn_n=10000
  sigma <-matrix(c(1,0.210754443741029,0.145349600511535,0.528406352623827,0.367212674407502,0.377429548420058,
                   0.210754443741029,1,0.601043770264608,0.233912234063957,0.156447757547805,0.27289232693919,
                   0.145349600511535,0.601043770264608,1,0.0580374512388783,0.166080151414515,0.357473226464959,
                   0.528406352623827,0.233912234063957,0.0580374512388783,1,0.210596654179787,0.222291349676154,
                   0.367212674407502,0.156447757547805,0.166080151414515,0.210596654179787,1,0.100172102842648,
                   0.377429548420058,0.27289232693919,0.357473226464959,0.222291349676154,0.100172102842648,1), nrow = 6)
  
  ncol(sigma) -> p
  
  synthetic_pop <- reactive({
    ifelse(input$imp>0, {
      syn_n*(input$imp/100) -> n2
      syn_n-n2 ->n
      data.frame((mvrnorm(n, mu = rep(0,p), Sigma = sigma,empirical = TRUE) * input$sd) +input$mean) -> synthetic_pop2
      rep(1,n) -> synthetic_pop2$label
      
      data.frame((mvrnorm(n2, mu = rep(0,p), Sigma = sigma,empirical = TRUE) * input$imp_sd) +input$imp_mean) ->synthetic_pop2.imp
      rep(2,n2) -> synthetic_pop2.imp$label
      rbind(synthetic_pop2, synthetic_pop2.imp) ->synthetic_pop2
      factor(synthetic_pop2$label) -> synthetic_pop2$label
    }, {
      data.frame((mvrnorm(syn_n, mu = rep(0,p), Sigma = sigma,empirical = TRUE) * input$sd) +input$mean) -> synthetic_pop2
      rep(1,syn_n) -> synthetic_pop2$label
      factor(synthetic_pop2$label) -> synthetic_pop2$label
    })
    apply(synthetic_pop2[1:p],1,function(x) sum(x<40)>1) -> HAND    # generates frascati
    means <- rep(50,p) # assuming a normative mean of 50 (like Frascatia and GDS)
    covar <- sigma * 100 # assumes a normative sd of 10 for each domain (like Frascatia and GDS)
    apply(synthetic_pop2[1:p],1, function(i) sign(mean(i)-50)*sqrt(mahalanobis(i,center = means, cov = covar))) -> mahal_dist # generates signed mahal distance 
    
    apply(synthetic_pop2[1:p],2,function(x) 
      ifelse(x>19,
             ifelse(x>24,
                    ifelse(x>29,
                           ifelse(x>34,
                                  ifelse(x>40,0,1),2),3),4),5)
    ) -> test
    apply(test,1,function(x) mean(x)>=0.5) ->GDS    # generates global deficit score
    
    input$alpha/50 ->alpha
    -sqrt(qbeta(1-alpha,p/2, (syn_n-p-1)/2 ) * ((syn_n-1)^2)/syn_n) -> critical # calculates critical value for MD
    
    list(synthetic_pop2=synthetic_pop2, HAND = HAND, GDS = GDS, mahal_dist = mahal_dist, critical = critical)
  })
  
  output$plot <- renderPlot({
    
    # get the dynamic values
    synthetic_pop2 <- synthetic_pop()$synthetic_pop2
    HAND <- synthetic_pop()$HAND
    GDS <- synthetic_pop()$GDS
    mahal_dist <- synthetic_pop()$mahal_dist
    critical <- synthetic_pop()$critical 
    
    rowMeans(synthetic_pop2[1:p]) -> global_score
    
    ifelse(input$select == "Histograms", {
      
      ggplot(data = synthetic_pop2,aes(x = global_score, ..count..,fill= factor(synthetic_pop2$label))) + 
        geom_histogram(alpha=0.5, position = "identity", binwidth = 3, colour= "black") + 
        scale_fill_discrete(name= "Population",labels=c("Normal","Impaired")) +
        theme_classic() +
        theme(
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          legend.position=c(0.9,0.9))+
        labs(x="Global T-score", y= "Count") +
        scale_y_continuous(limits = c(0,2000), expand = c(0, 0))+
        ggtitle("Population") -> plot.synthetic_imp.true
      
      ggplot(data = synthetic_pop2,aes(x = global_score, ..count..,fill= factor(HAND))) + geom_histogram(alpha=0.5, position = "identity", binwidth = 3, colour= "black") + 
        scale_fill_discrete(name= "Frascati",labels=c("Normal","Impaired")) +
        theme_classic() +
        theme(
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          legend.position=c(0.9,0.9))+
        labs(x="Global T-score", y= "Count") +
        scale_y_continuous(limits = c(0,2000), expand = c(0, 0))+
        ggtitle("Frascati") -> plot.synthetic_imp.HAND
      
      ggplot(data = synthetic_pop2,aes(x = global_score, ..count..,fill= factor(GDS))) + geom_histogram(alpha=0.5, position = "identity", binwidth = 3, colour= "black") + 
        scale_fill_discrete(name= "GDS",labels=c("Normal","Impaired")) +
        theme_classic() +
        theme(
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          legend.position=c(0.9,0.9))+
        labs(x="Global T-score", y= "Count") +
        scale_y_continuous(limits = c(0,2000), expand = c(0, 0))+
        ggtitle("Global deficit score") -> plot.synthetic_imp.GDS
      
      ggplot(data = synthetic_pop2,aes(x = global_score, ..count..,fill= factor(mahal_dist<critical))) + geom_histogram(alpha=0.5, position = "identity", binwidth = 3, colour= "black") + 
        scale_fill_discrete(name= "MD",labels=c("Normal","Impaired")) +
        theme_classic() +
        theme(
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          legend.position=c(0.9,0.9))+
        labs(x="Global T-score", y= "Count") +
        scale_y_continuous(limits = c(0,2000), expand = c(0, 0))+
        ggtitle("Mahalanobis Distance") -> plot.synthetic_imp.mahal
      
      print(plot_grid(plot.synthetic_imp.true, plot.synthetic_imp.HAND, plot.synthetic_imp.GDS, plot.synthetic_imp.mahal, align = "hv"))
    }, {
      princomp(synthetic_pop2[1:p])$score[,1] -> PC.1
      princomp(synthetic_pop2[1:p])$score[,2] -> PC.2
      
      qplot(PC.1,PC.2, col = synthetic_pop2$label, alpha = I(0.1))+
        scale_colour_discrete(name= "Population",labels=c("Normal","Impaired")) +
        theme_classic() +
        theme(
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          legend.position=c(0.9,0.9))+
        guides(colour = guide_legend(override.aes = list(alpha = 1)),alpha=FALSE)+
        labs(x="Principal component 1", y= "Principal component 2")+
        ggtitle("Population") -> plot.synthetic_imp.true
      
      qplot(PC.1,PC.2, col = factor(HAND), alpha = I(0.1))+
        scale_colour_discrete(name= "Frascati",labels=c("Normal","Impaired")) +
        theme_classic() +
        theme(
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          legend.position=c(0.9,0.9))+
        guides(colour = guide_legend(override.aes = list(alpha = 1)),alpha=FALSE)+
        labs(x="Principal component 1", y= "Principal component 2")+
        ggtitle("Frascati") -> plot.synthetic_imp.HAND
      
      qplot(PC.1,PC.2, col = factor(GDS), alpha = I(0.1))+
        scale_colour_discrete(name= "GDS",labels=c("Normal","Impaired")) +
        theme_classic() +
        theme(
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          legend.position=c(0.9,0.9))+
        guides(colour = guide_legend(override.aes = list(alpha = 1)),alpha=FALSE)+
        labs(x="Principal component 1", y= "Principal component 2") +
        ggtitle("GDS")-> plot.synthetic_imp.GDS
      
      qplot(PC.1,PC.2, col = factor(mahal_dist<critical), alpha = I(0.1))+
        scale_colour_discrete(name= "MD",labels=c("Normal","Impaired")) +
        theme_classic() +
        theme(
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          legend.position=c(0.9,0.9))+
        guides(colour = guide_legend(override.aes = list(alpha = 1)),alpha=FALSE)+
        labs(x="Principal component 1", y= "Principal component 2")+
        ggtitle("Mahalanobis Distance") -> plot.synthetic_imp.mahal
      
      print(plot_grid(plot.synthetic_imp.true, plot.synthetic_imp.HAND, plot.synthetic_imp.GDS, plot.synthetic_imp.mahal, align = "hv"))
    })
    
    
  })
  
  output$table <- renderTable({ 
    
    # get the dynamic values
    synthetic_pop2 <- synthetic_pop()$synthetic_pop2
    HAND <- synthetic_pop()$HAND
    GDS <- synthetic_pop()$GDS
    mahal_dist <- synthetic_pop()$mahal_dist
    critical <- synthetic_pop()$critical 
    
    #crunch the numbers
    sum(HAND)/syn_n ->HAND.prev
    sum(GDS)/syn_n -> GDS.prev
    sum(mahal_dist<critical)/syn_n -> mahal.prev
    sum(HAND == TRUE & synthetic_pop2$label==2)/sum(HAND == TRUE) -> HAND.ppv
    sum(HAND == FALSE & synthetic_pop2$label == 1)/sum(HAND == FALSE) -> HAND.npv
    sum(GDS == TRUE & synthetic_pop2$label==2)/sum(GDS == TRUE) -> GDS.ppv
    sum(GDS == FALSE & synthetic_pop2$label == 1)/sum(GDS == FALSE) -> GDS.npv
    sum(mahal_dist<critical & synthetic_pop2$label==2)/sum(mahal_dist<critical) -> mahal.ppv
    sum(mahal_dist>critical & synthetic_pop2$label==1)/sum(mahal_dist>critical) -> mahal.npv
    (sum(HAND== TRUE & synthetic_pop2$label==2) + sum(HAND==FALSE & synthetic_pop2$label==1))/syn_n -> HAND.acc   
    (sum(GDS==TRUE & synthetic_pop2$label==2) + sum(GDS==FALSE & synthetic_pop2$label==1))/syn_n -> GDS.acc
    (sum(mahal_dist<critical & synthetic_pop2$label==2) + sum(mahal_dist>critical & synthetic_pop2$label==1))/syn_n -> mahal.acc
    
    table <- rbind("Frascati", "GDS", "Mahalanobis distance")
    table <- cbind(table, rbind(paste(round(HAND.prev, digits = 2)*100, "%",sep = ""), paste(round(GDS.prev, digits = 2)*100, "%",sep = ""), paste(round(mahal.prev, digits = 2)*100, "%",sep = "")))
    table <- cbind(table, rbind(paste(round(HAND.ppv, digits = 2)*100, "%",sep = ""), paste(round(GDS.ppv, digits = 2)*100, "%",sep = ""), paste(round(mahal.ppv, digits = 2)*100, "%",sep = "")))
    table <- cbind(table, rbind(paste(round(HAND.npv, digits = 2)*100, "%",sep = ""), paste(round(GDS.npv, digits = 2)*100, "%",sep = ""), paste(round(mahal.npv, digits = 2)*100, "%",sep = "")))
    table <- cbind(table, rbind(paste(round(HAND.acc, digits = 2)*100, "%",sep = ""), paste(round(GDS.acc, digits = 2)*100, "%",sep = ""), paste(round(mahal.acc, digits = 2)*100, "%",sep = "")))
    colnames(table) <- c("Diagnostic method", "Prevalence of 'impairment'", "PPV", "NPV", "Accuracy")
    
    #return the table
    table
  }, striped = TRUE, bordered = TRUE, align = "lrrrr")
  
})

shinyApp(ui = ui, server = server)