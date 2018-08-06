#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws charts resulting from SD ERCOT cap growth model
ui <- fluidPage(
   
   # Application title
   titlePanel("System Dynamics Model for Capacity Expansion in ERCOT"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         sliderInput("carbonTax",
                     "Price on Carbon, in $ per ton:",
                     min = 0,
                     max = 100,
                     value = 0),
         
         sliderInput("windImprovRate",
                     "Wind Technology Improvement Rate, in percent per year:",
                     min = 0.01,
                     max = 0.1,
                     value = 0.05),
         
         sliderInput("marginalCostCCNG",
                     "Fuel cost for CCNG, in $ per MWh:",
                     min = 40,
                     max = 100,
                     value = 50)
      ),
      
      # Show a plot of the generated distribution
      mainPanel('Model Outputs',
         fluidRow(
           plotOutput('gencap'),
           plotOutput('genMix'),
           plotOutput('profitability'),
           plotOutput('emissions'),
           plotOutput('price')
         )
      )
   )
)

# Define server logic 
# Contains the model, and outputs
server <- function(input, output) {
  
  {myColors = c(Wind='#006600',Solar='#66FF00',FuelOil='#999999',NaturalGas='#FF9900',Hydro='#0000FF',
                Coal = '#222222',Coal.Market='#666666',Coal.Self='#444444',
                Nuclear='#FF0000',Waste.Disposal.Services='#00CCFF',
                Waste.Heat='#996600',Other='#993333',CCNG='#FF9900',CTNG='#999999',
                'Other Fuel' = '#993333',Biomass = '#993333','Natural Gas' = '#FF9900',
                'Petroleum Products' = '#999999',Uranium = '#FF0000', Water = '#0000FF')} # Plot colors
  
  
  require(deSolve,quietly = T)
  require(data.table,quietly = T)
  require(ggplot2,quietly = T)
   
  START <- 2017; FINISH <- 2037; STEP<- 1/32
  simtime <- seq(START, FINISH, by=STEP)
  
  # Define Stocks 
  {
    aInitPeakDemand <- 71250
    
    genTypes <- c('Wind','Nuclear','Coal','CCNG','CTNG')
    
    sCap <- c(Wind    = 19000,
              Nuclear = 0.17 * aInitPeakDemand,
              Coal    = 0.11 * aInitPeakDemand,
              CCNG    = 0.55 * aInitPeakDemand,
              CTNG    = 0.15 * aInitPeakDemand)
    
    aConstTime <- c(Wind    = 2,
                    Nuclear = 10,
                    Coal    = 7,
                    CCNG    = 3,
                    CTNG    = 1.5)
    
    aAssetLife <- c(Wind    = 25,
                    Nuclear = 35,
                    Coal    = 30,
                    CCNG    = 25,
                    CTNG    = 25)
    
    sCapUnder = sCap * aConstTime / aAssetLife
    
  }
  
  stocks <- c(SCap = sCap,SCapUnder = sCapUnder)
  
  # Model functions
  
  calcNetLoadParam <- function(windPct){
    P0 <- -0.0604 * windPct + 0.9968
    P20 <- -0.1921* windPct + 0.6467
    P75 <- -0.4511 * windPct + 0.4876
    P100 <- -0.7362*windPct + 0.3772
    
    return(c(P0=P0,P20=P20,P75=P75,P100=P100))
  }
  
  calcPTC <- function(tm){
    ptc <- pmin(30,pmax(0,30 * (START + 10 - tm)/10 ))
    return(ptc)
  }
  
  calcPriceByFuel <- function(sCap,aMarginalCost,tm,vPeakDemand,aEmissions,aPriceOnCarbon,aScarcityTable) {
    desiredCTNG <- vPeakDemand - sum(sCap[names(sCap) %in% c('Nuclear','Coal','CCNG')])
    
    ctngPriceMultiplier <- aScarcityTable(sCap['CTNG']/desiredCTNG)
    
    priceByFuel <- aMarginalCost 
    
    priceByFuel['CTNG'] <- ctngPriceMultiplier * priceByFuel['CTNG']
    
    priceByFuel['Wind'] <- priceByFuel['Wind'] - calcPTC(tm)
    
    priceByFuel <- priceByFuel + aEmissions * aPriceOnCarbon
    
    return(priceByFuel)
    
  }
  
  calcEnergyByFuel <- function (sCap,vPeakDemand, aQuantiles,aLoadDurationCurve) {
    windFraction <- as.numeric(sCap['Wind']/vPeakDemand)
    
    capFraction <- sCap[names(sCap)!='Wind']/vPeakDemand
    
    netLoadParam <- calcNetLoadParam(windFraction)
    
    hourPct <- approx(x = netLoadParam,
                      y = aQuantiles,
                      c(Wind = 0,cumsum(capFraction)),
                      rule = 2)
    hourPct <- setNames(hourPct$y,names(hourPct$x))
    
    clearingHourFct <- data.table::shift(hourPct,fill = 1) - hourPct 
    
    energyGenByFuel <- sapply(seq(1,length(genTypes)),function(k){
      fuel<-genTypes[k]
      if ( k > 1) {
        cheaperFuel <- genTypes[k-1]
        fromPct <- hourPct[k]
        toPct <- min(hourPct[k-1],
                     approx(x = netLoadParam,
                            y = aQuantiles,
                            xout = 0,rule = 2)$y)
        en <- MESS::auc(x = aQuantiles,y = netLoadParam,
                        from = fromPct, to = toPct) - ifelse(k>2,cumsum(capFraction)[cheaperFuel]*(toPct - fromPct),0)
        en <- en + hourPct[fuel] * capFraction[fuel]
      } else {
        totalDemand <- MESS::auc(x = aQuantiles, y = aLoadDurationCurve,
                                 from = 0 ,to = 1)
        h0 <- approx(x = netLoadParam,y=aQuantiles,0,rule = 2)$y
        totalNetDemand <- MESS::auc(x = aQuantiles, y = netLoadParam,
                                    from = 0 ,to = h0)
        en <- c(Wind = totalDemand - totalNetDemand)
      }
      return(en)
    })
    
    return(energyGenByFuel)
    
  }
  
  calcRevPerMwhByFuel <- function(sCap,aMarginalCost,tm,vPeakDemand,aEmissions,aPriceOnCarbon, aQuantiles,aLoadDurationCurve ,aScarcityTable) {
    
    windFraction <- as.numeric(sCap['Wind']/vPeakDemand)
    
    capFraction <- sCap[names(sCap)!='Wind']/vPeakDemand
    
    netLoadParam <- calcNetLoadParam(windFraction)
    
    hourPct <- approx(x = netLoadParam,
                      y = aQuantiles,
                      c(Wind = 0,cumsum(capFraction)),
                      rule = 2)
    hourPct <- setNames(hourPct$y,names(hourPct$x))
    
    clearingHourFct <- data.table::shift(hourPct,fill = 1) - hourPct 
    
    prices <- calcPriceByFuel(sCap,aMarginalCost,tm,vPeakDemand,aEmissions,aPriceOnCarbon,aScarcityTable)
    
    energyGenByFuel <- calcEnergyByFuel(sCap,vPeakDemand, aQuantiles,aLoadDurationCurve)
    
    revenuePerMwh <- sapply(seq(1,length(genTypes)),function(k){
      
      fuel<-genTypes[k]
      
      if (k == length(genTypes)){ # The most expensive fuel sets the price whenever it operates
        rvn <- prices[fuel]
      } else if (k > 1) { # For fuels different than wind
        
        
        # Revenue from hours cleared by more expensive fuels
        revHoursMoreExp <- sapply(seq(length(genTypes),k+1), function(p){
          expFuel <- genTypes[p]
          prices[expFuel]*clearingHourFct[expFuel]*capFraction[fuel] # This value is in 'dollars', ie, $/MWh * h * MW
        })
        
        revHoursMoreExp <- as.numeric(sum(revHoursMoreExp)) # This value is in $. As numeric to un-name it
        
        # Revenue from hours when the gen type is clearing
        expFuel <- genTypes[k+1]
        enGenWhenClearing <- energyGenByFuel[fuel] - hourPct[fuel] * capFraction[fuel]
        revHoursClearing <- enGenWhenClearing * prices[fuel] # This value is in dollars
        
        # Sum of revenues divided by total energy will give $/MWh for that gen type
        
        totalRev <- revHoursMoreExp + revHoursClearing
        
        rvn <- totalRev/energyGenByFuel[fuel]
        
      } else if (k == 1) { # This is the loop to calculate the revenue for wind
        
        # Revenue for wind: calculate how much wind energy is generated for each clearing fuel
        # Total Revenue divided by total wind gen equals revenue per MWh of wind
        # Wind revenue when it's clearing is equal to it's marginal cost, which could be negative
        
        h0 <- approx(x = netLoadParam,y=aQuantiles,0,rule = 2)$y # Indicates hourPct when wind sets the price
        
        # Revenue from hours cleared by more expensive fuels
        revHoursMoreExp <- sapply(seq(length(genTypes),2), function(p){
          expFuel <- genTypes[p]
          
          
          fromH <- hourPct[p]
          toH <- ifelse(p > 2 , hourPct[p-1], h0)
          
          totalDemand <- MESS::auc(x = aQuantiles, y = aLoadDurationCurve,
                                   from = fromH ,to = toH)
          
          totalNetDemand <- MESS::auc(x = aQuantiles, y = netLoadParam,
                                      from = fromH ,to = toH)
          
          # Difference between the two is generated by wind
          windGen <- totalDemand - totalNetDemand
          
          prices[expFuel] * windGen
          
        })
        
        revHoursMoreExp <- as.numeric(sum(revHoursMoreExp))
        
        # Revenue when wind is clearing the price
        # Calculate energy when not clearing and subtract from total energy
        totalDemand <- MESS::auc(x = aQuantiles, y = aLoadDurationCurve,
                                 from = 0 ,to = h0)
        
        totalNetDemand <- MESS::auc(x = aQuantiles, y = netLoadParam,
                                    from = 0 ,to = h0)
        
        # Difference between the two is generated by wind
        windGenNotClearing <- totalDemand - totalNetDemand
        
        windGenClearing <- energyGenByFuel['Wind'] - windGenNotClearing
        
        revHoursClearing <- windGenClearing * prices['Wind'] # This value is in dollars
        
        # Sum of revenues divided by total energy will give $/MWh for that gen type
        
        totalRev <- revHoursMoreExp + revHoursClearing
        
        rvn <- totalRev/energyGenByFuel['Wind']
      } else {stop('You shouldnt be here')}
      
      return(rvn)
      
    })
    
    return(revenuePerMwh)
  }
  
  calcCfByFuel <- function(sCap,vPeakDemand, aQuantiles,aLoadDurationCurve) {
    
    energyGenByFuel <- calcEnergyByFuel(sCap,vPeakDemand, aQuantiles,aLoadDurationCurve)
    
    capFraction <- sCap/vPeakDemand
    
    cf <- energyGenByFuel / capFraction
    
    return(cf)
  }
  
  calcSpecificConstCost <- function(aInitSpecificConstCost,aTechImprovRate,tm,START = 2017) {
    constCost <- aInitSpecificConstCost * (1-aTechImprovRate)^(tm - START)
    return(constCost)
  }
  
  calcCrf <- function (aInterestRate,aInvestHorizon) {
    crf <- ( aInterestRate * ( (1+aInterestRate) ^ (aInvestHorizon) ) )/(  ((1+aInterestRate)^(aInvestHorizon)) -1 )
    return(crf)
    
  }
  
  calcLcoe <- function(cf,aMarginalCost,constCost,aOperatingCost,crf,aEmissions,aPriceOnCarbon,tm) {
    
    vMarginalCost <- aMarginalCost + aEmissions * aPriceOnCarbon
    
    vMarginalCost['Wind'] <- vMarginalCost['Wind'] - calcPTC(tm)
    
    lcoe <- ((constCost * crf + aOperatingCost)/(pmax(cf,0.005)*8760)) + aMarginalCost
    return(lcoe)
  }
  
  model <- function(time,stocks,auxs) {
    if (!'windImprovRate' %in% names(auxs)) stop('Variable not flowing in')
    
    with(as.list(stocks,auxs),{
      
      {
        aConstTime <- c(Wind    = 2,
                        Nuclear = 10,
                        Coal    = 7,
                        CCNG    = 3,
                        CTNG    = 1.5)
        
        aAssetLife <- c(Wind    = 25,
                        Nuclear = 35,
                        Coal    = 30,
                        CCNG    = 25,
                        CTNG    = 25)
        
        aInitPeakDemand <- 71250
        
        genTypes <- c('Wind','Nuclear','Coal','CCNG','CTNG')
        
        aLoadDurationCurve <- c(D0 = 1, D20 = .6467, D75 = .4876, D100 = 0.3772)
        
        tDesiredCapGrowth <- approxfun(x = c(-.5,-.4,-.32,-.25,-.2,-.2,0,.1,.2,.25,.32,.4,.5),
                                       y = c(-.35,-.28,-.18,-.1,-.06,-.03,0,.05,.16,.22,.28,.32,.35),
                                       rule = 2,ties = 'ordered')
        
        aQuantiles <- c(0,.2,.75,1)
        
        ademandGrothRate <- 0.01
        
        aTimeToAdjustSupplyLine <- 2
        
        aEmissions <- c(Wind    = 0.02,
                        Nuclear = 0.012,
                        Coal    = 0.992,
                        CCNG    = 0.394,
                        CTNG    = 0.563)
        
        aInterestRate = 0.1
        aInvestHorizon = 20
        
        aOperatingCost <- c(Wind = 40000,
                            Nuclear = 90000,
                            Coal = 30000,
                            CCNG = 15370,
                            CTNG = 7040)
        
        aInitSpecificConstCost <- c(Wind = 1.8e+006,
                                    Nuclear = 4.8e+006,
                                    Coal = 3.5e+006,
                                    CCNG = 917000,
                                    CTNG = 625000)
        
        aTechImprovRate <- c(Wind = windImprovRate,
                             Nuclear = 0,
                             Coal = 0.005,
                             CCNG = 0.02,
                             CTNG = 0.02)
        
        aMarginalCost <- c(Wind = 0,
                           Nuclear = 11,
                           Coal = 31,
                           CCNG = marginalCostCCNG,
                           CTNG = 7 * marginalCostCCNG / 5)
        
        aPriceOnCarbon  <- carbonTax
        
        aScarcityTable <- approxfun(x = c(1,1.3),
                                    y = c(15,1),
                                    rule = 2)
      } # All static data inputs
      
      SCapUnder <- stocks[grepl('SCapUnder',names(stocks))]
      
      SCap <- stocks[!grepl('SCapUnder',names(stocks))]
      
      SCap <- setNames(SCap,nm = genTypes)
      SCapUnder <- setNames(SCapUnder,nm = genTypes)
      
      
      # Variables
      vPeakDemand <- aInitPeakDemand * (ademandGrothRate + 1) ^(time - START)
      
      # Profitability Module
      vRevByFuel <- calcRevPerMwhByFuel(SCap,aMarginalCost,time,vPeakDemand,aEmissions,aPriceOnCarbon, aQuantiles,aLoadDurationCurve ,aScarcityTable)
      
      vCostByFuel <- calcLcoe(cf = calcCfByFuel(sCap,vPeakDemand, aQuantiles,aLoadDurationCurve),
                              aMarginalCost,
                              constCost = calcSpecificConstCost(aInitSpecificConstCost,aTechImprovRate,tm = time,START),
                              aOperatingCost,
                              crf = calcCrf(aInterestRate,aInvestHorizon),tm = time,aEmissions,aPriceOnCarbon)
      
      vProf <- (vRevByFuel / vCostByFuel) - 1
      
      vDesiredGrowth <- tDesiredCapGrowth(vProf)
      
      
      # Flows = construction start rate, construction finish rate, retirement
      fConstStart <- pmax(0,(SCap / aAssetLife) + (SCap * vDesiredGrowth) + ((aConstTime/aAssetLife)*SCap - SCapUnder))
      
      fConstCompletion <- SCapUnder / aConstTime
      fRetirement <- SCap / aAssetLife
      
      # Net flows into the stocks
      dsCap_dt      <- fConstCompletion - fRetirement
      dsCapUnder_dt <- fConstStart      - fConstCompletion
      
      # Convert flows into matrix
      flows <- c(dsCap_dt,dsCapUnder_dt)
      
      # Desired outputs
      oGenMix        <- calcEnergyByFuel(SCap,vPeakDemand, aQuantiles,aLoadDurationCurve)  / (sum(calcEnergyByFuel(SCap,vPeakDemand, aQuantiles,aLoadDurationCurve)))
      oEmissions     <- sum(oGenMix * aEmissions) / (sum(calcEnergyByFuel(SCap,vPeakDemand, aQuantiles,aLoadDurationCurve)))
      oProfitability <- vProf
      oPrice <- weighted.mean(vRevByFuel,oGenMix)
      
      return(list(flows,
                  oGenMix = oGenMix,
                  oEmissions = oEmissions,
                  oProfitability = oProfitability,
                  oPrice = oPrice))
      
    })
  }
  
  results <- reactive({
    # input <- list(carbonTax = 0, windImprovRate = 0.05 , marginalCostCCNG = 50)
    
    carbonTax <- input$carbonTax
    windImprovRate <- input$windImprovRate
    marginalCostCCNG <- input$marginalCostCCNG
    
    auxs <- list(carbonTax = carbonTax, windImprovRate = windImprovRate , marginalCostCCNG = marginalCostCCNG)
    
    results <- data.table(deSolve::ode(y = stocks,times = simtime,func = model,parms = auxs, method = 'euler'))
    setnames(results,'time','Time')
    
    return(results)
    
  })
  
  
  
  
  output$genCap <- renderPlot({
    genCap <- results()[, colnames(results()) == 'Time' | grepl(x = colnames(results()),pattern = 'SCap') & !grepl(x = colnames(results()),pattern = 'SCapUnder'),with = F]
    setnames(genCap,colnames(genCap),gsub('SCap.','',colnames(genCap)))
    
    genCap <- melt.data.table(genCap,id.vars = 'Time',variable.name = 'FuelType',value.name = 'Capacity')
    
    ggplot(genCap,aes(Time,Capacity,color = FuelType))+
      geom_line(size = 0.8)+
      theme_minimal()+
      scale_color_manual(values = myColors)
    })
  
  
  
  output$genMix <- renderPlot({
    genMix <- results()[, colnames(results()) == 'Time' | grepl(x = colnames(results()),pattern = 'oGenMix') ,with = F]
    setnames(genMix,colnames(genMix),gsub('oGenMix.','',colnames(genMix)))
    
    genMix <- melt.data.table(genMix,id.vars = 'Time',variable.name = 'FuelType',value.name = 'Generation')
    
    ggplot(genMix,aes(Time,Generation,color = FuelType))+
      geom_line(size = 0.8)+
      theme_minimal()+
      scale_color_manual(values = myColors)
    })
  
  output$profitability <- renderPlot({
    profitability <- results()[, colnames(results()) == 'Time' | grepl(x = colnames(results()),pattern = 'oProfitability.') ,with = F]
    setnames(profitability,colnames(profitability),gsub('oProfitability.','',colnames(profitability)))
    
    profitability <- melt.data.table(profitability,id.vars = 'Time',variable.name = 'FuelType',value.name = 'Profitability')
    
    ggplot(profitability,aes(Time,Profitability,color = FuelType))+
      geom_line(size = 0.8)+
      theme_minimal()+
      scale_color_manual(values = myColors)
    })
  
  output$price <- renderPlot({ggplot(results(),aes(Time,oPrice))+geom_line()+
    theme_minimal()})
  
  output$emissions <- renderPlot({ggplot(results(),aes(Time,oEmissions))+geom_line()+
    theme_minimal()})
  
}

# Run the application 
shinyApp(ui = ui, server = server)

