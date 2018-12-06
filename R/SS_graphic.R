#' Graphical User Interface for Sample Size Calculation
#'
#' @import shiny
#' @import ggplot2
#'
#' @return nothing. But while processing, a graphic shell appears, which does all the calculations implemented here
#' @export
#'
#' @examples
#'\dontrun{SampleSize_graphic()}
#'
SampleSize_graphic <- function(){
  # Define User Interface
  ui <- shiny::fluidPage(
    # title
    shiny::titlePanel("Sample Size Calculation for Exponentially Distributed, Right Censored Two Sample Tests"),
    # Sidebar
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        # Static Parameters
        shiny::sliderInput(inputId = "invparam",
                           label = "Median Survival in Reference Group:",
                           min = 0,
                           max = 10,
                           value = 2,
                           ticks = FALSE,
                           step = 1/4),
        shiny::sliderInput(inputId = "effect",
                           label = "Ratio of Median Survival Times (= Experimental Group / Reference Group):",
                           min = 0,
                           max = 5,
                           value = 1.4,
                           ticks = FALSE,
                           step = .1),
        shiny::sliderInput(inputId = "alpha",
                           label = "Error of Type 1:",
                           min = 0,
                           max = .2,
                           value = .05,
                           ticks = FALSE,
                           step = .005),

        shiny::sliderInput(inputId = "beta",
                           label = "Error of Type 2:",
                           min = 0,
                           max = .4,
                           value = .2,
                           ticks = FALSE,
                           step = .01),
        shiny::sliderInput(inputId = "groupratio",
                           label = "Randomization Probability for Reference Group:",
                           min = 0,
                           max = 1,
                           value = .5,
                           ticks = FALSE,
                           step = .01),

        shiny::selectInput(inputId = "censoring",
                           label = "Censoring Scheme:",
                           choices = c("Constant Follow-Up",
                                       "Constant Follow-Up with Random Dropout",
                                       "Follow-Up with Accrual",
                                       "Follow-Up with Accrual and Random Dropout"),
                           selected="Follow-Up with Accrual and Random Dropout"),
        # Optional parameters
        shiny::conditionalPanel(
          condition = "input.censoring == 'Constant Follow-Up'",
          shiny::sliderInput(inputId = "time1",
                             label = "Study Time:",
                             min = 0,
                             max = 10,
                             value = 5,
                             ticks = TRUE,
                             step = 1/4)
        ),
        shiny::conditionalPanel(
          condition = "input.censoring == 'Constant Follow-Up with Random Dropout'",
          shiny::sliderInput(inputId = "time2",
                             label = "Study Time:",
                             min = 0,
                             max = 10,
                             value = 5,
                             ticks = TRUE,
                             step = 1/4),
          shiny::sliderInput(inputId = "dropout2",
                             label = "Dropout Ratio per Time Unit:",
                             min = 0,
                             max = .1,
                             value = .05,
                             ticks = TRUE,
                             step = .005)
        ),
        shiny::conditionalPanel(
          condition = "input.censoring == 'Follow-Up with Accrual'",
          shiny::sliderInput(inputId = "time3",
                             label = "Study Time:",
                             min = 0,
                             max = 10,
                             value = 5,
                             ticks = TRUE,
                             step = 1/4),
          shiny::sliderInput(inputId = "accrual3",
                             label = "Accrual Time:",
                             min = 0,
                             max = 10,
                             value = 2,
                             ticks = TRUE,
                             step = 1/4)
        ),
        shiny::conditionalPanel(
          condition = "input.censoring == 'Follow-Up with Accrual and Random Dropout'",
          shiny::sliderInput(inputId = "time4",
                             label = "Study Time:",
                             min = 0,
                             max = 10,
                             value = 5,
                             ticks = TRUE,
                             step = 1/4),
          shiny::sliderInput(inputId = "accrual4",
                             label = "Accrual Time:",
                             min = 0,
                             max = 10,
                             value = 2,
                             ticks = TRUE,
                             step = 1/4),
          shiny::sliderInput(inputId = "dropout4",
                             label = "Dropout Ratio per Time Unit:",
                             min = 0,
                             max = .1,
                             value = .05,
                             ticks = TRUE,
                             step = .005)
        )
      ),
      # Main panel for displaying outputs
      shiny::mainPanel(
        shiny::h3("Estimated Sample Sizes"),
        shiny::tableOutput(outputId = "tbl"),
        shiny::textOutput(outputId  = "schoenfeld_n"),
        shiny::h3("Densities of Event and Censoring Times"),
        shiny::plotOutput(outputId = "censPlot")
      )
    )
  )
  # Server Design
  server <- function(input, output) {
    # Plot Densities for better understanding of "what is simulated here"
    output$censPlot <- shiny::renderPlot({
      # create plot data
      futime  <- switch (input$censoring,
                         "Constant Follow-Up" = input$time1,
                         "Constant Follow-Up with Random Dropout" = input$time2,
                         "Follow-Up with Accrual" = input$time3,
                         "Follow-Up with Accrual and Random Dropout" = input$time4)
      actime  <- switch (input$censoring,
                         "Constant Follow-Up" = .1,
                         "Constant Follow-Up with Random Dropout" = .1,
                         "Follow-Up with Accrual" = input$accrual3,
                         "Follow-Up with Accrual and Random Dropout" = input$accrual4)
      dropout <- switch (input$censoring,
                         "Constant Follow-Up" = 0,
                         "Constant Follow-Up with Random Dropout" = input$dropout2,
                         "Follow-Up with Accrual" = 0,
                         "Follow-Up with Accrual and Random Dropout" = input$dropout4)
      dout_frac <- dropout * futime
      cens_frac <- 1 - dout_frac

      dt <- data.frame(ref = seq(from = 0, to = 10, by = .02))
      dt$x <- dexp(x = dt$ref, rate = 1 * log(2) / input$invparam)
      dt$y <- dexp(x = dt$ref, rate = 1 * log(2) / input$invparam / input$effect)
      dt$c <- dropout * (dt$ref <= futime) + cens_frac / actime * ((dt$ref >= (futime - actime)) & dt$ref <= futime)

      # plot object
      ggplot2::ggplot() +
        ggplot2::geom_area(data = dt,
                           mapping = ggplot2::aes(x=ref,
                                                  y=x,
                                                  fill = "Event in Reference Group"),
                           alpha = .5) +
        ggplot2::geom_area(data = dt,
                           mapping = ggplot2::aes(x=ref,
                                                  y=y,
                                                  fill = "Event in Experimental Group"),
                           alpha = .5) +
        ggplot2::geom_area(data = dt,
                           mapping = ggplot2::aes(x=ref,
                                                  y=c,
                                                  fill = "Censoring"),
                           alpha = .5) +
        ggplot2::xlim(0, 10) +
        ggplot2::labs(x = "Passed Time",
                      y = "Density",
                      title = "") +
        ggplot2::scale_fill_discrete(name = "Densities") +
        ggplot2::theme(legend.position="bottom",
                       axis.text=element_text(size=14),
                       axis.title=element_text(size=14),
                       legend.text=element_text(size=14),
                       legend.title=element_text(size=14, face="bold"))
    })
    # Calculate Sample Sizes
    output$tbl <- shiny::renderTable({
      # Define Accuracy
      nSamples <- 1e6
      # Fit randomization probabilities for the groups into this packages framework
      gf <- (1-input$groupratio)/ input$groupratio
      # Simulate Random numbers to gain censoring probabilities
      dist_x  <- "exp"
      param_x <- 1 * log(2) / input$invparam
      dist_y  <- "exp"
      param_y <- 1 * log(2) / input$invparam / input$effect

      dist_c  <- switch (input$censoring,
                         "Constant Follow-Up" = "followup",
                         "Constant Follow-Up with Random Dropout" = "followup",
                         "Follow-Up with Accrual" = "accrual-followup",
                         "Follow-Up with Accrual and Random Dropout" = "accrual-followup")
      param_c <- switch (input$censoring,
                         "Constant Follow-Up" = c(0, input$time1),
                         "Constant Follow-Up with Random Dropout" = c(input$dropout2, input$time2),
                         "Follow-Up with Accrual" = c(0, input$time3, input$accrual3),
                         "Follow-Up with Accrual and Random Dropout" = c(input$dropout4, input$time4, input$accrual4))

      dx <- observation_probability(dist_E = dist_x, param_E = param_x,
                                    dist_C = dist_c, param_C = param_c,
                                    N = nSamples)
      dy <- observation_probability(dist_E = dist_y, param_E = param_y,
                                    dist_C = dist_c, param_C = param_c,
                                    N = nSamples)
      #
      # Call Sample Size Formulas implemented here
      SS_Wald     <- ExpSurvCI::SampleSize_wald(k1 = 1/input$effect, dx = dx, dy = dy, alpha = input$alpha, beta = input$beta, c = gf)
      SS_Gradient <- ExpSurvCI::SampleSize_gradient(k1 = 1/input$effect, dx = dx, dy = dy, alpha = input$alpha, beta = input$beta, c = gf)
      SS_Schoenfeld <- ExpSurvCI::SampleSize_schoenfeld(k1 = 1/input$effect, dx = dx, dy = dy, alpha = input$alpha, beta = input$beta, c = gf)
      # return Results
      data.frame(Number = c("Overall", "Reference Group", "Experimental Group"),
                 Wald = as.integer(unlist(SS_Wald)),
                 Gradient = as.integer(unlist(SS_Gradient)),
                 Schoenfeld = as.integer(unlist(SS_Schoenfeld)),
                 "P.Event" = c(dx / (1+gf) + gf * dy / (1+gf), dx, dy))
    })

    output$schoenfeld_n <- shiny::renderText({
      z1  <- qnorm(1 - input$alpha / 2)
      z2  <- qnorm(1 - input$beta)
      p_r <- input$groupratio * (1 - input$groupratio)
      l_e <- log(input$effect)^2
      sprintf("According to Schoenfeld formula, %.0f events are required.",
              ceiling((z1 + z2)^2 / (p_r * l_e)))
    })
  }
  # final call
  shiny::shinyApp(ui, server)
}
