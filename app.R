library(shiny)
library(bslib)
library(actuar)
library(readxl)
library(ggplot2)
library(extraDistr)
library(DT)
library(jsonlite)

options(scipen = 999)

# ==============================================================================
# Norwegian Oil Extraction Dataset (n = 128)
# ==============================================================================
norwegian_oil_data <- c(
  5409000, 9560000, 16228000, 8177000, 30065000, 27990000, 27990000, 36984000,
  47299000, 15662000, 15662000, 20190000, 34908000, 4843000, 7674000, 1069000,
  755000, 58684000, 41072000, 39752000, 314000, 9749000, 9749000, 2076000,
  72333000, 41890000, 172089000, 23775000, 9623000, 9623000, 189000, 9120000,
  10693000, 190204000, 52834000, 5661000, 5472000, 90259000, 190895000, 174165000,
  2264000, 36166000, 63338000, 63338000, 126000, 12580000, 7862000, 12768000,
  2076000, 75163000, 99064000, 60822000, 32896000, 81453000, 25851000, 314000,
  15850000, 15976000, 188694000, 29877000, 184480000, 6227000, 38494000, 559792000,
  21826000, 12831000, 2893000, 65540000, 22769000, 189000, 2579000, 1195000,
  128249000, 31763000, 7296000, 46796000, 2705000, 2264000, 2264000, 49249000,
  3900000, 54533000, 41000000, 9498000, 1887000, 22958000, 4780000, 3900000,
  45916000, 189000, 8177000, 481862000, 69125000, 23398000, 27864000, 19247000,
  22580000, 16542000, 2327000, 47991000, 47991000, 1447000, 45287000, 39940000,
  27298000, 12517000, 18178000, 503000, 189000, 12014000, 28304000, 6101000,
  6353000, 755000, 1006000, 200204000, 15410000, 29562000, 29562000, 47111000,
  6667000, 19184000, 692000, 4214000, 4088000, 410409000, 451545000, 39500000
)

# ==============================================================================
# Document 1 Theme & Custom Styling
# ==============================================================================
app_theme <- bs_theme(
  version = 5,
  bootswatch = "flatly",
  primary = "#002855",    # UP Navy
  secondary = "#41748d",  # Slate Blue
  success = "#2a9d8f",    # Muted Teal
  warning = "#e9c46a",    # Amber
  danger = "#e63946",     # Crimson
  base_font = "Inter, sans-serif"
)

# Custom Helper Functions for App 1 (Pareto Simulation)
rpareto_custom <- function(v, a, b) {
  b * ((1 - runif(v))^(-1 / a))
}

scale_inflated_rpareto <- function(v, a, b, lambda) {
  lambda * b * ((1 - runif(v))^(-1 / a))
}

# ==============================================================================
# Helper Functions for App 2 (MCEM & Gibbs Sampler)
# ==============================================================================
mcmc_app2 <- function(x, hyper_parms, pareto_parms, R = 100, S = 100) {
  n      <- length(x)
  a1     <- hyper_parms[1]
  a2     <- hyper_parms[2]
  b_star <- hyper_parms[3]
  lambda <- hyper_parms[4]
  b1     <- hyper_parms[5]
  b2     <- hyper_parms[6]
  t1     <- hyper_parms[7]
  t2     <- hyper_parms[8] 
  
  alpha   <- pareto_parms[1] 
  beta    <- pareto_parms[2]
  epsilon <- pareto_parms[3]  
  theta   <- pareto_parms[4]
  
  alphas   <- numeric(R)
  betas    <- numeric(R)
  deltas   <- matrix(NA_real_, nrow = n, ncol = R)
  epsilons <- numeric(R)
  thetas   <- numeric(R)
  
  delta <- x > quantile(x, probs = 1 - epsilon)
  k     <- sum(delta)
  
  for (r in seq_len(R)) {
    for (s in seq_len(S)) {
      # 1. Update Alpha
      alpha_shape <- a1 + n
      alpha_rate  <- a2 + sum(log(x)) - n * log(theta) - k * log(beta)
      alpha       <- qgamma(runif(1) * (pgamma(1, shape = alpha_shape, rate = alpha_rate, lower.tail = FALSE)), 
                            shape = alpha_shape, rate = alpha_rate, lower.tail = FALSE)
      
      # 2. Update Epsilon
      epsilon     <- rbeta(1, shape1 = b1 + k, shape2 = b2 + n - k)
      
      # 3. Update Delta
      p_outlier   <- (x >= beta * theta) * (beta^alpha * epsilon) / (1 - epsilon + beta^alpha * epsilon)
      delta       <- rbinom(n, size = 1, prob = p_outlier)
      k           <- sum(delta)
      
      # 4. Update Beta
      beta_upper <- if (k > 0) { 
        max(c(b_star, min(x[delta == 1]) / theta)) 
      } else { 
        max(c(1.1, b_star + min(x) / theta)) 
      }
      beta <- qgamma(
        pgamma(b_star, shape = alpha * k + 1, rate = lambda) + 
          runif(1) * (pgamma(beta_upper, shape = alpha * k + 1, rate = lambda) - 
                        pgamma(b_star, shape = alpha * k + 1, rate = lambda)),
        shape = alpha * k + 1,
        rate = lambda
      )
      beta <- ifelse(is.infinite(beta), beta_upper, beta)
      
      # 5. Update Theta
      theta_upper <- min(x / beta^delta)
      theta       <- qgamma(pgamma(theta_upper, shape = t1 + alpha * n, rate = t2) * runif(1), 
                            shape = t1 + alpha * n, rate = t2)
      theta       <- ifelse(theta < .Machine$double.eps, min(x), theta)
    }
    
    alphas[r]   <- alpha
    betas[r]    <- beta
    deltas[, r] <- delta
    epsilons[r] <- pmax(epsilon, 10 * .Machine$double.eps)
    thetas[r]   <- theta
  }
  
  return(list(
    alphas   = alphas,
    betas    = betas,
    epsilons = epsilons,
    thetas   = thetas,
    deltas   = deltas
  ))
}

MCEM_app2 <- function(x, max.iter = 100, max.tol = 1e-2){
  diff        <- 10
  iter        <- 0
  prev_hyper_parms <- rep(100, 8)
  n           <- length(x)
  
  hyper_parms <- c(20, 2, quantile(x, 1 - 0.05) / min(x), 1, 1, 1, min(x), 1)
  pareto_parms <- c(1 / mean(log(x / min(x))), max(x) / min(x), 0.05, min(x))
  
  while(diff > max.tol & iter < max.iter){
    # E-step
    gibbs_sample <- mcmc_app2(x, hyper_parms, pareto_parms, R = 100, S = 100) 
    pareto_parms <- sapply(X = gibbs_sample, FUN = 'mean', simplify = TRUE)[1:4]
    
    # M-step
    alphas   <- gibbs_sample$alphas
    betas    <- gibbs_sample$betas
    epsilons <- gibbs_sample$epsilons
    thetas   <- gibbs_sample$thetas
    
    # Update prior for alpha
    s  <- log(mean(alphas)) - mean(log(alphas))
    a1 <- (3 - s + sqrt((s - 3)^2 + 24 * s)) / (12 * s)
    a2 <- a1 / mean(alphas)
    
    # Update prior for Beta
    b_star <- max(1.1, c(min(betas)))
    lambda <- 1 / mean(betas - b_star)
    
    # Update prior for epsilon
    ebar <- mean(epsilons)
    v    <- var(epsilons)
    b1   <- ebar * (ebar * (1 - ebar) / v - 1)
    b2   <- (1 - ebar) * (ebar * (1 - ebar) / v - 1)
    
    # Update prior for theta
    s  <- log(mean(thetas)) - mean(log(thetas))
    t1 <- (3 - s + sqrt((s - 3)^2 + 24 * s)) / (12 * s)
    t2 <- t1 / mean(thetas)
    
    hyper_parms <- c(a1, a2, b_star, lambda, b1, b2, t1, t2)
    iter        <- iter + 1
    diff        <- max(abs(log(hyper_parms) - log(prev_hyper_parms)))
    prev_hyper_parms <- hyper_parms
  }
  
  hyper_parms_list <- list(
    a1     = hyper_parms[1],
    a2     = hyper_parms[2],
    b_star = hyper_parms[3],
    lambda = hyper_parms[4],
    b1     = hyper_parms[5],
    b2     = hyper_parms[6],
    t1     = hyper_parms[7],
    t2     = hyper_parms[8]
  )
  return(list('Hyper.parameters' = hyper_parms_list, 'Gibbs.sample' = gibbs_sample))
}

ll_app2 <- function(log_parms, x) {
  alpha   <- exp(log_parms[1])
  beta    <- max(1.1, exp(log_parms[2]))
  epsilon <- plogis(log_parms[3])
  theta   <- min(x)
  ref  <- (1 - epsilon) * extraDistr::dpareto(x, a = alpha, b = theta)
  cont <- epsilon * extraDistr::dpareto(x, a = alpha, b = beta * theta)
  sum(log(ref + cont))
}

# ==============================================================================
# UI DEFINITION
# ==============================================================================
ui <- page_navbar(
  theme = app_theme,
  title = NULL,
  id = "nav",
  fillable = FALSE,
  
  header = tags$head(
    tags$link(
      rel = "stylesheet", 
      href = "https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap"
    ),
    tags$style(HTML("
      .sidebar .card, .sidebar .card-body {
        overflow: visible !important;
      }
      .selectize-dropdown {
        z-index: 1050 !important;
      }
      .compact-sidebar {
        max-height: calc(100vh - 100px) !important;
        overflow-y: auto !important;
        padding-right: 5px;
      }
      .compact-card {
        margin-bottom: 8px !important;
      }
      /* Center main navbar tab headings */
      .navbar-nav {
        justify-content: center !important;
        width: 100% !important;
      }
      /* Center nested tabsetPanel headings */
      .nav-tabs {
        justify-content: center !important;
      }
    "))
  ),
  
  # --------------------------------------------------------------------------
  # TAB 1: COVER PAGE
  # --------------------------------------------------------------------------
  nav_panel(
    title = "Cover Page",
    icon = icon("home"),
    
    fluidRow(
      column(10, offset = 1,
             div(
               style = "background-color: #ffffff; border: 1px solid #e2e8f0; border-left: 8px solid #002855; border-radius: 10px; padding: 30px; box-shadow: 0 4px 12px rgba(0,0,0,0.05); margin-top: 15px; margin-bottom: 25px;",
               
               h2("Identifying Outliers in the Norwegian Oil Dataset", 
                  style = "color: #002855; font-weight: 700; margin-bottom: 5px;"),
               h5("A Bayesian Pareto Scale-Inflated Mixture Framework", 
                  style = "color: #41748d; font-weight: 500; margin-bottom: 20px;"),
               
               hr(style = "border-top: 1px solid #e2e8f0; margin: 15px 0;"),
               
               div(
                 style = "font-size: 0.95rem; line-height: 1.6; color: #2c3e50; margin-bottom: 20px;",
                 fluidRow(
                   column(6, 
                          p(strong("Author: "), "Manfred Tristan Spiwak (u23665000)"),
                          p(strong("Degree: "), "BSc (Hons) Mathematical Statistics")
                   ),
                   column(6, 
                          p(strong("Supervisors: "), "Dr A Otto, Prof A Bekker & J Pillay"),
                          p(strong("Institution: "), "Department of Statistics, University of Pretoria")
                   )
                 )
               ),
               
               div(
                 style = "background-color: #f8fafc; padding: 15px 20px; border-radius: 8px; border: 1px solid #cbd5e1; margin-bottom: 25px;",
                 p(style = "margin: 0; font-size: 1rem; color: #1e293b; line-height: 1.6;",
                   strong("Overview: "), 
                   "This application provides an interactive suite for fitting and validating Bayesian Pareto scale-inflated mixture models. Designed to isolate extreme data contamination from heavy-tailed distributions, it leverages Gibbs sampling and Empirical Bayes (MCEM) hyperparameter tuning to deliver robust parameter estimates.")
               ),
               
               h4("Key Capabilities", style = "color: #002855; font-weight: 600; margin-bottom: 15px;"),
               
               layout_columns(
                 col_widths = c(6, 6, 6, 6),
                 card(
                   style = "border: 1px solid #e2e8f0; box-shadow: none;",
                   card_body(
                     style = "padding: 12px;",
                     h6(icon("database", class = "text-primary"), " Flexible Data Input", style = "font-weight: bold; margin-bottom: 5px;"),
                     p(style = "font-size: 0.85rem; color: #475569; margin: 0;", "Analyse the built-in 2023 Norwegian Oil dataset, generate synthetic scale-inflated samples, or upload custom Excel/CSV files.")
                   )
                 ),
                 card(
                   style = "border: 1px solid #e2e8f0; box-shadow: none;",
                   card_body(
                     style = "padding: 12px;",
                     h6(icon("gears", class = "text-primary"), " MCMC & MCEM Engines", style = "font-weight: bold; margin-bottom: 5px;"),
                     p(style = "font-size: 0.85rem; color: #475569; margin: 0;", "Run parallel or continuous Gibbs chains paired with Monte Carlo EM for Empirical Bayes prior estimation.")
                   )
                 ),
                 card(
                   style = "border: 1px solid #e2e8f0; box-shadow: none;",
                   card_body(
                     style = "padding: 12px;",
                     h6(icon("chart-line", class = "text-primary"), " Model Validation", style = "font-weight: bold; margin-bottom: 5px;"),
                     p(style = "font-size: 0.85rem; color: #475569; margin: 0;", "Evaluate model fidelity using Kolmogorov-Smirnov (KS) test statistics and empirical vs. predictive cumulative distribution plots.")
                   )
                 ),
                 card(
                   style = "border: 1px solid #e2e8f0; box-shadow: none;",
                   card_body(
                     style = "padding: 12px;",
                     h6(icon("bullseye", class = "text-primary"), " Outlier Profiling", style = "font-weight: bold; margin-bottom: 5px;"),
                     p(style = "font-size: 0.85rem; color: #475569; margin: 0;", "Isolate contaminated points probabilistically with posterior outlier classifications (\\(P(\\delta_i = 1)\\)) and interactive tables.")
                   )
                 )
               )
             )
      )
    )
  ),
  
  # --------------------------------------------------------------------------
  # TAB 2: METHODOLOGY
  # --------------------------------------------------------------------------
  nav_panel(
    title = "Methodology",
    icon = icon("graduation-cap"),
    
    withMathJax(),
    
    layout_columns(
      col_widths = c(6, 6, 12),
      
      # CARD 1: MODEL SPECIFICATION
      card(
        card_header(class = "bg-primary text-white", h6(icon("square-root-variable"), " 1. Scale-Inflated Pareto Mixture", style = "margin:0;")),
        card_body(
          style = "font-size: 0.9rem;",
          p("The model separates bulk observations from extreme scale-inflated contamination via a two-component mixture:"),
          p("$$f(x_i | \\alpha, \\theta, \\beta, \\epsilon) = (1 - \\epsilon) \\frac{\\alpha \\theta^\\alpha}{x_i^{\\alpha+1}} I(x_i \\ge \\theta) + \\epsilon \\frac{\\alpha (\\beta \\theta)^\\alpha}{x_i^{\\alpha+1}} I(x_i \\ge \\beta \\theta)$$"),
          p("Introducing a latent Bernoulli outlier indicator \\(\\delta_i \\sim \\text{Bernoulli}(\\epsilon)\\) yields the compact scale-inflated form:"),
          p("$$f(x_i | \\alpha, \\theta, \\beta, \\epsilon, \\delta_i) = \\frac{\\alpha (\\beta^{\\delta_i} \\theta)^\\alpha}{x_i^{\\alpha+1}} I\\left(x_i \\ge \\beta^{\\delta_i} \\theta\\right)$$")
        )
      ),
      
      # CARD 2: PRIORS & GIBBS SAMPLER
      card(
        card_header(class = "bg-primary text-white", h6(icon("diagram-project"), " 2. Priors & Full Conditional Posteriors", style = "margin:0;")),
        card_body(
          style = "padding: 10px; font-size: 0.85rem;",
          tags$table(
            class = "table table-sm table-striped align-middle",
            style = "margin-bottom: 0;",
            tags$thead(
              class = "table-light",
              tags$tr(
                tags$th(style = "width: 22%; font-weight: 600;", "Parameter"),
                tags$th(style = "width: 33%; font-weight: 600;", "Prior Distribution"),
                tags$th(style = "width: 45%; font-weight: 600;", "Full Conditional Posterior")
              )
            ),
            tags$tbody(
              tags$tr(
                tags$td(strong("Shape (\\(\\alpha\\))")),
                tags$td("\\(\\alpha \\sim \\text{Gamma}(a_1, a_2)\\)"),
                tags$td("\\(\\text{Gamma}\\left(a_1 + n, a_2 + \\sum_{i=1}^n \\ln x_i - n \\ln \\theta - k \\ln \\beta\\right)\\)")
              ),
              tags$tr(
                tags$td(strong("Threshold (\\(\\theta\\))")),
                tags$td("\\(\\theta \\sim \\text{Gamma}(t_1, t_2)\\)"),
                tags$td("Truncated \\(\\text{Gamma}(t_1 + \\alpha n, t_2)\\) on \\(0 \\le \\theta \\le \\min_i \\left(\\frac{x_i}{\\beta^{\\delta_i}}\\right)\\)")
              ),
              tags$tr(
                tags$td(strong("Contamination (\\(\\epsilon\\))")),
                tags$td("\\(\\epsilon \\sim \\text{Beta}(b_1, b_2)\\)"),
                tags$td("\\(\\text{Beta}(b_1 + k, b_2 + n - k)\\)")
              ),
              tags$tr(
                tags$td(strong("Inflation (\\(\\beta\\))")),
                tags$td("\\(\\beta \\sim \\text{Shifted-Exp}(\\beta^*, \\lambda)\\)"),
                tags$td("Truncated \\(\\text{Gamma}(\\alpha k + 1, \\lambda)\\) on \\(\\beta^* \\le \\beta \\le \\min_{i: \\delta_i=1}\\left(\\frac{x_i}{\\theta}\\right)\\)")
              ),
              tags$tr(
                tags$td(strong("Indicator (\\(\\delta_i\\))")),
                tags$td("\\(\\delta_i | \\epsilon \\sim \\text{Bernoulli}(\\epsilon)\\)"),
                tags$td("\\(P(\\delta_i = 1) = \\frac{\\beta^\\alpha \\epsilon}{1 - \\epsilon + \\beta^\\alpha \\epsilon}\\) for \\(x_i \\ge \\beta \\theta\\); \\(0\\) for \\(\\theta \\le x_i < \\beta \\theta\\)")
              )
            )
          )
        )
      ),
      
      # CARD 3: EMPIRICAL BAYES & DENSITY ESTIMATION
      card(
        card_header(class = "bg-primary text-white", h6(icon("gears"), " 3. Empirical Bayes (MCEM) & Diagnostics", style = "margin:0;")),
        card_body(
          style = "font-size: 0.88rem;",
          p(strong("Monte Carlo Expectation-Maximisation (MCEM):")),
          p("When prior parameters are unknown, Empirical Bayes tuning alternates between an ", strong("E-step"), " (sampling via Gibbs draws) and an ", strong("M-step"), " (updating hyperparameters \\(a_1, a_2, b^*, \\lambda, b_1, b_2, t_1, t_2\\) until convergence)."),
          hr(style = "margin: 8px 0;"),
          p(strong("Rao-Blackwellised Density Estimator:")),
          p("Reduces sampling variance by averaging conditional posterior distributions over MCMC draws \\(R\\):"),
          p("$$\\hat{f}_{\\text{RB}}(\\alpha) = \\frac{1}{R}\\sum_{r=1}^R \\text{Gamma}\\left(\\alpha; a_1 + n, a_2 + \\sum \\ln x_i - n\\ln \\theta^{(r)} - k^{(r)}\\ln \\beta^{(r)}\\right)$$")
        )
      )
    )
  ),
  
  # --------------------------------------------------------------------------
  # TAB 3: DATA & MCMC CONTROL
  # --------------------------------------------------------------------------
  nav_panel(
    title = "Data & Setup",
    icon = icon("sliders"),
    
    layout_sidebar(
      sidebar = sidebar(
        width = 360,
        class = "compact-sidebar",
        title = h6(icon("gear"), " Control Panel", style = "margin:0; color:#002855; font-weight:bold;"),
        
        card(
          class = "compact-card",
          card_header(h6("1. Data Input Strategy", style = "margin:0; font-size:0.85rem;")),
          card_body(
            style = "padding: 8px;",
            selectInput("data_source", "Select Data Source:",
                        choices = c("Norwegian Oil Dataset" = "norwegian",
                                    "Simulated Pareto Data" = "simulated",
                                    "Upload Excel File (.xlsx)" = "upload")),
            
            conditionalPanel(
              condition = "input.data_source == 'upload'",
              fileInput("file", "Upload Excel File (.xlsx)", accept = c(".xlsx")),
              textInput("col_name", "Column Name:", value = "x")
            ),
            
            conditionalPanel(
              condition = "input.data_source == 'simulated'",
              fluidRow(
                column(6, numericInput("n_sim", "Obs n (Std):", value = 160, min = 1)),
                column(6, numericInput("m_sim", "Obs m (Inf):", value = 40, min = 0))
              ),
              fluidRow(
                column(6, numericInput("alpha_true", "True Alpha:", value = 2.5, step = 0.1)),
                column(6, numericInput("beta_true", "True Beta:", value = 3.0, step = 0.1))
              ),
              actionButton("gen_sim", "Regenerate Synthetic Data", class = "btn-outline-primary btn-sm w-100 mt-1", icon = icon("rotate"))
            )
          )
        ),
        
        card(
          class = "compact-card",
          card_header(h6("2. MCMC Sampling Settings", style = "margin:0; font-size:0.85rem;")),
          card_body(
            style = "padding: 8px;",
            checkboxInput("use_rs", "Use R x S Parallel Chains Scheme?", value = TRUE),
            
            conditionalPanel(
              condition = "input.use_rs == true",
              fluidRow(
                column(6, numericInput("R_chains", "Chains (R):", value = 100, min = 10, step = 50)),
                column(6, numericInput("S_steps", "Steps/Chain (S):", value = 100, min = 1, step = 5))
              )
            ),
            
            conditionalPanel(
              condition = "input.use_rs == false",
              fluidRow(
                column(6, numericInput("N_iter", "Total Steps:", value = 200000, min = 100, step = 500)),
                column(6, numericInput("burn_in", "Burn-in Steps:", value = 10000, min = 0, step = 100))
              )
            ),
            
            numericInput("beta_star", "Beta Star (beta_star):", value = 1000, min = 1.0, step = 0.1)
          )
        ),
        
        card(
          class = "compact-card",
          card_header(h6("3. Hyperparameters & Settings", style = "margin:0; font-size:0.85rem;")),
          card_body(
            style = "padding: 8px;",
            checkboxInput("use_rb", "Use Rao-Blackwellized Densities (α & θ)?", value = TRUE),
            checkboxInput("fix_theta", "Fix Theta Threshold (θ)?", value = FALSE),
            numericInput("theta_input", "Theta (Value or Initial):", value = 1000000, min = 1, step = 1000),
            
            conditionalPanel(
              condition = "!input.fix_theta",
              fluidRow(
                column(6, numericInput("t1", "Theta t1:", value = 10)),
                column(6, numericInput("t2", "Theta t2:", value = 0.0001))
              )
            ),
            
            fluidRow(
              column(6, numericInput("a1", "Alpha a1:", value = 1000)),
              column(6, numericInput("a2", "Alpha a2:", value = 1000))
            ),
            fluidRow(
              column(6, numericInput("b1", "Eps b1:", value = 0.1842)),
              column(6, numericInput("b2", "Eps b2:", value = 3.5))
            ),
            numericInput("lambda_prior", "Beta Prior lambda:", value = 0.001)
          )
        ),
        
        actionButton("run_mcmc", "Run Gibbs Sampler", class = "btn-primary w-100 mt-1 mb-2", icon = icon("play"))
      ),
      
      card(
        card_header(class = "bg-light", h6(icon("chart-column"), " Active Data Overview", style = "margin:0;")),
        card_body(
          style = "padding: 15px;",
          plotOutput("data_hist", height = "220px"),
          hr(style = "margin: 10px 0;"),
          verbatimTextOutput("data_summary")
        )
      )
    )
  ),
  
  # --------------------------------------------------------------------------
  # TAB 4: MODEL 1 RESULTS
  # --------------------------------------------------------------------------
  nav_panel(
    title = "Model (With Outliers)",
    icon = icon("chart-line"),
    
    div(style = "margin-bottom: 10px;", uiOutput("m1_metrics_ui")),
    
    layout_columns(
      col_widths = c(6, 6),
      card(
        card_header(class = "bg-primary text-white", style = "padding: 8px 12px;", h6(icon("chart-area"), " Posterior Distributions", style = "margin:0;")),
        card_body(style = "padding: 8px;", plotOutput("post_plots_m1", height = "350px"))
      ),
      card(
        card_header(class = "bg-secondary text-white", style = "padding: 8px 12px;", h6(icon("wave-square"), " MCMC Convergence (ACF)", style = "margin:0;")),
        card_body(style = "padding: 8px;", plotOutput("acf_plots_m1", height = "350px"))
      )
    )
  ),
  
  # --------------------------------------------------------------------------
  # TAB 5: KS TEST EVALUATION
  # --------------------------------------------------------------------------
  nav_panel(
    title = "KS Test Evaluation",
    icon = icon("check-double"),
    
    layout_columns(
      col_widths = c(6, 6),
      
      div(
        uiOutput("ks_output")
      ),
      
      card(
        card_header(class = "bg-primary text-white", style = "padding: 8px 12px;", 
                    h6(icon("chart-line"), " Empirical vs Model Cumulative Distribution (eCDF)", style = "margin:0;")),
        card_body(
          style = "padding: 12px;",
          plotOutput("ks_ecdf_plot", height = "430px")
        )
      )
    )
  ),
  
  # --------------------------------------------------------------------------
  # TAB 6: COMPARATIVE ANALYSIS
  # --------------------------------------------------------------------------
  nav_panel(
    title = "Outliers vs. No Outliers",
    icon = icon("code-compare"),
    
    layout_columns(
      col_widths = c(4, 8),
      
      card(
        card_header(class = "bg-primary text-white", h5(icon("table"), " Parameter Comparison Table", style="margin:0;")),
        card_body(
          tableOutput("comp_estimates_tab"),
          hr(),
          p(class = "text-muted", style = "font-size:0.88rem;",
            "Model 1 separates scale-inflated observations using ", span("\\(\\beta > 1\\)", style="color:#002855; font-weight:bold;"), 
            " while Model 2 forces all points into a single Pareto curve (", span("\\(\\beta = 1\\)", style="color:#e63946; font-weight:bold;"), ").")
        )
      ),
      
      div(
        layout_columns(
          col_widths = c(6, 6),
          card(
            card_header(class = "bg-primary text-white", h6(icon("chart-line"), " Alpha Parameter Shift", style="margin:0;")),
            card_body(plotOutput("alpha_comparison_plot", height = "220px"))
          ),
          card(
            card_header(class = "bg-primary text-white", h6(icon("chart-line"), " Theta Parameter Shift", style="margin:0;")),
            card_body(plotOutput("theta_comparison_plot", height = "220px"))
          )
        ),
        card(
          style = "margin-top: 12px;",
          card_header(class = "bg-primary text-white", h6(icon("chart-bar"), " Posterior Outlier Probabilities vs Observation Values", style="margin:0;")),
          card_body(plotOutput("outlier_probs_plot", height = "220px"))
        )
      )
    )
  ),
  
  # --------------------------------------------------------------------------
  # TAB 7: INDEPENDENT MCEM ANALYSIS MODULE
  # --------------------------------------------------------------------------
  nav_panel(
    title = " MCEM Outlier Analysis",
    icon = icon("calculator"),
    
    div(
      style = "padding: 10px 0;",
      h3("Contaminated Pareto Mixture Modeling & Outlier Estimation (MCEM Algorithm)", 
         style = "color: #002855; font-weight: 700; margin-bottom: 20px;"),
      
      sidebarLayout(
        sidebarPanel(
          width = 3,
          selectInput("app2_data_source", "Data Input Source", 
                      choices = c("Norwegian Oil Dataset", 
                                  "Simulated Pareto Mixture", 
                                  "Upload Dataset (.xlsx / .csv)")),
          
          conditionalPanel(
            condition = "input.app2_data_source == 'Simulated Pareto Mixture'",
            tags$h4("Simulation Parameters"),
            numericInput("app2_sim_n", "Sample Size (n):", value = 128, min = 10, step = 10),
            numericInput("app2_sim_alpha", "Alpha (Shape):", value = 1.0, min = 0.1, step = 0.1),
            numericInput("app2_sim_theta", "Theta (Scale):", value = 125500, min = 1),
            numericInput("app2_sim_beta", "Beta (Contamination multiplier):", value = 2800, min = 1.1),
            numericInput("app2_sim_epsilon", "Epsilon (Outlier proportion):", value = 0.05, min = 0.001, max = 0.5, step = 0.01)
          ),
          
          conditionalPanel(
            condition = "input.app2_data_source == 'Upload Dataset (.xlsx / .csv)'",
            tags$h4("Upload Dataset"),
            fileInput("app2_file_upload", "Select File", accept = c(".xlsx", ".xls", ".csv")),
            uiOutput("app2_column_select_ui")
          ),
          
          tags$hr(),
          tags$h4("MCEM Controls"),
          numericInput("app2_max_tol", "Convergence Tolerance (max.tol):", value = 0.01, min = 0.0001, step = 0.005),
          numericInput("app2_max_iter", "Max Iterations:", value = 50, min = 5, step = 5),
          
          tags$br(),
          actionButton("app2_run_btn", "Execute Analysis", class = "btn-primary btn-block w-100")
        ),
        
        mainPanel(
          width = 9,
          tabsetPanel(
            tabPanel("Summary & Hyperparameters", 
                     tags$br(),
                     wellPanel(
                       tags$h4("Parameter Estimates Comparison", style = "margin-top: 0; font-weight: bold; color: #2c3e50;"),
                       tableOutput("app2_summary_table")
                     ),
                     tags$br(),
                     wellPanel(
                       tags$h4("Final Estimated Prior Hyperparameters", style = "margin-top: 0; font-weight: bold; color: #2c3e50;"),
                       tags$p("M-step converged prior hyperparameters for empirical Bayes inference:", style = "color: #7f8c8d; font-size: 13px;"),
                       tableOutput("app2_hyperparameters_table")
                     )
            ),
            tabPanel("Posterior Distributions", 
                     plotOutput("app2_plot_histograms", height = "650px")),
            tabPanel("Outlier Probability Plot", 
                     plotOutput("app2_plot_outliers", height = "500px")),
            tabPanel("Observations & Posterior Probabilities", 
                     DTOutput("app2_data_table"))
          )
        )
      )
    )
  )
)

# ==============================================================================
# SERVER DEFINITION
# ==============================================================================
server <- function(input, output, session) {
  
  # --------------------------------------------------------------------------
  # SERVER LOGIC - APP 1
  # --------------------------------------------------------------------------
  
  active_dataset <- reactiveVal(norwegian_oil_data)
  sim_seed <- reactiveVal(123)
  
  observeEvent(c(input$gen_sim, input$data_source), {
    if (input$data_source == "simulated") {
      set.seed(sim_seed())
      n <- input$n_sim
      m <- input$m_sim
      alpha_val <- input$alpha_true
      theta_val <- input$theta_input
      beta_val  <- input$beta_true
      
      sim_val <- rpareto_custom(n, alpha_val, theta_val)
      inflated_sim_val <- scale_inflated_rpareto(m, alpha_val, theta_val, beta_val)
      active_dataset(round(c(sim_val, inflated_sim_val), 2))
    } else if (input$data_source == "norwegian") {
      active_dataset(norwegian_oil_data)
    }
  }, ignoreInit = FALSE)
  
  observeEvent(input$file, {
    req(input$file)
    df <- read_excel(input$file$datapath, progress = FALSE)
    raw_vec <- if (input$col_name %in% colnames(df)) df[[input$col_name]] else df[[1]]
    
    clean_vec <- as.numeric(unname(raw_vec))
    clean_vec <- clean_vec[!is.na(clean_vec) & clean_vec > 0]
    
    if (length(clean_vec) > 0) {
      active_dataset(clean_vec)
    }
  })
  
  observeEvent(input$gen_sim, {
    sim_seed(sample.int(1e6, 1))
  })
  
  comp_sim_val_data <- reactive({
    active_dataset()
  })
  
  output$data_hist <- renderPlot({
    data_val <- comp_sim_val_data()
    par(mar = c(3, 3, 1, 1))
    hist(data_val, breaks = 30, main = "", xlab = "Observations", freq = FALSE, col = "#d0e1fd", border = "#002855")
    lines(density(data_val), col = "#e63946", lwd = 2.5)
  })
  
  output$data_summary <- renderPrint({
    val <- comp_sim_val_data()
    cat("Summary of Active Data:\n")
    print(summary(val))
    cat("\nSample Size (n):", length(val))
  })
  
  mcmc_results <- eventReactive(input$run_mcmc, {
    comp_sim_val <- comp_sim_val_data()
    n2 <- length(comp_sim_val)
    
    use_rs <- input$use_rs
    
    if (use_rs) {
      R_chains <- input$R_chains
      S_steps  <- input$S_steps
      num_samples <- R_chains
    } else {
      N_iter  <- input$N_iter
      burn_in <- input$burn_in
      num_samples <- max(1, N_iter - burn_in)
    }
    
    a_1 <- input$a1
    a_2 <- input$a2
    b_1 <- input$b1
    b_2 <- input$b2
    lambda <- input$lambda_prior
    beta_star <- input$beta_star
    
    fix_theta <- input$fix_theta
    theta_fixed <- min(input$theta_input, min(comp_sim_val) * 0.9999)
    t_1 <- if (!fix_theta) input$t1 else 10
    t_2 <- if (!fix_theta) input$t2 else 0.0001
    use_rb <- input$use_rb
    
    alphas   <- numeric(num_samples) 
    thetas   <- numeric(num_samples) 
    epsilons <- numeric(num_samples) 
    betas    <- numeric(num_samples) 
    deltas   <- matrix(0, nrow = num_samples, ncol = n2) 
    p_mat    <- matrix(0, nrow = num_samples, ncol = n2)
    
    if (use_rs) {
      withProgress(message = 'Sampling Model 1 (Parallel R x S Chains)...', value = 0, {
        for (r in 1:R_chains) {
          theta_ini   <- if (fix_theta) theta_fixed else runif(1, min(comp_sim_val) * 0.8, min(comp_sim_val) * 0.999)
          alpha_ini   <- runif(1, 0.5, 5.0)
          epsilon_ini <- rbeta(1, b_1, b_2)
          beta_ini    <- runif(1, 1.1, 10.0)
          
          for (s in 1:S_steps) {
            prob_inflate <- ((beta_ini^alpha_ini) * epsilon_ini) / (1 - epsilon_ini + (beta_ini^alpha_ini) * epsilon_ini)
            prob_vec_inflated <- (comp_sim_val >= theta_ini * beta_ini) * prob_inflate
            delta_ini <- rbinom(n2, 1, prob_vec_inflated)
            
            k <- sum(delta_ini) 
            epsilon_ini <- rbeta(1, b_1 + k, b_2 + n2 - k) 
            
            shape_a <- a_1 + n2
            rate_a  <- a_2 + sum(log(comp_sim_val)) - n2 * log(theta_ini) - k * log(beta_ini)
            alpha_ini <- rgamma(1, shape = shape_a, rate = rate_a) 
            
            if (fix_theta) {
              theta_ini <- theta_fixed
            } else {
              theta_ub <- min(c(comp_sim_val[delta_ini == 0], comp_sim_val[delta_ini == 1] / beta_ini))
              post_shape_t <- t_1 + n2 * alpha_ini
              post_rate_t  <- t_2
              max_theta_ub_log <- pgamma(theta_ub, post_shape_t, post_rate_t, log.p = TRUE)
              unif_theta_log   <- max_theta_ub_log + log(runif(1))
              theta_ini        <- qgamma(unif_theta_log, post_shape_t, post_rate_t, log.p = TRUE)
            }
            
            if (k > 0) {
              beta_ub <- max(beta_star + 1e-5, min(comp_sim_val[delta_ini == 1] / theta_ini))
              post_shape_b  <- alpha_ini * k + 1
              post_rate_b   <- lambda
              min_prob_beta <- pgamma(beta_star, post_shape_b, post_rate_b)
              max_prob_beta <- pgamma(beta_ub, post_shape_b, post_rate_b)
              unif_beta     <- runif(1, min_prob_beta, max_prob_beta)
              beta_ini      <- qgamma(unif_beta, post_shape_b, post_rate_b)
            } else {
              beta_ini <- beta_star + rexp(1, lambda)
            }
          }
          
          alphas[r]   <- alpha_ini 
          thetas[r]   <- theta_ini 
          epsilons[r] <- epsilon_ini 
          betas[r]    <- beta_ini 
          deltas[r, ] <- delta_ini 
          p_mat[r, ]  <- prob_vec_inflated
          
          if (r %% 20 == 0) incProgress(0.5 * (r / R_chains))
        }
      })
    } else {
      theta_ini   <- if (fix_theta) theta_fixed else runif(1, min(comp_sim_val) * 0.8, min(comp_sim_val) * 0.999)
      alpha_ini   <- runif(1, 0.5, 5.0)
      epsilon_ini <- rbeta(1, b_1, b_2)
      beta_ini    <- runif(1, 1.1, 10.0)
      
      idx <- 1
      withProgress(message = 'Sampling Model 1 (Single Continuous Chain)...', value = 0, {
        for (i in 1:N_iter) {
          prob_inflate <- ((beta_ini^alpha_ini) * epsilon_ini) / (1 - epsilon_ini + (beta_ini^alpha_ini) * epsilon_ini)
          prob_vec_inflated <- (comp_sim_val >= theta_ini * beta_ini) * prob_inflate
          delta_ini <- rbinom(n2, 1, prob_vec_inflated)
          
          k <- sum(delta_ini) 
          epsilon_ini <- rbeta(1, b_1 + k, b_2 + n2 - k) 
          
          shape_a <- a_1 + n2
          rate_a  <- a_2 + sum(log(comp_sim_val)) - n2 * log(theta_ini) - k * log(beta_ini)
          alpha_ini <- rgamma(1, shape = shape_a, rate = rate_a) 
          
          if (fix_theta) {
            theta_ini <- theta_fixed
          } else {
            theta_ub <- min(c(comp_sim_val[delta_ini == 0], comp_sim_val[delta_ini == 1] / beta_ini))
            post_shape_t <- t_1 + n2 * alpha_ini
            post_rate_t  <- t_2
            max_theta_ub_log <- pgamma(theta_ub, post_shape_t, post_rate_t, log.p = TRUE)
            unif_theta_log   <- max_theta_ub_log + log(runif(1))
            theta_ini        <- qgamma(unif_theta_log, post_shape_t, post_rate_t, log.p = TRUE)
          }
          
          if (k > 0) {
            beta_ub <- max(beta_star + 1e-5, min(comp_sim_val[delta_ini == 1] / theta_ini))
            post_shape_b  <- alpha_ini * k + 1
            post_rate_b   <- lambda
            min_prob_beta <- pgamma(beta_star, post_shape_b, post_rate_b)
            max_prob_beta <- pgamma(beta_ub, post_shape_b, post_rate_b)
            unif_beta     <- runif(1, min_prob_beta, max_prob_beta)
            beta_ini      <- qgamma(unif_beta, post_shape_b, post_rate_b)
          } else {
            beta_ini <- beta_star + rexp(1, lambda)
          }
          
          if (i > burn_in) {
            alphas[idx]   <- alpha_ini
            thetas[idx]   <- theta_ini
            epsilons[idx] <- epsilon_ini
            betas[idx]    <- beta_ini
            deltas[idx, ] <- delta_ini
            p_mat[idx, ]  <- prob_vec_inflated
            idx <- idx + 1
          }
          
          if (i %% 500 == 0) incProgress(0.5 * (i / N_iter))
        }
      })
    }
    
    alpha_grid <- NULL
    rb_alpha_density <- NULL
    theta_grid <- NULL
    rb_theta_density <- NULL
    
    if (use_rb) {
      alpha_grid <- seq(min(alphas) * 0.8, max(alphas) * 1.2, length.out = 300)
      rb_alpha_density <- numeric(length(alpha_grid))
      for (r in 1:num_samples) {
        k_r <- sum(deltas[r, ])
        shape_r <- a_1 + n2
        rate_r  <- a_2 + sum(log(comp_sim_val)) - n2 * log(thetas[r]) - k_r * log(betas[r])
        rb_alpha_density <- rb_alpha_density + dgamma(alpha_grid, shape = shape_r, rate = rate_r) / num_samples
      }
      
      if (!fix_theta) {
        t_min <- min(thetas)
        t_max <- max(thetas)
        if (t_min == t_max) {
          t_min <- t_min * 0.99
          t_max <- t_max * 1.01
        } else {
          t_min <- t_min * 0.95
          t_max <- t_max * 1.05
        }
        
        theta_grid <- seq(t_min, t_max, length.out = 300)
        rb_theta_density <- numeric(length(theta_grid))
        
        for (r in 1:num_samples) {
          delta_r <- deltas[r, ]
          beta_r  <- betas[r]
          theta_ub_r <- min(c(comp_sim_val[delta_r == 0], comp_sim_val[delta_r == 1] / beta_r))
          post_shape_t <- t_1 + n2 * alphas[r]
          post_rate_t  <- t_2
          
          log_prob_trunc <- pgamma(theta_ub_r, shape = post_shape_t, rate = post_rate_t, log.p = TRUE)
          valid_idx <- theta_grid <= theta_ub_r
          
          if (any(valid_idx) && is.finite(log_prob_trunc)) {
            log_dens <- dgamma(theta_grid[valid_idx], shape = post_shape_t, rate = post_rate_t, log = TRUE) - log_prob_trunc
            rb_theta_density[valid_idx] <- rb_theta_density[valid_idx] + exp(log_dens) / num_samples
          }
        }
      }
    }
    
    alphas_2   <- numeric(num_samples)
    thetas_2   <- numeric(num_samples)
    epsilons_2 <- numeric(num_samples)
    betas_2    <- numeric(num_samples)
    
    if (use_rs) {
      withProgress(message = 'Sampling Model 2 (Parallel R x S Chains)...', value = 0.5, {
        for (r in 1:R_chains) {
          theta_ini_2   <- if (fix_theta) theta_fixed else min(comp_sim_val) * 0.999
          alpha_ini_2   <- n2 / max(1e-8, (sum(log(comp_sim_val)) - n2 * log(theta_ini_2)))
          epsilon_ini_2 <- 0.05
          beta_ini_2    <- 1
          
          for (s in 1:S_steps) {
            k_2 <- 0
            epsilon_ini_2 <- rbeta(1, b_1 + k_2, b_2 + n2 - k_2)
            alpha_ini_2   <- rgamma(1, a_1 + n2, a_2 + sum(log(comp_sim_val)) - n2 * log(theta_ini_2))
            
            if (fix_theta) {
              theta_ini_2 <- theta_fixed
            } else {
              theta_ub_2 <- min(comp_sim_val)
              post_shape_t_2 <- t_1 + n2 * alpha_ini_2
              post_rate_t_2  <- t_2
              max_theta_ub_log_2 <- pgamma(theta_ub_2, post_shape_t_2, post_rate_t_2, log.p = TRUE)
              unif_theta_log_2   <- max_theta_ub_log_2 + log(runif(1))
              theta_ini_2        <- qgamma(unif_theta_log_2, post_shape_t_2, post_rate_t_2, log.p = TRUE)
            }
          }
          
          alphas_2[r]   <- alpha_ini_2
          thetas_2[r]   <- theta_ini_2
          epsilons_2[r] <- epsilon_ini_2
          betas_2[r]    <- 1
          
          if (r %% 20 == 0) incProgress(0.5 * (r / R_chains))
        }
      })
    } else {
      theta_ini_2   <- if (fix_theta) theta_fixed else min(comp_sim_val) * 0.999
      alpha_ini_2   <- n2 / max(1e-8, (sum(log(comp_sim_val)) - n2 * log(theta_ini_2)))
      epsilon_ini_2 <- 0.05
      beta_ini_2    <- 1
      
      idx2 <- 1
      withProgress(message = 'Sampling Model 2 (Single Continuous Chain)...', value = 0.5, {
        for (i in 1:N_iter) {
          k_2 <- 0
          epsilon_ini_2 <- rbeta(1, b_1 + k_2, b_2 + n2 - k_2)
          alpha_ini_2   <- rgamma(1, a_1 + n2, a_2 + sum(log(comp_sim_val)) - n2 * log(theta_ini_2))
          
          if (fix_theta) {
            theta_ini_2 <- theta_fixed
          } else {
            theta_ub_2 <- min(comp_sim_val)
            post_shape_t_2 <- t_1 + n2 * alpha_ini_2
            post_rate_t_2  <- t_2
            max_theta_ub_log_2 <- pgamma(theta_ub_2, post_shape_t_2, post_rate_t_2, log.p = TRUE)
            unif_theta_log_2   <- max_theta_ub_log_2 + log(runif(1))
            theta_ini_2        <- qgamma(unif_theta_log_2, post_shape_t_2, post_rate_t_2, log.p = TRUE)
          }
          
          if (i > burn_in) {
            alphas_2[idx2]   <- alpha_ini_2
            thetas_2[idx2]   <- theta_ini_2
            epsilons_2[idx2] <- epsilon_ini_2
            betas_2[idx2]    <- 1
            idx2 <- idx2 + 1
          }
          
          if (i %% 500 == 0) incProgress(0.5 * (i / N_iter))
        }
      })
    }
    
    list(
      alphas = alphas, thetas = thetas, epsilons = epsilons, betas = betas, deltas = deltas, p_mat = p_mat,
      use_rb = use_rb, alpha_grid = alpha_grid, rb_alpha_density = rb_alpha_density,
      theta_grid = theta_grid, rb_theta_density = rb_theta_density,
      alphas_2 = alphas_2, thetas_2 = thetas_2, epsilons_2 = epsilons_2, betas_2 = betas_2,
      est_alpha = mean(alphas), est_theta = mean(thetas), est_epsilon = mean(epsilons), est_beta = mean(betas),
      est_alpha_2 = mean(alphas_2), est_theta_2 = mean(thetas_2), est_epsilon_2 = mean(epsilons_2), est_beta_2 = mean(betas_2),
      comp_sim_val = comp_sim_val,
      is_theta_fixed = fix_theta,
      is_simulated = (input$data_source == "simulated"),
      true_alpha = input$alpha_true,
      true_beta = input$beta_true,
      true_theta = input$theta_input,
      true_epsilon = input$m_sim / max(1, (input$n_sim + input$m_sim))
    )
  })
  
  output$m1_metrics_ui <- renderUI({
    res <- mcmc_results()
    
    alpha_sub <- if (res$is_simulated) span(icon("bullseye"), paste0(" True: ", res$true_alpha), style = "font-size: 0.8rem; opacity: 0.9; font-weight: 600;") else NULL
    beta_sub  <- if (res$is_simulated) span(icon("bullseye"), paste0(" True: ", res$true_beta), style = "font-size: 0.8rem; opacity: 0.9; font-weight: 600;") else NULL
    eps_sub   <- if (res$is_simulated) span(icon("bullseye"), paste0(" True: ", round(res$true_epsilon, 4)), style = "font-size: 0.8rem; opacity: 0.9; font-weight: 600;") else NULL
    theta_sub <- if (res$is_simulated) span(icon("bullseye"), paste0(" True: ", format(res$true_theta, big.mark = ",")), style = "font-size: 0.8rem; opacity: 0.9; font-weight: 600;") else NULL
    
    layout_columns(
      col_widths = c(3, 3, 3, 3),
      value_box(
        title = "Posterior Alpha",
        value = round(res$est_alpha, 4),
        alpha_sub,
        showcase = icon("chart-area"),
        theme = "primary"
      ),
      value_box(
        title = "Posterior Beta",
        value = round(res$est_beta, 4),
        beta_sub,
        showcase = icon("up-right-and-down-left-from-center"),
        theme = "info"
      ),
      value_box(
        title = "Posterior Epsilon",
        value = round(res$est_epsilon, 4),
        eps_sub,
        showcase = icon("percent"),
        theme = "warning"
      ),
      value_box(
        title = "Posterior Theta",
        value = format(round(res$est_theta, 2), big.mark = ","),
        theta_sub,
        showcase = icon("ruler-horizontal"),
        theme = "success"
      )
    )
  })
  
  output$post_plots_m1 <- renderPlot({
    res <- mcmc_results()
    par(mfrow = c(2, 2), mar = c(4.5, 4.5, 1.5, 1), mgp = c(2.4, 0.7, 0))
    
    hist(res$alphas, breaks = 20, probability = TRUE, col = "#bfdbfe", border = "#1d4ed8", main = "Posterior Alpha", xlab = "Alpha")
    if (isTRUE(res$use_rb) && !is.null(res$rb_alpha_density)) {
      lines(res$alpha_grid, res$rb_alpha_density, col = "#1e3a8a", lwd = 2)
    } else {
      lines(density(res$alphas), col = "#1e3a8a", lwd = 2)
    }
    abline(v = res$est_alpha, col = "#e63946", lty = 2, lwd = 2)
    if (res$is_simulated) abline(v = res$true_alpha, col = "#000000", lty = 3, lwd = 2)
    
    if (res$is_theta_fixed) {
      plot(res$thetas, type = "l", col = "#065f46", lwd = 2, main = expression(paste("Fixed Theta (", theta, ") Constant")),
           ylab = "Theta", xlab = "Chain Index", ylim = c(res$est_theta * 0.9, res$est_theta * 1.1))
      abline(h = res$est_theta, col = "#e63946", lty = 2, lwd = 2)
      if (res$is_simulated) abline(h = res$true_theta, col = "#000000", lty = 3, lwd = 2)
    } else {
      hist(res$thetas, breaks = 30, probability = TRUE, col = "#a7f3d0", border = "#047857", main = "Posterior Theta", xlab = "Theta")
      if (isTRUE(res$use_rb) && !is.null(res$rb_theta_density)) {
        lines(res$theta_grid, res$rb_theta_density, col = "#065f46", lwd = 2)
      } else {
        lines(density(res$thetas), col = "#065f46", lwd = 2)
      }
      abline(v = res$est_theta, col = "#e63946", lty = 2, lwd = 2)
      if (res$is_simulated) abline(v = res$true_theta, col = "#000000", lty = 3, lwd = 2)
    }
    
    hist(res$epsilons, breaks = 20, probability = TRUE, col = "#fef08a", border = "#a16207", main = "Posterior Epsilon", xlab = "Epsilon")
    lines(density(res$epsilons), col = "#713f12", lwd = 2)
    abline(v = res$est_epsilon, col = "#e63946", lty = 2, lwd = 2)
    if (res$is_simulated) abline(v = res$true_epsilon, col = "#000000", lty = 3, lwd = 2)
    
    hist(res$betas, breaks = 20, probability = TRUE, col = "#fbcfe8", border = "#be185d", main = "Posterior Beta", xlab = "Beta")
    lines(density(res$betas), col = "#831843", lwd = 2)
    abline(v = res$est_beta, col = "#e63946", lty = 2, lwd = 2)
    if (res$is_simulated) abline(v = res$true_beta, col = "#000000", lty = 3, lwd = 2)
  })
  
  output$acf_plots_m1 <- renderPlot({
    res <- mcmc_results()
    par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1), mgp = c(2.4, 0.7, 0), cex.main = 1.0)
    acf(res$alphas, main = "ACF Alpha", col = "#1e3a8a")
    
    if (res$is_theta_fixed) {
      plot(1, type = "n", xlab = "", ylab = "", main = "ACF Theta (Fixed Constant)", axes = FALSE)
      text(1, 1, paste0("Theta fixed at\n", format(res$est_theta, big.mark = ",")), cex = 1.1, col = "#065f46")
    } else {
      acf(res$thetas, main = "ACF Theta", col = "#065f46")
    }
    
    acf(res$epsilons, main = "ACF Epsilon", col = "#713f12")
    acf(res$betas, main = "ACF Beta", col = "#831843")
  })
  
  ks_test_data <- reactive({
    res <- mcmc_results()
    comp_sim_val <- res$comp_sim_val
    sim_number <- max(1000, length(comp_sim_val) * 10)
    
    set.seed(123)
    is_inflated <- rbinom(sim_number, 1, res$est_epsilon)
    
    n_inf <- sum(is_inflated == 1)
    n_std <- sim_number - n_inf
    
    sim_val_new <- numeric(sim_number)
    if (n_inf > 0) {
      sim_val_new[is_inflated == 1] <- scale_inflated_rpareto(n_inf, res$est_alpha, res$est_theta, res$est_beta)
    }
    if (n_std > 0) {
      sim_val_new[is_inflated == 0] <- rpareto_custom(n_std, res$est_alpha, res$est_theta)
    }
    
    KS_test <- ks.test(comp_sim_val, sim_val_new, alternative = "greater")
    
    list(
      obs = comp_sim_val,
      sim = sim_val_new,
      ks = KS_test,
      stat_val = round(unname(KS_test$statistic), 4),
      p_val = KS_test$p.value,
      n = length(comp_sim_val)
    )
  })
  
  output$ks_output <- renderUI({
    ks_res <- ks_test_data()
    stat_val <- ks_res$stat_val
    p_val <- ks_res$p_val
    p_val_fmt <- if (p_val < 0.0001) "< 0.0001" else sprintf("%.4f", p_val)
    is_good_fit <- p_val > 0.05
    
    tagList(
      if (is_good_fit) {
        div(
          style = "background: linear-gradient(135deg, #2a9d8f 0%, #1e7268 100%); color: white; padding: 18px 20px; border-radius: 10px; margin-bottom: 16px; box-shadow: 0 4px 6px -1px rgba(0,0,0,0.1);",
          div(style = "display: flex; align-items: center; justify-content: space-between;",
              div(
                h5(icon("circle-check"), " Excellent Model Fit", style = "margin: 0; font-weight: 700; letter-spacing: 0.5px;"),
                p("Fail to reject H₀ at α = 0.05 level. Model accurately matches observed data.", style = "margin: 4px 0 0 0; font-size: 0.88rem; opacity: 0.95;")
              ),
              span(style = "background: rgba(255,255,255,0.22); padding: 6px 14px; border-radius: 20px; font-weight: bold; font-size: 0.85rem;", "p > 0.05")
          )
        )
      } else {
        div(
          style = "background: linear-gradient(135deg, #e63946 0%, #b81d28 100%); color: white; padding: 18px 20px; border-radius: 10px; margin-bottom: 16px; box-shadow: 0 4px 6px -1px rgba(0,0,0,0.1);",
          div(style = "display: flex; align-items: center; justify-content: space-between;",
              div(
                h5(icon("triangle-exclamation"), " Significant Divergence", style = "margin: 0; font-weight: 700; letter-spacing: 0.5px;"),
                p("Reject H₀ at α = 0.05 level. Significant difference detected.", style = "margin: 4px 0 0 0; font-size: 0.88rem; opacity: 0.95;")
              ),
              span(style = "background: rgba(255,255,255,0.22); padding: 6px 14px; border-radius: 20px; font-weight: bold; font-size: 0.85rem;", "p ≤ 0.05")
          )
        )
      },
      
      div(
        style = "display: grid; grid-template-columns: repeat(3, 1fr); gap: 10px; margin-bottom: 16px;",
        
        div(
          style = "background: #ffffff; padding: 14px 10px; border-radius: 8px; border: 1px solid #e2e8f0; border-top: 4px solid #002855; box-shadow: 0 1px 3px rgba(0,0,0,0.05); text-align: center;",
          span(style = "display: block; font-size: 0.72rem; font-weight: 700; color: #64748b; text-transform: uppercase; letter-spacing: 0.5px;", "KS Distance (D)"),
          span(style = "font-size: 1.35rem; font-weight: 800; color: #002855;", stat_val)
        ),
        
        div(
          style = "background: #ffffff; padding: 14px 10px; border-radius: 8px; border: 1px solid #e2e8f0; border-top: 4px solid #2a9d8f; box-shadow: 0 1px 3px rgba(0,0,0,0.05); text-align: center;",
          span(style = "display: block; font-size: 0.72rem; font-weight: 700; color: #64748b; text-transform: uppercase; letter-spacing: 0.5px;", "p-Value"),
          span(style = "font-size: 1.35rem; font-weight: 800; color: #2a9d8f;", p_val_fmt)
        ),
        
        div(
          style = "background: #ffffff; padding: 14px 10px; border-radius: 8px; border: 1px solid #e2e8f0; border-top: 4px solid #41748d; box-shadow: 0 1px 3px rgba(0,0,0,0.05); text-align: center;",
          span(style = "display: block; font-size: 0.72rem; font-weight: 700; color: #64748b; text-transform: uppercase; letter-spacing: 0.5px;", "Sample Size"),
          span(style = "font-size: 1.2rem; font-weight: 800; color: #2c3e50; line-height: 1.8;", paste0("n = ", ks_res$n))
        )
      ),
      
      card(
        style = "border: 1px solid #e2e8f0; border-radius: 8px; box-shadow: none;",
        card_header(class = "bg-light", style = "padding: 10px 14px;", 
                    h6(icon("circle-info"), " Understanding the KS Test Results", style = "margin: 0; font-size: 0.88rem; color: #002855; font-weight: 600;")),
        card_body(
          style = "padding: 12px 14px; font-size: 0.88rem; color: #334155;",
          div(style = "margin-bottom: 10px;",
              strong("What are we testing?"),
              p(style = "margin-top: 4px; margin-bottom: 0;",
                "The Kolmogorov-Smirnov (KS) test measures the maximum distance between the actual empirical cumulative distribution of your dataset and the predictive distribution generated by the fitted Scale-Inflated Pareto model.")
          ),
          hr(style = "margin: 8px 0;"),
          p(style = "margin: 0;",
            strong("Statistical Conclusion: "),
            if (is_good_fit) {
              paste0("Because the p-value (", p_val_fmt, ") is greater than 0.05, there is no statistically significant difference between the observed data and the model. The scale-inflated Pareto model provides a solid fit.")
            } else {
              paste0("Because the p-value (", p_val_fmt, ") is 0.05 or lower, the model distribution differs significantly from the observed data.")
            }
          )
        )
      )
    )
  })
  
  output$ks_ecdf_plot <- renderPlot({
    ks_res <- ks_test_data()
    obs <- ks_res$obs
    sim <- ks_res$sim
    
    f_obs <- ecdf(obs)
    f_sim <- ecdf(sim)
    
    x_grid <- seq(min(c(obs, sim)), max(c(obs, sim)), length.out = 500)
    
    par(mar = c(4, 4.5, 2, 1), mgp = c(2.4, 0.7, 0))
    
    plot(x_grid, f_obs(x_grid), type = "l", col = "#002855", lwd = 2.5,
         xlab = "Observation Value (x)", ylab = "Cumulative Probability F(x)",
         main = "Empirical CDF (Observed Data) vs Model CDF (Predictive)",
         cex.main = 1.05, font.main = 2, frame.plot = TRUE)
    
    lines(x_grid, f_sim(x_grid), col = "#2a9d8f", lwd = 2.5, lty = 2)
    grid(col = "#e2e8f0", lty = "dotted")
    
    legend("bottomright", 
           legend = c("Observed Data (eCDF)", "Model Fit (Predictive CDF)"),
           col = c("#002855", "#2a9d8f"),
           lty = c(1, 2), lwd = 2.5, bty = "o", bg = "white", box.col = "#cbd5e1")
  })
  
  output$comp_estimates_tab <- renderTable({
    res <- mcmc_results()
    data.frame(
      Parameter = c("Alpha (Shape)", "Beta (Scale Inflation)", "Epsilon (Outlier Prob)", "Theta (Threshold)"),
      `Model 1 (With Outliers)` = c(round(res$est_alpha, 4), round(res$est_beta, 4), round(res$est_epsilon, 4), format(round(res$est_theta, 2), big.mark = ",")),
      `Model 2 (Standard Pareto)` = c(round(res$est_alpha_2, 4), "1.0000 (Fixed)", round(res$est_epsilon_2, 4), format(round(res$est_theta_2, 2), big.mark = ",")),
      check.names = FALSE
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  output$alpha_comparison_plot <- renderPlot({
    res <- mcmc_results()
    par(mar = c(4.5, 4.5, 1.5, 1), mgp = c(2.4, 0.7, 0))
    
    if (isTRUE(res$use_rb) && !is.null(res$rb_alpha_density)) {
      plot(res$alpha_grid, res$rb_alpha_density, col = "#002855", type = "l", lwd = 2.5,
           main = "", xlab = "Estimated Alpha", ylab = "Density")
    } else {
      d_a1 <- density(res$alphas)
      plot(d_a1, col = "#002855", lwd = 2.5, main = "", xlab = "Estimated Alpha", ylab = "Density")
    }
    lines(density(res$alphas_2), col = "#e63946", lwd = 2.5, lty = 2)
    legend("topright", legend = c("Model 1 (With Outliers)", "Model 2 (No Outliers)"),
           col = c("#002855", "#e63946"), lty = c(1, 2), lwd = 2, bty = "n")
    mtext(round(res$est_alpha, 2), side = 1, at = res$est_alpha, cex = 0.8, col = "#002855")
  })
  
  output$theta_comparison_plot <- renderPlot({
    res <- mcmc_results()
    par(mar = c(4.5, 4.5, 1.5, 1), mgp = c(2.4, 0.7, 0))
    
    if (isTRUE(res$use_rb) && !is.null(res$rb_theta_density)) {
      plot(res$theta_grid, res$rb_theta_density, col = "#002855", type = "l", lwd = 2.5,
           main = "", xlab = "Estimated Theta", ylab = "Density",
           xlim = range(c(res$thetas, res$thetas_2)))
    } else {
      d1 <- density(res$thetas)
      plot(d1, col = "#002855", lwd = 2.5, main = "", xlab = "Estimated Theta", ylab = "Density",
           xlim = range(c(res$thetas, res$thetas_2)))
    }
    lines(density(res$thetas_2), col = "#e63946", lwd = 2.5, lty = 2)
    legend("topleft", legend = c("Model 1 (With Outliers)", "Model 2 (No Outliers)"),
           col = c("#002855", "#e63946"), lty = c(1, 2), lwd = 2, bty = "n")
    mtext(format(round(res$est_theta, 0), big.mark = ","), side = 1, at = res$est_theta, cex = 0.8, col = "#002855")
  })
  
  output$outlier_probs_plot <- renderPlot({
    res <- mcmc_results()
    post_outlier_probs <- colMeans(res$p_mat)
    ord <- order(res$comp_sim_val)
    
    par(mar = c(4.5, 4.5, 1.5, 1), mgp = c(2.4, 0.7, 0))
    plot(res$comp_sim_val[ord], post_outlier_probs[ord], type = "h", lwd = 2, col = "#002855", 
         ylim = c(0, 1.05),
         xlab = "Observation Value", ylab = "Posterior Outlier Probability P(δ_i = 1)",
         main = "", frame.plot = TRUE)
    points(res$comp_sim_val[ord], post_outlier_probs[ord], pch = 19, col = "#002855", cex = 0.8)
  })
  
  # --------------------------------------------------------------------------
  # SERVER LOGIC - APP 2
  # --------------------------------------------------------------------------
  
  app2_uploaded_df <- reactive({
    req(input$app2_file_upload)
    ext <- tools::file_ext(input$app2_file_upload$name)
    if (ext == "csv") {
      read.csv(input$app2_file_upload$datapath)
    } else if (ext %in% c("xls", "xlsx")) {
      read_excel(input$app2_file_upload$datapath)
    } else {
      showNotification("Invalid format. Upload an Excel (.xlsx/.xls) or CSV file.", type = "error")
      return(NULL)
    }
  })
  
  output$app2_column_select_ui <- renderUI({
    df <- app2_uploaded_df()
    req(df)
    selectInput("app2_selected_col", "Select Value Column (x):", choices = names(df))
  })
  
  app2_analysis_results <- eventReactive(input$app2_run_btn, {
    withProgress(message = 'Estimating Contaminated Pareto Model...', value = 0, {
      
      true_parms <- NULL
      if (input$app2_data_source == "Simulated Pareto Mixture") {
        n       <- input$app2_sim_n
        alpha   <- input$app2_sim_alpha
        theta   <- input$app2_sim_theta
        beta    <- input$app2_sim_beta
        epsilon <- input$app2_sim_epsilon
        
        x <- c(
          extraDistr::rpareto(n = round(n * (1 - epsilon)), a = alpha, b = theta),
          extraDistr::rpareto(n = round(n * epsilon), a = alpha, b = beta * theta)
        )
        true_parms <- c(alpha = alpha, beta = beta, epsilon = epsilon, theta = theta)
      } else if (input$app2_data_source == "Norwegian Oil Dataset") {
        x <- norwegian_oil_data
        true_parms <- NULL
      } else {
        df <- app2_uploaded_df()
        req(df, input$app2_selected_col)
        x <- as.numeric(df[[input$app2_selected_col]])
        x <- x[!is.na(x)]
        true_parms <- NULL
      }
      
      incProgress(0.3, detail = "Running MCEM Gibbs Sampler...")
      
      tst <- MCEM_app2(x, max.tol = input$app2_max_tol, max.iter = input$app2_max_iter)
      
      incProgress(0.7, detail = "Computing MLE Estimates...")
      
      alpha_init   <- 1 / mean(log(x / min(x)))
      beta_init    <- max(x) / min(x)
      epsilon_init <- 0.05
      
      parms_init <- c(log(alpha_init), log(beta_init), qlogis(epsilon_init))
      mle_parms  <- optim(par = parms_init, fn = ll_app2, x = x, control = list(fnscale = -1))
      mle_est    <- c(exp(mle_parms$par[1:2]), plogis(mle_parms$par[3]), min(x))
      names(mle_est) <- c("alpha", "beta", "epsilon", "theta")
      
      incProgress(1.0, detail = "Completed.")
      
      list(
        x = x,
        tst = tst,
        mle_est = mle_est,
        true_parms = true_parms
      )
    })
  })
  
  output$app2_plot_histograms <- renderPlot({
    res <- app2_analysis_results()
    tst <- res$tst
    mle_est <- res$mle_est
    true_parms <- res$true_parms
    
    alphas   <- tst$Gibbs.sample$alphas
    betas    <- tst$Gibbs.sample$betas
    epsilons <- tst$Gibbs.sample$epsilons
    thetas   <- tst$Gibbs.sample$thetas
    
    par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3.5, 3))
    
    plot_param_hist <- function(samples, mle_val, true_val, col_line, param_symbol, round_digits) {
      all_vals <- c(samples, mle_val)
      if (!is.null(true_val)) all_vals <- c(all_vals, true_val)
      x_range <- range(all_vals, na.rm = TRUE)
      
      hist(samples, breaks = 20, probability = TRUE, col = "white", border = "lightgrey",
           main = bquote("Histogram of " ~ .(param_symbol)), xlab = param_symbol, xlim = x_range,
           cex.lab = 1.5, cex.main = 1.5)
      lines(density(samples), col = col_line, lwd = 2)
      
      mcmc_mean <- mean(samples)
      abline(v = mcmc_mean, col = col_line, lty = 2, lwd = 1.8)
      mtext(paste0("MCMC: ", round(mcmc_mean, round_digits)), side = 1, at = mcmc_mean, cex = 1, col = col_line, line = 0.2)
      
      abline(v = mle_val, col = "forestgreen", lty = 4, lwd = 1.8)
      mtext(paste0("MLE: ", round(mle_val, round_digits)), side = 3, at = mle_val, cex = 1, col = "forestgreen", line = 0.2)
      
      if (!is.null(true_val)) {
        abline(v = true_val, col = "orange", lty = 3, lwd = 2)
        mtext(paste0("True: ", round(true_val, round_digits)), side = 3, at = true_val, cex = 1, col = "orange", line = 0.2)
      }
    }
    
    plot_param_hist(alphas, mle_est[1], true_parms["alpha"], "red", expression(alpha), 3)
    plot_param_hist(betas, mle_est[2], true_parms["beta"], "purple", expression(beta), 1)
    plot_param_hist(epsilons, mle_est[3], true_parms["epsilon"], "green", expression(epsilon), 4)
    plot_param_hist(thetas, mle_est[4], true_parms["theta"], "blue", expression(theta), 0)
    
    par(mfrow = c(1, 1))
  })
  
  output$app2_plot_outliers <- renderPlot({
    res <- app2_analysis_results()
    x <- res$x
    prob_outlying <- rowMeans(res$tst$Gibbs.sample$deltas)
    
    df_plot <- data.frame(
      log_x = log(x),
      prob  = prob_outlying,
      Status = factor(
        ifelse(prob_outlying > 0.5, "Classified Outlier (P > 0.5)", "Standard Observations"),
        levels = c("Classified Outlier (P > 0.5)", "Standard Observations")
      )
    )
    
    ggplot(df_plot, aes(x = log_x, y = prob, color = Status)) +
      geom_segment(
        data = subset(df_plot, prob > 0.5),
        aes(x = log_x, xend = log_x, y = 0, yend = prob),
        linewidth = 0.8,
        show.legend = FALSE
      ) +
      geom_point(size = 2.2) +
      scale_color_manual(values = c(
        "Classified Outlier (P > 0.5)" = "#D93B2B",
        "Standard Observations"        = "#3C6382"
      )) +
      scale_y_continuous(limits = c(-0.02, 1.02), breaks = seq(0, 1, 0.2)) +
      labs(
        x = "Log Observations (log(x))",
        y = "Posterior Outlier Probability",
        color = NULL
      ) +
      theme_minimal(base_size = 13) +
      theme(
        legend.position = "top",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = "#EAEAEA", linewidth = 0.5)
      )
  })
  
  output$app2_summary_table <- renderTable({
    res <- app2_analysis_results()
    tst <- res$tst
    
    mcmc_means <- c(
      mean(tst$Gibbs.sample$alphas),
      mean(tst$Gibbs.sample$betas),
      mean(tst$Gibbs.sample$epsilons),
      mean(tst$Gibbs.sample$thetas)
    )
    
    df_summary <- data.frame(
      Parameter = c("Alpha (Shape)", "Beta (Contamination multiplier)", "Epsilon (Outlier Proportion)", "Theta (Scale)"),
      `MCMC Posterior Mean` = mcmc_means,
      `MLE Estimate` = res$mle_est,
      check.names = FALSE
    )
    
    if (!is.null(res$true_parms)) { df_summary$`Ground Truth` <- res$true_parms
    }
    
    df_summary
  }, digits = 4)
  
  output$app2_hyperparameters_table <- renderTable({
    res <- app2_analysis_results()
    hp <- res$tst$Hyper.parameters
    
    data.frame(
      `Target Parameter` = c("Alpha (α)", "Alpha (α)", 
                             "Beta (β)", "Beta (β)", 
                             "Epsilon (ε)", "Epsilon (ε)", 
                             "Theta (θ)", "Theta (θ)"),
      `Symbol` = c("a1", "a2", "b*", "λ", "b1", "b2", "t1", "t2"),
      `Prior Parameter Role` = c("Truncated Gamma Shape", "Truncated Gamma Rate", 
                                 "Shifted Distribution Lower Bound", "Shifted Exponential Rate", 
                                 "Beta Distribution Shape 1", "Beta Distribution Shape 2", 
                                 "Gamma Distribution Shape", "Gamma Distribution Rate"),
      `Value` = c(hp$a1, hp$a2, hp$b_star, hp$lambda, hp$b1, hp$b2, hp$t1, hp$t2),
      check.names = FALSE
    )
  }, digits = 4)
  
  output$app2_data_table <- renderDT({
    res <- app2_analysis_results()
    x <- res$x
    prob_outlying <- rowMeans(res$tst$Gibbs.sample$deltas)
    
    data_df <- data.frame(
      Observation_ID = seq_along(x),
      x = x,
      Log_x = log(x),
      Posterior_Outlier_Prob = round(prob_outlying, 4),
      Is_Outlier = prob_outlying > 0.5
    )
    
    datatable(data_df, options = list(pageLength = 15, autoWidth = TRUE))
  })
}

# ==============================================================================
# RUN APP
# ==============================================================================
shinyApp(ui = ui, server = server)
