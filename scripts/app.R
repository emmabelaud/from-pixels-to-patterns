library(shiny)
library(tidyverse)
library(plotly)

# ══════════════════════════════════════════════════════════════════════════════
# DATA — adjust paths to match your project structure
# ══════════════════════════════════════════════════════════════════════════════
microclimate_index <- readLines("../output/microclimate_index.txt")
mc1_pred <- microclimate_index[1]
mc2_pred <- microclimate_index[2]
mc1_resp <- "microclimate_1"
mc2_resp <- "microclimate_2"

data_sem_results <- read.csv("../output/SEM_results_database.csv",
                             na.strings = c("", "NA")) %>%
  mutate(
    date_start = as.Date(date_start),
    taxon      = if_else(taxon == "all_taxa", "Community", taxon),
    predictor  = case_when(
      predictor == "microclimate_1" ~ microclimate_index[1],
      predictor == "microclimate_2" ~ microclimate_index[2],
      TRUE                          ~ predictor
    )
  )

indirect_long <- read.csv("../output/indirect_effects.csv",
                          na.strings = c("", "NA")) %>%
  mutate(date_start = as.Date(date_start),
         taxon      = if_else(taxon == "all_taxa", "Community", taxon)) %>%
  pivot_longer(
    cols      = c(via_mc1, via_mc2, via_root,
                  via_mc1_root, via_mc2_root,
                  total_indirect, total_effect),
    names_to  = "path_id",
    values_to = "std_estimate"
  ) %>%
  filter(!is.na(std_estimate)) %>%
  mutate(predictor = NA_character_, response = NA_character_)

# ── Path metadata ─────────────────────────────────────────────────────────────
path_meta <- bind_rows(
  path_meta,
  tribble(
    ~path_id,         ~beta,  ~label,             ~colour,   ~dash,
    "via_mc1",        "β₁×β₆","lu→mc1→fauna",     "#c2cab2", "dot",
    "via_mc2",        "β₂×β₇","lu→mc2→fauna",     "#555d45", "dot",
    "via_root",       "β₅×β₈","lu→root→fauna",    "#955f42", "dot",
    "via_mc1_root",   "β₁×β₃×β₈","lu→mc1→root→fauna","#dea88b","dot",
    "via_mc2_root",   "β₂×β₄×β₈","lu→mc2→root→fauna","#d1855c","dot",
    "total_indirect", "Σ indirect","Total indirect","#999999", "solid",
    "total_effect",   "Total","Total effect",      "#444444", "solid"
  )
)

# ── Prepare smoothed data ─────────────────────────────────────────────────────
# Direct SEM paths
direct_paths <- data_sem_results %>%
  filter(!str_detect(predictor, "~~"), !is.na(std_estimate)) %>%
  mutate(path_id = paste0(predictor, "__", response)) %>%
  filter(path_id %in% path_meta$path_id)

# Computed indirect paths
indirect_paths <- indirect_long %>%
  filter(path_id %in% path_meta$path_id)

plot_data <- bind_rows(direct_paths, indirect_paths) %>%
  left_join(path_meta, by = "path_id") %>%
  group_by(path_id) %>%
  mutate(
    lo           = quantile(std_estimate, 0.025, na.rm = TRUE),
    hi           = quantile(std_estimate, 0.975, na.rm = TRUE),
    std_estimate = pmax(lo, pmin(hi, std_estimate))
  ) %>%
  ungroup()

taxa_list <- c("Community",
               sort(setdiff(unique(plot_data$taxon), "Community")))

smooth_series <- function(df, span = 0.3) {
  df <- arrange(df, date_start)
  if (nrow(df) < 4) return(mutate(df, smoothed = std_estimate))
  fit <- tryCatch(
    loess(std_estimate ~ as.numeric(date_start), data = df, span = span),
    error = function(e) NULL
  )
  df$smoothed <- if (!is.null(fit)) predict(fit) else df$std_estimate
  df
}

smoothed_data <- plot_data %>%
  group_by(path_id, taxon) %>%
  group_modify(~ smooth_series(.x)) %>%
  ungroup()

# ══════════════════════════════════════════════════════════════════════════════
# UI
# ══════════════════════════════════════════════════════════════════════════════
ui <- fluidPage(
  theme = bslib::bs_theme(bootswatch = "flatly", base_font_size = "13px"),
  tags$head(tags$style(HTML("
    .well { padding: 10px 12px; background: #fafafa;
            border: 0.5px solid #ddd; border-radius: 8px; }
    .cb-group label { font-size: 12px; }
    .section-title {
      font-size: 10px; font-weight: 600; color: #666;
      text-transform: uppercase; letter-spacing: .05em;
      border-bottom: 0.5px solid #eee; padding-bottom: 4px;
      margin-bottom: 6px; margin-top: 10px;
    }
    .section-title:first-child { margin-top: 0; }
    #plot { border-radius: 8px; }
  "))),
  
  fluidRow(
    # ── Sidebar ──────────────────────────────────────────────────────────────
    column(2,
           wellPanel(
             div(class = "section-title", "Paths"),
             actionLink("paths_all",  "select all",  style = "font-size:11px"),
             " · ",
             actionLink("paths_none", "clear all", style = "font-size:11px"),
             tags$br(), tags$br(),
             # One checkbox per path, pre-coloured via inline HTML
             div(class = "cb-group",
                 checkboxGroupInput(
                   "sel_paths", label = NULL,
                   choiceNames  = purrr::map(seq_len(nrow(path_meta)), function(i) {
                     pm <- path_meta[i, ]
                     tags$span(
                       tags$span(style = paste0(
                         "display:inline-block;width:10px;height:10px;",
                         "border-radius:2px;background:", pm$colour,
                         ";margin-right:5px;vertical-align:middle;"
                       )),
                       paste0(pm$beta, ": ", pm$label)
                     )
                   }),
                   choiceValues = path_meta$path_id,
                   selected     = path_meta$path_id
                 )
             ),
             
             div(class = "section-title", "Taxa"),
             actionLink("taxa_all",  "select all",  style = "font-size:11px"),
             " · ",
             actionLink("taxa_none", "clear all", style = "font-size:11px"),
             tags$br(), tags$br(),
             div(class = "cb-group",
                 checkboxGroupInput(
                   "sel_taxa", label = NULL,
                   choiceNames  = purrr::map(taxa_list, function(t) {
                     col <- if (t == "Community") "#3d5a3e" else "#aaaaaa"
                     tags$span(
                       tags$span(style = paste0(
                         "display:inline-block;width:8px;height:8px;",
                         "border-radius:50%;background:", col,
                         ";margin-right:5px;vertical-align:middle;"
                       )),
                       t
                     )
                   }),
                   choiceValues = taxa_list,
                   selected     = c("Community", taxa_list[2])
                 )
             ),
             
             tags$hr(style = "margin: 10px 0;"),
             sliderInput("span", "Smoothing span",
                         min = 0.1, max = 0.8, value = 0.3, step = 0.05),
             checkboxInput("show_zero", "Show zero line", value = TRUE)
           )
    ),
    
    # ── Plot ─────────────────────────────────────────────────────────────────
    column(10,
           plotlyOutput("plot", height = "580px")
    )
  )
)

# ══════════════════════════════════════════════════════════════════════════════
# SERVER
# ══════════════════════════════════════════════════════════════════════════════
server <- function(input, output, session) {
  
  # Toggle all / clear all — paths
  observeEvent(input$paths_all,  {
    updateCheckboxGroupInput(session, "sel_paths", selected = path_meta$path_id)
  })
  observeEvent(input$paths_none, {
    updateCheckboxGroupInput(session, "sel_paths", selected = character(0))
  })
  
  # Toggle all / clear all — taxa
  observeEvent(input$taxa_all,  {
    updateCheckboxGroupInput(session, "sel_taxa", selected = taxa_list)
  })
  observeEvent(input$taxa_none, {
    updateCheckboxGroupInput(session, "sel_taxa", selected = character(0))
  })
  
  # Re-smooth when span changes
  smoothed <- reactive({
    plot_data %>%
      group_by(path_id, taxon) %>%
      group_modify(~ smooth_series(.x, span = input$span)) %>%
      ungroup()
  })
  
  output$plot <- renderPlotly({
    req(length(input$sel_paths) > 0, length(input$sel_taxa) > 0)
    
    df <- smoothed() %>%
      filter(path_id %in% input$sel_paths,
             taxon   %in% input$sel_taxa)
    
    if (nrow(df) == 0) {
      return(plotly_empty(type = "scatter") %>%
               layout(title = "No data for current selection"))
    }
    
    fig <- plot_ly()
    
    for (tx in intersect(taxa_list, input$sel_taxa)) {
      is_community <- tx == "Community"
      lw      <- if (is_community) 3 else 2.5
      opacity <- if (is_community) 1.0 else 0.8
      
      for (pi in seq_len(nrow(path_meta))) {
        pm <- path_meta[pi, ]
        if (!pm$path_id %in% input$sel_paths) next
        
        sub <- df %>% filter(path_id == pm$path_id, taxon == tx)
        if (nrow(sub) < 2) next
        
        fig <- add_trace(fig,
                         data          = sub,
                         x             = ~date_start,
                         y             = ~smoothed,
                         type          = "scatter",
                         mode          = "lines",
                         name          = paste0(pm$beta, ": ", pm$label),
                         legendgroup   = pm$path_id,
                         showlegend    = (tx == input$sel_taxa[1] ||
                                            tx == "Community" && "Community" %in% input$sel_taxa),
                         line          = list(
                           color = pm$colour,
                           width = lw,
                           dash  = if (pm$dash == "dash") "dash" else "solid"
                         ),
                         opacity       = opacity,
                         hovertemplate = paste0(
                           "<b>", tx, "</b><br>",
                           pm$beta, ": ", pm$label, "<br>",
                           "β = %{y:.3f}<extra></extra>"
                         )
        )
      }
    }
    
    fig %>%
      layout(
        xaxis = list(
          title      = "",
          tickformat = "%b %Y",
          gridcolor  = "rgba(200,200,200,0.3)"
        ),
        yaxis = list(
          title         = "Standardised β",
          gridcolor     = "rgba(200,200,200,0.3)",
          zeroline      = input$show_zero,
          zerolinecolor = "rgba(80,80,80,0.5)",
          zerolinewidth = 1
        ),
        legend = list(
          title   = list(text = "<b>Path</b>"),
          x = 1.01, y = 1,
          bgcolor = "rgba(255,255,255,0.9)",
          bordercolor = "rgba(180,180,180,0.4)",
          borderwidth = 0.5
        ),
        hovermode    = "x unified",
        plot_bgcolor = "white",
        paper_bgcolor= "white",
        margin       = list(t = 20, r = 20, l = 60, b = 50)
      ) %>%
      config(
        displaylogo = FALSE,
        modeBarButtonsToRemove = c("zoom2d","pan2d","select2d",
                                   "lasso2d","autoScale2d")
      )
  })
}

shinyApp(ui, server)
