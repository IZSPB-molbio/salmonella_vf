#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# https://community.rstudio.com/t/failing-to-deploy-shinyapp-depending-on-bioconductor-packages/6970/3
library(BiocManager)
options(repos = BiocManager::repositories())
library(heatmaply)
library(shiny)
library(shinydashboard)
library(DT)
library(tidyverse)
library(plotly)
library(cetcolor)
library(ggtree)
library(ape)
library(here)
getwd()
read.abricate <- function(x){
    t <- data.frame(read_delim(x, delim = "\t", col_names = TRUE, col_types = "ffnnfffffnnffff"))
    colnames(t) <- c("isolate","sequence", "start", "end",
                     "strand", "gene", "coverage", "coverage_map", "gaps", "percent_coverage",
                     "percent_identity", "db", "db_accession", "product", "resistance")
    t$resistance <- gsub(";", "; ", t$resistance)
    t$product <- gsub(":", ": ", t$product)
    return(t)
}

get.heatmap.gene.list <- function(my_isolate, my_db, abricate_results){
    gene.list <- abricate_results %>%
        subset(isolate == my_isolate) %>%
        subset(db == my_db) %>%
        select(gene, product, resistance)
    gene.list
}

# get.heatmap.compound.gene.list(my_isolate = isolate,
#   my_db = db,
#   abricate_results = abricate.results)
get.heatmap.compound.gene.list <- function(my_isolate, my_comp, abricate_results){
    aaa <- abricate_results %>%
        mutate(across(where(is.factor), as.character))
    p <- aaa %>%
        dplyr::filter(grepl(my_comp, .$resistance)) %>%
        dplyr::filter(isolate == my_isolate) %>%
        inner_join(aaa, by = c("sequence", "isolate"), suffix = c(".x", ""), keep = FALSE) %>%
        select(isolate, sequence, gene, db, product, resistance)
    p
}

get.barchart.gene.list <- function(my_isolate, my_comp, abricate_results){
    aaa <- abricate_results %>%
        mutate(across(where(is.factor), as.character))
    p <- aaa %>%
        dplyr::filter(grepl(my_comp, .$resistance)) %>%
        dplyr::filter(isolate == my_isolate) %>%
        inner_join(aaa, by = c("sequence", "isolate"), suffix = c(".x", ""), keep = FALSE) %>%
        select(isolate, sequence, db, product, resistance)
    # subset(grepl(my_comp, .$resistance)) %>%
    # subset(isolate == my_isolate) #%>%
    #     inner_join(aaa, by = "sequence") %>%
    #     select(isolate, sequence, db, product, resistance)
    p
}

# craft the table
abricate.results.files <- list.files(path = "data/annotation/abricate", pattern = "*_*.out", full.names = TRUE)
# abricate.results.files <- list.files(path = here("data", "annotation", "abricate"), pattern = "*_*.out", full.names = TRUE)
abricate.results <- data.frame(do.call("rbind", lapply(abricate.results.files, read.abricate))) %>% mutate(resistance=tolower(resistance))

# heatmap with number of genes/sample/db
## build matrix of occurrences
abricate.results.counts <- abricate.results %>%
    group_by(isolate, db) %>%
    summarise(n_hits = n_distinct(sequence)) %>%
    ungroup() %>% #subset(isolate == "100658")
    spread(isolate, n_hits, fill = 0)

abricate.results.counts.rownames <- abricate.results.counts$db
abricate.results.counts <- select(abricate.results.counts, -db)
rownames(abricate.results.counts) <- abricate.results.counts.rownames
abricate.results.counts <- as.matrix(abricate.results.counts) %>% t()

heatmap_counts <- plot_ly(x = colnames(abricate.results.counts),
                          y = as.character(rownames(abricate.results.counts)),
                          z = abricate.results.counts,
                          type="heatmap") %>% 
    layout(yaxis = list(type = "category"))

# number of occurrences of resistance per sample
# https://stackoverflow.com/a/38108559
resistance <- strsplit(as.character(abricate.results$resistance), "; ")
l_resistance <- lengths(resistance)
max_l_resistance <- max(l_resistance)
resistance <- t(sapply(resistance[as.logical(l_resistance)], function(a) c(a, rep("",max_l_resistance-length(a)))))

compound.counts <- abricate.results %>%
    select(isolate, db) %>%
    cbind(resistance) %>%
    subset(!is.na(`1`)) %>% #head(20) %>%
    pivot_longer(3:ncol(.), values_to = "compound") %>%
    group_by(isolate, compound, db) %>% count(compound) %>% ungroup() %>% subset(compound != "") #%>% dplyr::mutate(across(compound), factor)

# phylogenetic tree
phy <- ape::read.tree('data/pangenome/roary_1617111546/core_gene_alignment.tree.rooted.Svictoria.newick')
# getwd()
# the actual shinyapp starts here
ui <- dashboardPage(
    dashboardHeader(title = "Salmonella virulence factors: dashboard",
                    # tags$li(img(src = "data/izs_logo.jpeg"),
                    #         class = "dropdown")
                    tags$li(a(href = 'http://www.izsfg.it',
                              img(src = "data/izs_logo.jpeg",
                                  title = "IZSPB", height = "60px"),
                              style = "padding-top:10px; padding-bottom:10px;"),
                            class = "dropdown")
                    ),
    dashboardSidebar(
        sidebarMenu(
            menuItem("Introduction", tabName = "intro", icon = icon("lightbulb")),
            menuItem("ABRicate", icon = icon("binoculars"),
                     menuSubItem("ABRicate: overview", tabName = "overview", icon = icon("binoculars")),
                     menuSubItem("ABRicate: results", icon = icon("th"), tabName = "results")),
            menuItem("Phylogeny", icon = icon("tree"), tabName = "phylogeny")
            # menuItem("ABRicate: overview", tabName = "overview", icon = icon("binoculars")),
            # menuItem("ABRicate: results", icon = icon("th"), tabName = "results"),
            # menuItem("Phylogeny", icon = icon("tree"), tabName = "phylogeny")
        )
    ),
    dashboardBody(
        # Boxes need to be put in a row (or column)
        tabItems(
            tabItem(tabName = "intro", h2("Introduction"),
                    fluidRow(
                        column(12,
                               includeMarkdown("data/introduction.md"))
                    )),
            tabItem(tabName = "overview", h2("Overview of results"),
                    fluidRow(
                        box(# wanna align? https://stackoverflow.com/questions/29738975/how-to-align-a-group-of-checkboxgroupinput-in-r-shiny
                            checkboxGroupInput("show_dbs", "DBs to show results for:",
                                               unique(abricate.results$db),
                                               inline = TRUE,
                                               selected = unique(abricate.results$db))),
                        box(
                            radioButtons(inputId = "show_resistance_dbs",
                                         label = "DB to show resistance for:",
                                         inline = TRUE,
                                         choices = unique(compound.counts$db),
                                         selected = unique(compound.counts$db)[1])
                        )
                    ),
                    fluidRow(
                        box(title = "Number of hits per sample (click on a cell to show details about the hits in the table below)",
                            plotlyOutput("heatmap"), height = 860),
                        box(title = "Number of resistance hits per sample (click on a cell to display all hits for a cds)",
                            plotlyOutput("heatmap.compounds"), height = 860),
                        fluidRow(
                            box(DT::dataTableOutput("heatmap.genes")),
                            box(DT::dataTableOutput("heatmap.compound.genes")),
                        ),
                    )),
            tabItem(tabName = "results", h2("Results"),
                    fluidRow(
                        column(7,
                               includeMarkdown("data/description.md")),
                        column(4,
                               checkboxGroupInput("show_vars", "Columns to show:",
                                                  names(abricate.results),
                                                  selected = c("isolate", "gene", "percent_identity",
                                                               "db", "resistance")))),
                    DT::dataTableOutput("salmonella.abricate.table"), height = 800),
            tabItem(tabName = "phylogeny", h2("Phylogeny"),
                    # mainPanel(
                    #     plotOutput("phy"
                    #     )
                    fluidRow(box(plotOutput("phy"), width = 12, height = "100%"))               
                    )
        ),
    )
)

server <- function(input, output) {
    set.seed(122)
    output$heatmap <- renderPlotly({
        plot1 <- abricate.results.counts[, input$show_dbs, drop=FALSE] %>%
            plot_ly(x = colnames(.),
                    y = as.character(rownames(.)),
                    z = .,
                    xgap = 0.4, ygap = 0.4, 
                    hovertemplate = paste(
                        "isolate: %{y}\ndb: %{x}\nnumber of hits: %{z}<extra></extra>"),
                    source = "heatmap",
                    type="heatmap", height = 800) %>%
            layout(yaxis = list(title = "isolate", type = "category"),
                   xaxis = list(title = "abricate db"))
        #plot1 %>% event_register("plotly_click")
    })
    output$heatmap.compounds <- renderPlotly({
        compound.counts.matrix <- compound.counts %>%
            # subset(db == "card", select = -db) %>%
            subset(db == input$show_resistance_dbs, select = -db) %>%
            # group_by(isolate, db) %>%
            # summarise(n_hits = n_distinct(sequence)) %>%
            # ungroup() %>% #subset(isolate == "100658")
            spread(compound, n, fill = 0)
        compound.counts.matrix.rownames <- compound.counts.matrix$isolate
        #mamt <- dplyr::select(compound.counts.matrix, -c(isolate))
        compound.counts.matrix <- as.matrix(select(compound.counts.matrix, -isolate))
        rownames(compound.counts.matrix) <- compound.counts.matrix.rownames
        ### Heatmaply plot (but there are issues with exposing source)
        plot3_custom_hovertext <- compound.counts.matrix %>% as.data.frame()
        plot3_custom_hovertext[] <- paste("isolate: ", rownames(plot3_custom_hovertext))
        plot3_custom_hovertext[] <- lapply(colnames(plot3_custom_hovertext), function(colname) {
            paste0(plot3_custom_hovertext[, colname], ", ", colname)
        })
        # plot3_custom_hovertext[] <- paste("isolate:", rownames(compound.counts.matrix),
        #                                   "\ncompound:", colnames(compound.counts.matrix),
        #                                   "\nn of hits:", compound.counts.matrix)
        plot3 <- compound.counts.matrix %>%
            heatmaply(., source = "heatmap.compounds.2",
                      xgap = 0.9, ygap = 0.9,
                      colors = cet_pal(2, name = "r1"),
                      plot_method = "plotly",
                      # source = "heatmap.compounds",
                      height = 800,
                      custom_hovertext = plot3_custom_hovertext
            ) %>%
            layout(yaxis = list(title = "isolate", type = "category"),
                   xaxis = list(title = "compound", tickangle = -45))
        ### End of heatmaply data
        #
        ### plotly heatmap
        # # add dendrogram https://stackoverflow.com/questions/43794870/plotting-a-clustered-heatmap-with-dendrograms-using-rs-plotly
        # #
        # # dendrogram data
        # dd.col <- as.dendrogram(hclust(dist(compound.counts.matrix)))
        # dd.row <- as.dendrogram(hclust(dist(t(compound.counts.matrix))))
        # dx <- dendro_data(dd.row)
        # dy <- dendro_data(dd.col)
        # # helper function for creating dendograms
        # ggdend <- function(df) {
        #   ggplot() +
        #     geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend)) +
        #     labs(x = "", y = "") + theme_minimal() +
        #     theme(axis.text = element_blank(), axis.ticks = element_blank(),
        #           panel.grid = element_blank())
        # }
        # # x/y dendograms
        # px <- ggdend(dx$segments)
        # py <- ggdend(dy$segments) + coord_flip()
        
        plot2 <- compound.counts.matrix %>%
            plot_ly(x = colnames(.),
                    y = as.character(rownames(.)),
                    z = .,
                    xgap = 0.4, ygap = 0.4,
                    colors = colorRamp(c("red", "green")),
                    hovertemplate = paste(
                        "isolate: %{y}\ncompound: %{x}\nnumber of hits: %{z}<extra></extra>"),
                    source = "heatmap.compounds",
                    type="heatmap", height = 800) %>%
            layout(yaxis = list(title = "isolate", type = "category"),
                   xaxis = list(title = "compound", tickangle = -45))
        #subplot(px, plot2, )
        #plot1 %>% event_register("plotly_click")
    })
    # https://webinars.cpsievert.me/20180220/#10
    # try to get the list of genes hits
    # get.heatmap.gene.list <- function(my_isolate, my_db, abricate_results)
    output$heatmap.genes <- DT::renderDataTable({
        isolate <- event_data("plotly_click", source = "heatmap")$y
        db <- event_data("plotly_click", source = "heatmap")$x
        #isolate
        get.heatmap.gene.list(my_isolate = isolate,
                              my_db = db,
                              abricate_results = abricate.results)
        # "List of genes!"
    }, rownames = FALSE)
    output$heatmap.compound.genes <- DT::renderDataTable({
        if (is.null(event_data("plotly_click", source = "heatmap.compounds"))) {data.frame(Message=c("Click on a point in the heatmap!"))} else {
            isolate <- event_data("plotly_click", source = "heatmap.compounds")$y
            compound <- event_data("plotly_click", source = "heatmap.compounds")$x
            p <- get.heatmap.compound.gene.list(my_isolate = isolate,
                                                my_comp = compound,
                                                abricate_results = abricate.results)
            if(nrow(p) != 0){
                # print(unique(p$sequence))
                # print(length(unique(p$sequence)))
                datatable(p) %>% formatStyle(
                    'sequence',
                    target = 'row',
                    backgroundColor = styleEqual(c(unique(p$sequence)),
                                                 paste0(
                                                     colorRampPalette(RColorBrewer::brewer.pal(8, 'Set2'))(length(unique(p$sequence))), "33"
                                                 ))
                )} else {data.frame(Message=c("No data to be shown!"))} }
        # "List of genes!"
    }, rownames = FALSE)
    output$barchart.genes <- DT::renderDataTable({
        my_isolate <- event_data("plotly_click", source = "barchart")$y
        comps <- compound.counts %>%
            # subset(db %in% c("card")) %>%
            as_tibble() %>%
            subset(db == input$show_resistance_dbs) %>%
            mutate(across(where(is.factor), as.character))
        print(head(comps, 60))
        comp <- (comps %>%
                     dplyr::pull(compound))[as.integer(event_data("plotly_click", source = "barchart")$curveNumber)+1]# %>%
        # unique())[as.integer(event_data("plotly_click", source = "barchart")$curveNumber)+1]
        print(comps %>%
                  dplyr::pull(compound))
        print(event_data("plotly_click", source = "barchart"))
        print(paste("isolate is ", my_isolate))
        print(paste("index is ", event_data("plotly_click", source = "barchart")$curveNumber+1))
        print(paste("comp is ", comp, " which is ", typeof(comp)))
        # if(comp == ""){comp <- "qwerty"}
        print(paste("comp is ", comp, " which is ", typeof(comp)))
        get.barchart.gene.list(my_isolate = my_isolate, my_comp = comp, abricate_results = abricate.results)
    }, rownames = FALSE)
    output$barchart <- renderPlotly(
        plot2 <- compound.counts %>%
            subset(db %in% input$show_resistance_dbs) %>%
            plot_ly(., x = ~n, y = ~isolate, hovertemplate = paste('compound:', 
                                                                   .$compound, 
                                                                   '\nisolate: %{y}\nnumber of hits: %{x}\ndb:',
                                                                   .$db,
                                                                   '<extra></extra>'),
                    source = "barchart",
                    orientation = 'h', type = 'bar', color = ~compound, height = 800) %>%
            layout(yaxis = list(title = 'isolate'),
                   xaxis = list(title = 'number of hits'),
                   barmode = 'stack')
    )
    output$salmonella.abricate.table = DT::renderDataTable({
        datatable(abricate.results[, input$show_vars, drop=FALSE],
                  # abricate.results[, input$show_vars, drop=FALSE],
                  filter = list(position="top",
                                clear=TRUE,
                                plain=FALSE),
                  extensions = c("Buttons",
                                 "ColReorder",
                                 "FixedHeader"),
                  options = list(
                      autoWidth = TRUE,
                      dom = "Blfrtip",
                      colReorder = TRUE,
                      fixedHeader = TRUE,
                      pageLength = 50,
                      lengthMenu = c(5, 10, 25, 50, 100, 200, 500, 1000),
                      buttons = list(
                          # list(extend = "colvis", columns = 1:ncol(.)),
                          c('copy', 'csv', 'excel'))
                  ),
                  rownames = FALSE
        )
    })
    output$phy <- renderPlot(ggtree(phy) + geom_tiplab(size=3), width = "auto", height = "auto")
}

shinyApp(ui, server)

### old version

# ui <- fluidPage(
#     
#     # Application title
#     #titlePanel("TE Biogas: bin annotations"),
#     # tags$img(src = "SLUBI_logo_smaller.png", align="right"),
#     
#     fluidRow(
#         column(7, 
#                includeMarkdown("data/description.md")),
#         column(4,
#                checkboxGroupInput("show_vars", "Columns to show:",
#                                   names(abricate.results), selected = c("isolate", "gene", "percent_identity", "db", "resistance"))),
#         # column(4,
#         #        checkboxGroupInput("show_vars", "Databases to show results for:",
#         #                           unique(abricate.results$db), selected = unique(abricate.results$db)),
#         # )
#         
#         # box(title = "",
#         #     status = "primary",
#         #     p = "BWE",
#         #     solidHeader = F,
#         #     collapsible = F,
#         #     width = 12)
#         # fluidRow(column(width = 10, box(p=instructions)),
#         #          column(width = 2, align = "center",
#         #                 img(src="SLUBI_logo_smaller.png", width=100))))
#         # ),
#         # Sidebar with checkboxes for displaying columns 
#         # sidebarLayout(
#         #     sidebarPanel(
#         #         checkboxGroupInput("show_vars", "Columns/annotations to show:",
#         #                            names(gene.go_ko), selected = names(gene.go_ko))
#         # sliderInput("bins",
#         #             "Number of bins:",
#         #             min = 1,
#         #             max = 50,
#         #             value = 30)
#     ),
#     
#     DT::dataTableOutput("salmonella.abricate.table")
# )
# 
# # # Define UI for application that draws a histogram
# # ui <- fluidPage(
# # 
# #     # Application title
# #     titlePanel("Old Faithful Geyser Data"),
# # 
# #     # Sidebar with a slider input for number of bins 
# #     sidebarLayout(
# #         sidebarPanel(
# #             sliderInput("bins",
# #                         "Number of bins:",
# #                         min = 1,
# #                         max = 50,
# #                         value = 30)
# #         ),
# # 
# #         # Show a plot of the generated distribution
# #         mainPanel(
# #            plotOutput("distPlot")
# #         )
# #     )
# # )
# 
# server <- function(input, output) {
#     
#     output$salmonella.abricate.table = DT::renderDataTable({
#         DT::datatable(abricate.results[, input$show_vars, drop=FALSE],
#                       filter = list(position="top",
#                                     clear=TRUE,
#                                     plain=FALSE),
#                       extensions = c("Buttons",
#                                      "ColReorder",
#                                      "FixedHeader"),
#                       options = list(
#                           autoWidth = TRUE,
#                           dom = "Blfrtip",
#                           colReorder = TRUE,
#                           fixedHeader = TRUE,
#                           pageLength = 50,
#                           lengthMenu = c(5, 10, 25, 50, 100, 200, 500, 1000),
#                           buttons = list(
#                               # list(extend = "colvis", columns = 1:ncol(.)),
#                               c('copy', 'csv', 'excel'))
#                       )
#         )
#     })
# }
# 
# # Run the application 
# shinyApp(ui = ui, server = server)
