

# ---------- load data ----------
source(here("script/import_data.R"))
source(here("script/utils.R"))
load(here("data/vpf.RData"))

ctgmeta <- ctganno_merged
# filter(!is.na(virsorter_category))
smeta <- data.frame(sample_data(pseq)) %>% 
  mutate(samplename=rownames(.))
abundance <- data.frame(otu_table(pseq))

minlen <- min(ctganno_merged$contig_length)
maxlen <- max(ctganno_merged$contig_length)

# ------- Server -------
shinyServer(function(input, output) {
  # ------- render UI ---------
  output$contig_len <- renderUI({
    sliderInput("contig_len", "Contig length",
                min=minlen, 
                max=maxlen, 
                step = 500,
                ticks = FALSE,
                post = " bp",
                value=c(minlen,maxlen))
  })
  
  output$checkv_quality <- renderUI({
    selectInput("checkv_quality", "CheckV quality",
                sort(unique(ctgmeta$checkv_quality)),
                selected = c("Complete", "High-quality"),
                multiple = TRUE)
  })
  
  output$vs2cat <- renderUI({
    selectInput("vs2cat", "VirSorter2 category",
                sort(unique(ctgmeta$virsorter_category)),
                selected = c("1", "2"),
                multiple = TRUE)
  })
  
  output$sample_metadata <- renderUI({
    selectInput("sample_metadata", "Sample metadata",
                colnames(smeta),
                selected = "samplename")
  })
  
  output$ctg_feature_col <- renderUI({
    selectInput("ctg_feature_col", "Feature column",
                colnames(ctgmeta),
                selected = "Family")
  })
  
  output$genomeSel <- renderUI({
    selectInput("genomeSel", "Select viral contig",
                fgbks,
                selected = "DB_S2__NODE_13_length_111139_cov_10.702044-cat_1")
  })
  
  output$geneFeature <- renderUI({
    selectInput("geneFeature", "Select gene annotation",
                colnames(ganno),
                selected = "auxiliary_score")
  })
  
  output$geneLabel <- renderUI({
    selectInput("geneLabel", "Show gene label",
                colnames(ganno),
                selected = "kegg_hit")
  })
  
  # --------- reactive --------------
  
  # -------- filter contigs ---------
  ctgmeta_filtered <- reactive({
    if (is.null(input$checkv_quality)) {
      return(NULL)
    }
    if (input$logic=="AND") {
      ctgmeta %>% 
        filter(contig_length >= input$contig_len[1],
               contig_length <= input$contig_len[2]) %>% 
        filter(virsorter_category %in% input$vs2cat) %>%
        filter(checkv_quality %in% input$checkv_quality)
    } else if (input$logic=="OR") {
      ctgmeta %>% 
        filter(contig_length >= input$contig_len[1],
               contig_length <= input$contig_len[2]) %>% 
        filter(virsorter_category %in% input$vs2cat | 
                 checkv_quality %in% input$checkv_quality)
    }
  })
  
  # --------- load genbank --------------
  gbk <- reactive({
    fgbk <- here(fpath_gbk, paste0(input$genomeSel, ".gbk"))
    gggenomes::read_gbk(fgbk) %>% 
      mutate(gene=name) %>% 
      left_join(ganno, by="gene")
  })
  
  # ------ checkv quality --------
  output$coolplot <- renderPlotly({
    if (is.null(ctgmeta_filtered())) {
      return()
    }
    p <- ctgmeta_filtered() %>% 
      ggplot(aes(x=contig_length, color=checkv_quality, fill=checkv_quality)) +
      geom_histogram(alpha=0.5, position = "identity", bins=30) +
      scale_x_continuous(trans=log10_trans(),
                         breaks = trans_breaks("log10",
                                               function(x) as.integer(10^x))) +
      theme_bw()
    ggplotly(p)
  })
  
  # ------ subset Phyloseq -------
  pseq_subset <- reactive({
    filter_tax_table(pseq, Contig %in% ctgmeta_filtered()$Contig)
    # filter_tax_table(pseq, Contig %in% ctgmeta_shared$key())
  })
  
  # --------- taxonomy bar plot ---------
  output$taxa_abundance <- renderPlotly({
    # if (is.null(ctgmeta_filtered())) {
    #   return()
    # }
    taxa_collapse <-  tax_glom(pseq_subset(), taxrank = input$ctg_feature_col)
    if (input$abundance_scale=="Relative") {
      taxa_collapse <- taxa_collapse %>% 
        transform_sample_counts(function(x) {x/sum(x)})
    }
    p <- taxa_collapse %>% 
      psmelt() %>% 
      filter(Abundance > input$abundance_min_threshold) %>%
      ggplot(aes(x = Sample, y = Abundance, fill=.data[[input$ctg_feature_col]])) +    # Color by Phylum
      geom_bar(stat = "identity", position="stack") +
      guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
      ylab(paste0(input$abundance_scale, " abundance (> ", input$abundance_min_threshold, ")")) +
      theme_bw()
    
    # theme_bw()
    ggplotly(p)
  })
  
  # ------- feature annotation table --------
  output$results <- renderReactable({
    reactable(ctgmeta_filtered(),
              selection = "multiple", 
              filterable = TRUE,
              onClick = "select",
              searchable = TRUE)
  })
  
  # -------- alpha diversity -----------
  output$alpha_diversity <- renderPlot({
    plot_richness(pseq_subset())
    # ggplotly(p)
  })
    
  # ----------- genome plot -----------
  output$genomePlot <- renderPlot({
    gbk() %>% 
      gggenomes() +
      geom_seq() +
      geom_seq_label(size=6) +
      geom_gene(aes(fill=.data[[input$geneFeature]], label=input$geneLabel), position="strand") +
      geom_gene_tag(aes(label=.data[[input$geneLabel]]), size=4) +
      geom_text(aes(label=input$geneLabel, x=-25, y=1.1), size=6) +
      theme(legend.position = "bottom",
            legend.text=element_text(size=16),
            text=element_text(size=16),
            axis.text.x = element_text(size=16))
  })
  
  # ----------- gene annotation ------------
  output$geneAnno <- renderReactable({
    reactable(ganno,
              selection = "multiple", 
              filterable = TRUE,
              onClick = "select")
  })
  
  # --------- abundance ---------
  output$count <- renderReactable({
    reactable(df_count, filterable = TRUE, onClick = "select", resizable = TRUE)
  })
  
  output$tpm <- renderReactable({
    reactable(df_tpm, filterable = TRUE, onClick = "select", resizable = TRUE)
  })

})
