library(shiny)
# library(shinythemes)
library(dplyr)
library(readr)
library(ggplot2)
# remotes::install_github('oganm/allenBrain')
library(allenBrain)
library(magick)
library(plotly)
# library(EnhancedVolcano)
# library(ggplotify)
library(BiocManager)
# library(Giotto)
library(Seurat)

options(repos = BiocManager::repositories())

# getwd()

#Reads in data from directory...
spatial.data.mast <- Read10X(data.dir = getwd())

#Reads in csv that labels each cell class
cell_class <- read.csv("cell_class.csv")

#rownames needs to be exact same as seurat for this to work
cell_class <- as.data.frame(cell_class[,-1], row.names = cell_class[,1])
colnames(cell_class) <- c('cell.class')


#create searat object for analysis... seems to take long time...
scrna.mast <- CreateSeuratObject(counts = spatial.data.mast, meta.data = cell_class)


# standard workflow
scrna.mast <- NormalizeData(object = scrna.mast)
scrna.mast <- FindVariableFeatures(object = scrna.mast)
scrna.mast <- ScaleData(object = scrna.mast)
scrna.mast <- RunPCA(object = scrna.mast)
scrna.mast <- FindNeighbors(object = scrna.mast)
scrna.mast <- FindClusters(object = scrna.mast)
scrna.mast <- RunTSNE(object = scrna.mast)



# get a list of structure names and ids
IDs = getStructureIDs()
IDs %>% head

brain_regions_names_comb <- c("Olfactory tubercle - OT" = "Olfactory tubercle", 
                              "Bed nuclei of the stria terminalis - BST" = "Bed nuclei of the stria terminalis",
                              "Bed nuclei of the stria terminalis, anterior division, anterolateral area - BSTal" = "Bed nuclei of the stria terminalis, anterior division, anterolateral area",
                              "Bed nuclei of the stria terminalis, anterior division, ventral nucleus - BSTv" = "Bed nuclei of the stria terminalis, anterior division, ventral nucleus",
                              "Bed nuclei of the stria terminalis, posterior division - BSTp" = "Bed nuclei of the stria terminalis, posterior division",
                              "Bed nucleus of the anterior commissure - BAC" = "Bed nucleus of the anterior commissure",
                              "Paraventricular nucleus of the thalamus - PVT" = "Paraventricular nucleus of the thalamus",
                              "Paraventricular hypothalamic nucleus - PVH" = "Paraventricular hypothalamic nucleus",
                              "Periventricular hypothalamic nucleus, anterior part - PVa" = "Periventricular hypothalamic nucleus, anterior part",
                              "Anteroventral periventricular nucleus - AVPV" = "Anteroventral periventricular nucleus",
                              "Median preoptic nucleus - MEPO" = "Median preoptic nucleus",
                              "Medial preoptic area - MPO" = "Medial preoptic area",
                              "Parastrial nucleus - PS" = "Parastrial nucleus",
                              "Periventricular hypothalamic nucleus, preoptic part - PVpo" = "Periventricular hypothalamic nucleus, preoptic part",
                              "Suprachiasmatic nucleus - SCH" = "Suprachiasmatic nucleus",
                              "Ventromedial preoptic nucleus - VMPO" = "Ventromedial preoptic nucleus",
                              "Ventrolateral preoptic nucleus - VLPO" = "Ventrolateral preoptic nucleus",
                              "Medial preoptic nucleus - MPN" = "Medial preoptic nucleus",
                              "Lateral preoptic area - LPO" = "Lateral preoptic area",
                              "hypothalamus related - mfsbshy" = "hypothalamus related")


oprd1.markers <- read_csv("oprd1.markers.csv")
oprm1.markers <- read_csv("oprm1.markers.csv")
oprk1.markers <- read_csv("oprk1.markers.csv")






ui = fluidPage(
  # Application title
  titlePanel("Gene Markers for Opioid Receptor Cell Types in Hypothalamus Preoptic"),
  tabsetPanel(tabPanel('Delta Opioid Receptor DGE Results', splitLayout(imageOutput('delta_diff_photo')),
                                                 dataTableOutput('delta_diff_table')),
              tabPanel('Kappa Opioid Receptor DGE Results', splitLayout(imageOutput('kappa_diff_photo')),
                       dataTableOutput('kappa_diff_table')), 
              tabPanel('Mu Opioid Receptor DGE Results', splitLayout(imageOutput('mu_diff_photo')),
                       dataTableOutput('mu_diff_table')),
              tabPanel('Analytical Visualization', 
                       tabsetPanel(tabPanel("Cell Spatial Mapping and Visualization", 
                                            sidebarLayout(sidebarPanel(textInput(inputId = "datasetID", 
                                                                                 "Enter Gene ID", 
                                                                                 placeholder = "Gene ID:", 
                                                                                 value = 'Oprk1'),
                                                                       selectInput(inputId = "granuleID", 
                                                                                   "Brain Region:", 
                                                                                   choices = brain_regions_names_comb, 
                                                                                   selected = "Medial preoptic area")), 
                                                          mainPanel(plotlyOutput('img_1'),
                                                                    plotlyOutput('img_2'),
                                                                    plotlyOutput('img_3'),
                                                                    plotlyOutput('img_4')))), 
                                    tabPanel("Cell Spatial Mapping and Visualization", sidebarLayout(sidebarPanel(textInput(inputId = "gene1", 
                                                                                                        "Enter Gene 1", 
                                                                                                        placeholder = "Gene ID 1", 
                                                                                                        value = 'Oprk1'), 
                                                                                              textInput(inputId = "gene2", 
                                                                                                        "Enter Gene 2", 
                                                                                                        placeholder = "Gene ID 2", 
                                                                                                        value = 'Gal')), 
                                                                                 mainPanel(plotOutput('cell_types_map'),
                                                                                           plotlyOutput('moffitt1'), 
                                                                                           plotlyOutput('moffitt2'),
                                                                                           plotlyOutput('moffitt3'),
                                                                                           plotlyOutput('moffitt4')))
                       )

  
  
))))

server <- function(input, output, session) {
  
  output$delta_diff_photo <- renderImage({
    list(
      src = "Differential_Delta_Opioid_Receptor.jpg",
      width = 650,
      height = 600
    )
  }, deleteFile = FALSE)
  
  output$delta_diff_table <- renderDataTable({
    oprd1.markers})
  
  
  output$kappa_diff_photo <- renderImage({
    list(
      src = "Differential_Kappa_Opioid_Receptor.jpg",
      width = 550,
      height = 500
    )
  }, deleteFile = FALSE)
  
  output$kappa_diff_table <- renderDataTable({
    oprk1.markers})
  
  output$mu_diff_photo <- renderImage({
    list(
      src = "Differential_Mu_Opioid_Receptor.jpg",
      width = 550,
      height = 500
    )
  }, deleteFile = FALSE)
  
  output$mu_diff_table <- renderDataTable({
    oprm1.markers})
  
  
  
  
  
  output$img_1 <- renderPlotly({
    # get the id of the desired region
    granuleID = IDs[input$granuleID == IDs$name,]$id
    
    
    # get the dataset for the desired gene (the first saggital experiment that did not fail)
    datasetID = getGeneDatasets(gene = input$datasetID,
                                planeOfSection = 'coronal',
                                probeOrientation = 'antisense')[1]
    
    
    # get the slide that has the desired brain region and coordinates of the center of the region
    imageID = structureToImage(datasetID = datasetID, regionIDs = granuleID)
    
    
    # get the closest atlas image. 
    atlasID = imageToAtlas(imageID$section.image.id,imageID$x,imageID$y,planeOfSection ='coronal')
    
    # decide how much to you wish to downsample
    downsample = 2
    
    downloadAtlas(imageID = atlasID$section.image.id,
                  outputFile = 'section_atlas_anno.jpg',
                  downsample = downsample)
    
    
    fig_2 <- plot_ly()
    
    fig_2 <- fig_2 %>% layout(images = list(source = base64enc::dataURI(file = "section_atlas_anno.jpg"),
                                            xref = "x",
                                            yref = "y",
                                            x= -1,
                                            y= 4,
                                            sizex = 6,
                                            sizey = 6,
                                            opacity = 0.8))
    
  })
    
    
  output$img_2 <- renderPlotly({
    
    # get the id of the desired region
    granuleID = IDs[input$granuleID == IDs$name,]$id
    
    
    # get the dataset for the desired gene (the first saggital experiment that did not fail)
    datasetID = getGeneDatasets(gene = input$datasetID,
                                planeOfSection = 'coronal',
                                probeOrientation = 'antisense')[1]
    
    
    # get the slide that has the desired brain region and coordinates of the center of the region
    imageID = structureToImage(datasetID = datasetID, regionIDs = granuleID)
    
    
    # get the closest atlas image. 
    atlasID = imageToAtlas(imageID$section.image.id,imageID$x,imageID$y,planeOfSection ='coronal')
    
    
    # decide how much to you wish to downsample
    downsample = 2
    
    # download the slide
    downloadImage(imageID = imageID$section.image.id, 
                  view = 'projection',
                  outputFile = 'section_atlas.jpg',
                  downsample = downsample)
    
    
    
    fig <- plot_ly()
    
    
    fig <- fig %>% layout(images = list(source = base64enc::dataURI(file = "section_atlas.jpg"),
                                        xref = "x",
                                        yref = "y",
                                        x= -1,
                                        y= 4,
                                        sizex = 6,
                                        sizey = 6,
                                        opacity = 0.8))
    
    
  })
    
    
  output$img_3 <- renderPlotly({
    
    # get the id of the desired region
    granuleID = IDs[input$granuleID == IDs$name,]$id
    
    
    # get the dataset for the desired gene (the first saggital experiment that did not fail)
    datasetID = getGeneDatasets(gene = input$datasetID,
                                planeOfSection = 'coronal',
                                probeOrientation = 'antisense')[1]
    
    
    # get the slide that has the desired brain region and coordinates of the center of the region
    imageID = structureToImage(datasetID = datasetID, regionIDs = granuleID)
    
    
    # get the closest atlas image.
    atlasID = imageToAtlas(imageID$section.image.id,imageID$x,imageID$y,planeOfSection ='coronal')
    
    
    # decide how much to you wish to downsample
    downsample = 2
    
    # download the slide
    downloadImage(imageID = imageID$section.image.id,
                  view = 'projection',
                  outputFile = 'section_atlas.jpg',
                  downsample = downsample)
    
    centerImage(image = 'section_atlas.jpg',
                x = imageID$x,
                y= imageID$y,
                xProportions = c(.3,.2),
                yProportions =c(.3,.8),
                outputFile = 'section_atlas_cropped.jpg',
                downsample = downsample)
    
    
    fig_3 <- plot_ly()
    
    fig_3 <- fig_3 %>% layout(images = list(source = base64enc::dataURI(file = "section_atlas_cropped.jpg"),
                                            xref = "x",
                                            yref = "y",
                                            x= -1,
                                            y= 4,
                                            sizex = 6,
                                            sizey = 6,
                                            opacity = 0.8))
    
  })
    
  output$img_4 <- renderPlotly({
    
    # get the id of the desired region
    granuleID = IDs[input$granuleID == IDs$name,]$id
    
    
    # get the dataset for the desired gene (the first saggital experiment that did not fail)
    datasetID = getGeneDatasets(gene = input$datasetID,
                                planeOfSection = 'coronal',
                                probeOrientation = 'antisense')[1]
    
    
    # get the slide that has the desired brain region and coordinates of the center of the region
    imageID = structureToImage(datasetID = datasetID, regionIDs = granuleID)
    
    
    # get the closest atlas image. 
    atlasID = imageToAtlas(imageID$section.image.id,imageID$x,imageID$y,planeOfSection ='coronal')
    
    # decide how much to you wish to downsample
    downsample = 2
    
    downloadAtlas(imageID = atlasID$section.image.id,
                  outputFile = 'section_atlas_anno.jpg',
                  downsample = downsample)
    
    centerImage(image = 'section_atlas_anno.jpg',
                x = imageID$x,
                y= imageID$y,
                xProportions = c(.3,.2),
                yProportions =c(.3,.8),
                outputFile = 'section_atlas_anno_cropped.jpg',
                downsample = downsample)
    
    
    fig_4 <- plot_ly()
    
    fig_4 <- fig_4 %>% layout(images = list(source = base64enc::dataURI(file = "section_atlas_anno_cropped.jpg"),
                                            xref = "x",
                                            yref = "y",
                                            x= -1,
                                            y= 4,
                                            sizex = 6,
                                            sizey = 6,
                                            opacity = 0.8))
    
    fig_4
    
    
  
  })
  
  
  output$moffitt1 <- renderPlotly({
    
    # Create ggplots
    feat <- FeaturePlot(object = scrna.mast, features = c(input$gene1, input$gene2), col = c('red', 'blue'), pt.size = 0.5, blend = TRUE) 
    
    feat1 <- feat[[1]] + NoLegend()

    ggplotly(feat1, height = 350, width=400)
    
    ply1 <- ggplotly(feat1, height = 350, width=400)

    
    })
  
  output$moffitt2 <- renderPlotly({
    
    # Create ggplots
    feat <- FeaturePlot(object = scrna.mast, features = c(input$gene1, input$gene2), col = c('red', 'blue'), pt.size = 0.5, blend = TRUE) 
    
    
    feat2 <- feat[[2]] + NoLegend()
    
    ggplotly(feat2, height = 350, width=400)
    
    ply2 <- ggplotly(feat2, height = 350, width=400)
    
  })
  
  output$moffitt3 <- renderPlotly({
    
    # Create ggplots
    feat <- FeaturePlot(object = scrna.mast, features = c(input$gene1, input$gene2), col = c('red', 'blue'), pt.size = 0.5, blend = TRUE) 
    
    
    feat3 <- feat[[3]] + NoLegend()
    
    ggplotly(feat3, height = 350, width=400)
    
    ply3 <- ggplotly(feat3, height = 350, width=400)
    
  })
  
  output$moffitt4 <- renderPlotly({
    
    # Create ggplots
    feat <- FeaturePlot(object = scrna.mast, features = c(input$gene1, input$gene2), col = c('red', 'blue'), pt.size = 0.5, blend = TRUE) 
    
    
    feat4 <- feat[[4]] + NoLegend()
    
    ggplotly(feat4, height = 350, width=400)
    
    ply4 <- ggplotly(feat4, height = 350, width=400)
    
  })
  
  
  
  output$cell_types_map <- renderPlot({
    
    mycolorcode = c('yellow', 'lightblue', 'yellowgreen','purple', 'black', 'magenta', 'mediumblue', 'red', 'gray', 'darkgreen', 
                            'orange', 'blue', 'pink', 'white')
                            
    
    DimPlot(scrna.mast, group.by = 'cell.class', repel = TRUE, cols = mycolorcode) + 
      ggtitle('Cell Types in Hypothalamus Preoptic Region') + 
      guides(color=guide_legend(title="Cell Types", override.aes = list(size=4)))
    

    
  })
  
  
  
  
  
  
  
  }



shinyApp(ui, server)

