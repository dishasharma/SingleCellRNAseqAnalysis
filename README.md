# SingleCellRNAseqAnalysis

1. OutputName <- scRNASeurat(QCedRDS = ObjectName, Dims = Dims ,saveDir = "OutputDir" ,Assay = "RNA",samplename = "SampleName", condition = "Condition")
2. OutputName <- IntegrationSeurat(QCedRDS = ObjectName, Dims = Dims ,saveDir = "OutputDir" ,Assay = "RNA", samplename = "SampleName", condition = "Condition", resolution = 0.5, splitby = "column to split by")
3. logisticregression(table = input table, labelcolumn = "label column", label1 = "class1", label2 = "class2", featurecolumn = "Gene to identify as marker")
