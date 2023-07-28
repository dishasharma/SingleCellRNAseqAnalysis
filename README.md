# SingleCellRNAseqAnalysis

OutputName <- scRNASeurat(QCedRDS = ObjectName, Dims = Dims ,saveDir = "OutputDir" ,Assay = "RNA",samplename = "SampleName", condition = "Condition")
OutputName <- IntegrationSeurat(QCedRDS = ObjectName, Dims = Dims ,saveDir = "OutputDir" ,Assay = "RNA", samplename = "SampleName", condition = "Condition", resolution = 0.5, splitby = "column to split by")
logisticregression(table = input table, labelcolumn = "label column", label1 = "class1", label2 = "class2", featurecolumn = "Gene to identify as marker")
