#A) Module score quantitatif – “Rupture-protective signature”
#1) Définition des modules (à copier-coller)
rupture_protective_modules <- list(
  ImmediateEarly = c("Nr4a1","Egr1","Ccn1","Btg2","Dusp1"),
  Contractile_SMC = c("Myocd","Acta2","Tagln","Cnn1","Cnn3","Myl9","Actg2"),
  Proteostasis = c("Hspa1a","Hspa1b","Hspa8","Hsp90aa1","Hsp90ab1","Hspb1","Bag3"),
  Mechanotransduction = c("Itga5","Itga8","Fermt2","Rock1","Piezo1"),
  MetabolicSupport = c("Pdk4","Pfkfb3","Ppargc1a")
)
rupture_protective_modules
# Astuce : adaptez la casse des gènes si nécessaire (mouse symbols).
dataBl6.Balbc
table(Idents(data.integrated2))
data.integrated2$celltype <- Idents(data.integrated2)
table(dataBl6.Balbc$celltype)
celltype_of_interest <- c("EC vas", "EC lymph/vas")
data_EC <- subset(dataBl6.Balbc, subset = celltype == celltype_of_interest)
data_EC
# Extraction de DEG exprimés entre KI et WT (chez les mâles)

Idents(data_EC) <- "genotype"  
levels(Idents(data_EC))
markers_male_KI_vs_WT <- FindMarkers(
  data_EC,
  ident.1 = "KI",
  ident.2 = "WT",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
head(markers_male_KI_vs_WT)
write.xlsx(markers_male_KI_vs_WT,
           "/home/afall/r_docker/Analyses_SnRNA_10_2025/Sans_Balbc/Version_nautilus/Analyse_différentielle_types_cellulaires/Mâles/EC/DEG_EC_KI_vs_WT_MALE.xlsx",
           rowNames = TRUE)

library(ggrepel)

markers_male_KI_vs_WT$gene <- rownames(markers_male_KI_vs_WT)

# Création colonne "signification" selon les seuils
markers_male_KI_vs_WT$signification <- "NS"
markers_male_KI_vs_WT$signification[markers_male_KI_vs_WT$p_val_adj < 0.05 &
                                      markers_male_KI_vs_WT$avg_log2FC > 0.25] <- "Up in KI"
markers_male_KI_vs_WT$signification[markers_male_KI_vs_WT$p_val_adj < 0.05 &
                                      markers_male_KI_vs_WT$avg_log2FC < -0.25] <- "Up in WT"
# Sélectionner les gènes à annoter (les 20 plus différentiels)
genes_to_label <- subset(markers_male_KI_vs_WT, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)
genes_to_label <- genes_to_label[order(abs(genes_to_label$avg_log2FC), decreasing = TRUE), ][1:20, ]
genes_to_label
# Tous les points = markers_femelle_KI_vs_WT
volcano <- ggplot(markers_male_KI_vs_WT, aes(x = avg_log2FC, y = -log10(p_val_adj), color = signification)) +
  geom_point(alpha = 0.7, size = 1.5) +
  
  # Annotations = seulement diff_genes
  geom_text_repel(data = genes_to_label, aes(label = gene), size = 3, max.overlaps = 20) +
  
  scale_color_manual(values = c("Up in KI" = "red", "Up in WT" = "blue", "NS" = "grey70")) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Volcano plot: KI vs WT (Male EC)",
       x = "log2 Fold Change",
       y = "-log10(p_adj)") +
  theme_classic(base_size = 14)
volcano
library(Seurat)

dataBl6.Balbc <- AddModuleScore(
  object = dataBl6.Balbc,
  features = rupture_protective_modules,
  name = names(rupture_protective_modules)
)
dataBl6.Balbc@meta.data
# Cela crée dans meta.data : ImmediateEarly1, Contractile_SMC2, Proteostasis3, Mechanotransduction4, MetabolicSupport5

# 3) Score global “Rupture-Protective Index” (RPI)
dataBl6.Balbc$Rupture_Protective_Index <- rowMeans(
  dataBl6.Balbc@meta.data[, c("ImmediateEarly1",
                      "Contractile_SMC2",
                      "Proteostasis3",
                      "Mechanotransduction4",
                      "MetabolicSupport5")],
  na.rm = TRUE
)
# 4) Visualisations recommandées (simples et robustes)
# a) Violin plot – par fond génétique et génotype
dataBl6.Balbc$genotype <- ifelse(grepl("WT", dataBl6.Balbc$orig.ident), "WT", "KI")
table(dataBl6.Balbc$genotype)
unique(dataBl6.Balbc$orig.ident)
V1 <- VlnPlot(
  dataBl6.Balbc,
  features = "Rupture_Protective_Index",
  group.by = "genotype",
  split.by = "orig.ident",
  pt.size = 0,
  cols = c("steelblue", "firebrick", "darkgreen", "orange")
)
ggsave(filename = "/home/afall/r-docker/Analyses_SnRNA_10_2025/Sans_Balbc/Version_nautilus/Analyses_SnRNA_01_2026/Demande_26_01_2026/Vlnplot_fond_génétique_B6Balbc_res005_new.png", plot = V1, width = 12, height = 7, dpi = 300)

# Attendu :
 # • BALB/c KI > BALB/c WT
 # • BALB/c KI >> BL6J KI
 # • BL6J KI peu ou pas augmenté

# b) Dot plot – par type cellulaire
colnames(dataBl6.Balbc@meta.data)
levels(Idents(dataBl6.Balbc))
dataBl6.Balbc$celltype <- Idents(dataBl6.Balbc)
#D1 <- DotPlot(
  #dataBl6.Balbc,
  #features = c(
    #"ImmediateEarly1","Contractile_SMC2",
   # "Proteostasis3","Mechanotransduction4",
   # "MetabolicSupport5"
  #),
 # group.by = "celltype",
  #split.by = "orig.ident",
  #cols = c("lightgrey", "blue", "red", "darkred")
#) + RotatedAxis()

D1 <- DotPlot(
  dataBl6.Balbc,
  features = c(
    "ImmediateEarly1","Contractile_SMC2",
    "Proteostasis3","Mechanotransduction4",
    "MetabolicSupport5"
  ),
  group.by = "celltype",
  split.by = "orig.ident",
  cols = c("lightgrey", "blue", "red", "darkred")
) +
  RotatedAxis() +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 80,
      hjust = 1,
      vjust = 1,
      margin = margin(t = 5)
    ),
    axis.text.y = element_text(size = 9),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10)
  )
D1
table(dataBl6.Balbc$orig.ident)
ggsave(filename = "/home/afall/r-docker/Analyses_SnRNA_10_2025/Sans_Balbc/Version_nautilus/Analyses_SnRNA_01_2026/Demande_12_01_2026/Dotplot_fond_génétique_B6Balbc_res005.png", plot = D1, width = 12, height = 7, dpi = 300)

# Attendu :
  # • Signal maximal dans les CML contractiles BALB/c
  # • Faible ou absent dans fibroblastes BL6J
dataBl6.Balbc
 # c) UMAP – score global
F1 <- FeaturePlot(
  dataBl6.Balbc,
  features = "Rupture_Protective_Index",
  split.by = "orig.ident",
  cols = c("white","red")
)
F2 <- FeaturePlot(
  dataBl6.Balbc,
  features = "Rupture_Protective_Index",
  split.by = "orig.ident",
  cols = c("navy", "yellow", "red")  # gradient froid → chaud
)
ggsave(filename = "/home/afall/r-docker/Analyses_SnRNA_10_2025/Sans_Balbc/Version_nautilus/Analyses_SnRNA_01_2026/Demande_12_01_2026/FeaturePlot_fond_génétique_B6Balbc_res005_esssai.png", plot = F2, width = 12, height = 7, dpi = 300)

# 5) Test statistique simple (publication-safe)
library(dplyr)
library(rstatix)

  sample_means <- dataBl6.Balbc@meta.data %>%
  group_by(orig.ident, genotype) %>%
  summarise(
    mean_RPI = mean(Rupture_Protective_Index),
    .groups = "drop"
  )

sample_means
wilcox.test(mean_RPI ~ genotype, data = sample_means)
library(lme4)
library(lmerTest)

model <- lmer(
  Rupture_Protective_Index ~ genotype + (1 | orig.ident),
  data = dataBl6.Balbc@meta.data
)
library(openxlsx)

# Extraire uniquement les coefficients fixes (estimate, SE, t-value, p-value)
T1_fixed <- as.data.frame(coef(summary(model)))

T1_fixed
#T1 <- dataBl6.Balbc@meta.data %>%
  #group_by(orig.ident, genotype) %>%
  #summarise(mean_RPI = mean(Rupture_Protective_Index),
            #sd_RPI = sd(Rupture_Protective_Index),
            #n = n())
write.xlsx(T1_fixed, "/home/afall/r-docker/Analyses_SnRNA_10_2025/Sans_Balbc/Version_nautilus/Analyses_SnRNA_01_2026/Demande_26_01_2026/Test_significativité_statistique_B6Balbc_res005_new.xlsx", row.names = FALSE)

file_path <- "/home/afall/r-docker/Analyses_SnRNA_10_2025/Sans_Balbc/Version_nautilus/Analyses_SnRNA_01_2026/Demande_26_01_2026/Summary_model_test_significativité_statistique_B6Balbc_res005.txt"

# Enregistrer tout le summary
sink(file_path)
summary(model)
sink()  # fermeture du fichier

##### Demande 26/01/2026 Xavier 

#1. Définir « signature endotheliale »
Endothelial_Markers <- c(
  "Pecam1",   # CD31
  "Cdh5",     # VE-cadherin
  "Kdr",      # VEGFR2
  "Flt1",     # VEGFR1
  "Nos3",     # eNOS
  "Vwf",
  "Klf2",
  "Klf4",
  "Tek",      # Tie2
  "Esam"
)
#2. Calculer un « Endothelium Module Score »

dataB6M <- AddModcore(
  object = dataB6M,
  features = list(Endothelial_Markers),
  name = "EC_Score"
)

#Cela crée une colonne :  EC_Score1
colnames(dataB6M@meta.data)
#3. Comparer Bl6 WT vs KI
V2 <- VlnPlot(
  dataB6M,
  features = "EC_Score1",
  group.by = "genotype",
  split.by = "orig.ident",
  pt.size = 0
)
ggsave(filename = "/home/afall/r-docker/Analyses_SnRNA_10_2025/Sans_Balbc/Version_nautilus/Analyses_SnRNA_01_2026/Demande_26_01_2026/Demandes_d'analyse_Xavier/Vlnplot_fond_génétique_EC_Bl6M.png", plot = V2, width = 12, height = 7, dpi = 300)
#Faire cette même comparaison dans le fond génétique BL6J entre WT et KI : On s’attend voir des différences !
#Normalement, on ne s’attend pas à des différences de signature endothéliale dans ce fond génétique BALBC entre WT et KI

# Test statistique sur le vlnplot dataB6M:
df <- FetchData(dataB6M, vars = c("EC_Score1", "genotype"))
wilcox.test(EC_Score1 ~ genotype, data = df)
       ###### Résultat du test: 

       #Wilcoxon rank sum test with continuity correction

       #data:  EC_Score1 by genotype
       #W = 39841872, p-value < 2.2e-16
       #alternative hypothesis: true location shift is not equal to 0

#### Meme analyse précédente avec les Balbc
dataBalcM <- AddModuleScore(
  object = dataBalcM,
  features = list(Endothelial_Markers),
  name = "EC_Score"
)
colnames(dataBalcM@meta.data)
#3. Comparer BALB/c WT vs KI
dataBalcM$genotype <- ifelse(grepl("WT", dataBalcM$orig.ident), "WT", "KI")
table(dataBalcM$genotype)
Idents(dataBalcM) <- "genotype" 
V3 <- VlnPlot(
  dataBalcM,
  features = "EC_Score1",
  group.by = "genotype",
  split.by = "orig.ident",
  pt.size = 0
)
ggsave(filename = "/home/afall/r-docker/Analyses_SnRNA_10_2025/Sans_Balbc/Version_nautilus/Analyses_SnRNA_01_2026/Demande_26_01_2026/Demandes_d'analyse_Xavier/Vlnplot_fond_génétique_EC_BalcM.png", plot = V3, width = 12, height = 7, dpi = 300)
# Test statistique sur le vlnplot dataBalcM:
df <- FetchData(dataBalcM, vars = c("EC_Score1", "genotype"))
wilcox.test(EC_Score1 ~ genotype, data = df)
           ###### Résultat du test: 

           #Wilcoxon rank sum test with continuity correction

           #data:  EC_Score1 by genotype
           #W = 52663334, p-value < 2.2e-16
           #alternative hypothesis: true location shift is not equal to 0

### Analyse des signatures
#CML Contractile Signature

SMC_Contractile <- c(
  "Myocd",
  "Acta2",
  "Tagln",
  "Cnn1",
  "Cnn3",
  "Myl9",
  "Actg2",
  "Pdlim3",
  "Ldb3",
  "Tpm1"
)

# CML Modulated / Synthetic Signature

SMC_Modulated <- c(
  "Klf4",
  "Spp1",
  "Fn1",
  "Col1a1",
  "Mmp2",
  "Lgals3",
  "Thbs1",
  "Serpine1",
  "Pdgfrb",
  "Ccn2"
)

# Fibroblast ECM Signature

Fibro_ECM <- c(
  "Col1a1",
  "Col1a2",
  "Col3a1",
  "Dcn",
  "Lum",
  "Fbln1",
  "Postn",
  "Lox",
  "Pdgfra",
  "Tcf21"
)


# Fibroblast Activated/Remodeling Signature

Fibro_Activated <- c(
  "Ccn2",
  "Crispld2",
  "Serpine1",
  "Adamts1",
  "Mmp14",
  "Thbs4",
  "Tgfb1i1",
  "Ctgf",
  "Sfrp4",
  "Bmp6"
)

#Les calculs dans Seurat pour l'analyse des signatures: pour B6M

dataB6M <- AddModuleScore(dataB6M, list(SMC_Contractile), name="SMC_Contractile")
colnames(dataB6M@meta.data)
dataB6M <- AddModuleScore(dataB6M, list(SMC_Modulated), name="SMC_Modulated")
dataB6M <- AddModuleScore(dataB6M, list(Fibro_ECM), name="Fibro_ECM")
dataB6M <- AddModuleScore(dataB6M, list(Fibro_Activated), name="Fibro_Activated")
table(dataB6M$genotype)
# Comparaison entre module
unique(dataB6M$orig.ident)
V4 <- VlnPlot(
  dataB6M,
  features = c(
    "SMC_Contractile1",
    "SMC_Modulated1",
    "Fibro_ECM1",
    "Fibro_Activated1"
  ),
  group.by = "genotype",
  pt.size = 0
)
ggsave(filename = "/home/afall/r-docker/Analyses_SnRNA_10_2025/Sans_Balbc/Version_nautilus/Analyses_SnRNA_01_2026/Demande_26_01_2026/Demandes_d'analyse_Xavier/Vlnplot_signatures_B6M.png", plot = V4, width = 12, height = 7, dpi = 300)

#Les calculs dans Seurat pour l'analyse des signatures: pour BalbcM

dataBalcM <- AddModuleScore(dataBalcM, list(SMC_Contractile), name="SMC_Contractile")
colnames(dataBalcM@meta.data)
dataBalcM <- AddModuleScore(dataBalcM, list(SMC_Modulated), name="SMC_Modulated")
dataBalcM <- AddModuleScore(dataBalcM, list(Fibro_ECM), name="Fibro_ECM")
dataBalcM <- AddModuleScore(dataBalcM, list(Fibro_Activated), name="Fibro_Activated")

# Comparaison entre module
unique(dataBalcM$orig.ident)
V5 <- VlnPlot(
  dataBalcM,
  features = c(
    "SMC_Contractile1",
    "SMC_Modulated1",
    "Fibro_ECM1",
    "Fibro_Activated1"
  ),
  group.by = "genotype",
  pt.size = 0
)
V5
ggsave(filename = "/home/afall/r-docker/Analyses_SnRNA_10_2025/Sans_Balbc/Version_nautilus/Analyses_SnRNA_01_2026/Demande_26_01_2026/Demandes_d'analyse_Xavier/Vlnplot_signatures_BalbcM.png", plot = V5, width = 12, height = 7, dpi = 300)

# D'autre analyse de module: demande du 01/02/2026

#1 BL6J VSMC Modulation Module Score

#This captures the pathological phenotypic switch seen in BL6J.

BL6J_VSMC_Modulation <- c(
  "Egr1","Atf3","Nr4a3",
  "Thbs1","Fst","Cytl1","Igfbp6","Retnla",
  "Ccl24","Fcna",
  "Sgpp2",
  "Crispld1"
)

#Calculate:
  
dataB6M <- AddModuleScore(dataB6M,
                          features = list(BL6J_VSMC_Modulation),
                          name = "BL6J_Modulation")

 # 2 BALB/c Protective VSMC Module Score

BALBc_Protective_VSMC <- c(
  "Myocd","Acta2","Tagln","Cnn1","Cnn3","Myl9",
  "Hspa1a","Hspa1b","Hspb1","Bag3",
  "Itga5","Itga8","Rock1","Piezo1"
)

dataBalcM <- AddModuleScore(dataBalcM,
                        features = list(BALBc_Protective_VSMC),
                        name = "BALBc_Protection")

# 3 Single Mechanistic Balance Index  (au final la partie 3 n'a pas été réalisé parce que ce sode ne peux pas marcher:   aorta$VSMC_Fate_Index <- dataBalcM$BALBc_Protection1 - dataB6M$BL6J_Modulation1)

#This is extremely powerful:
  
  aorta$VSMC_Fate_Index <- dataBalcM$BALBc_Protection1 - dataB6M$BL6J_Modulation1

#Interpretation
# Index value Meaning
#High (+)   Stable contractile compensation (BALB/c-like)
#Low (−)    Pathological modulation (BL6J-like)

  #Visualisation

#VlnPlot(aorta,
       # features = "VSMC_Fate_Index",
       # group.by = "genotype",
       # split.by = "strain",
        #pt.size = 0)

# Expected:
  
# BALB/c KI → strongly positive

# BL6J KI → strongly negative

# Pourriez-vous aussi  projeter le  VSMC_Fate_Index sur l’UMAP ?
  
################ Pour les EC entre bl6j et balbc 
  
# 1 Define the EC stress gene set
  
#Based on your BL6J EC DEG list (and vascular biology literature), use:
    
    EC_Stress_Genes <- c(
      "Fos",
      "Egr1",
      "Serpine1",
      "Jun",
      "Dusp1",
      "Atf3",
      "Klf6",
      "Nr4a1",
      "Ccn1"
    )
  
  #These represent:
    
    #Immediate-early activation
  
    #Mechanical/shear stress response
  
    #Endothelial remodeling signaling
    
    ### pour les B6M  
    
    #2 Compute the EC Stress Score (Seurat)
    
    dataB6M <- AddModuleScore(
      object = dataB6M,
      features = list(EC_Stress_Genes),
      name = "EC_Stress"
    )
    
    #This creates:
    
    dataB6M$EC_Stress1
    
    #3 Compare WT vs KI, BL6J vs BALB/c
    
    V6 <- VlnPlot(
      dataB6M,
      features = "EC_Stress1",
      group.by = "genotype",
      split.by = "orig.ident",
      pt.size = 0
    )
    ggsave(filename = "/home/afall/r-docker/Analyses_SnRNA_10_2025/Sans_Balbc/Version_nautilus/Analyses_SnRNA_01_2026/Demande_26_01_2026/Demandes_d'analyse_Xavier/Vlnplot_EC_Stress_B6M.png", plot = V6, width = 12, height = 7, dpi = 300)
    
    #️4 Visualize on UMAP
    
    F3 <- FeaturePlot(
      dataB6M,
      features = "EC_Stress1",
      split.by = "orig.ident",
      cols = c("navy", "yellow", "red")
    )
    ggsave(filename = "/home/afall/r-docker/Analyses_SnRNA_10_2025/Sans_Balbc/Version_nautilus/Analyses_SnRNA_01_2026/Demande_26_01_2026/Demandes_d'analyse_Xavier/Featureplot_EC_Stress_B6M.png", plot = F3, width = 12, height = 7, dpi = 300)
    #You should see signal restricted to endothelial clusters, especially in BL6J KI.
    #5 Statistical comparison
    
    library(dplyr)
    
    dataB6M@meta.data %>%
      group_by(orig.ident, genotype) %>%
      summarise(mean = mean(EC_Stress1),
                sd = sd(EC_Stress1))
### pour les Balbc 
  
#2 Compute the EC Stress Score (Seurat)
  
    dataBalcM <- AddModuleScore(
    object = dataBalcM,
    features = list(EC_Stress_Genes),
    name = "EC_Stress"
  )
  
  #This creates:
    
    dataBalcM$EC_Stress1
  
#3 Compare WT vs KI, BL6J vs BALB/c
  
  V7 <- VlnPlot(
    dataBalcM,
    features = "EC_Stress1",
    group.by = "genotype",
    split.by = "orig.ident",
    pt.size = 0
  )
ggsave(filename = "/home/afall/r-docker/Analyses_SnRNA_10_2025/Sans_Balbc/Version_nautilus/Analyses_SnRNA_01_2026/Demande_26_01_2026/Demandes_d'analyse_Xavier/Vlnplot_EC_Stress_BalbcM.png", plot = V7, width = 12, height = 7, dpi = 300)
  
#️4 Visualize on UMAP
  
 F4 <- FeaturePlot(
   dataBalcM,
    features = "EC_Stress1",
    split.by = "orig.ident",
   cols = c("navy", "yellow", "red")
  )
 
table(dataBalcM$genotype) # Pour connaitre le nbre de noyau qe contient chaque échantillon WT et KI
ggsave(filename = "/home/afall/r-docker/Analyses_SnRNA_10_2025/Sans_Balbc/Version_nautilus/Analyses_SnRNA_01_2026/Demande_26_01_2026/Demandes_d'analyse_Xavier/Featureplot_EC_Stress_BalbcM.png", plot = F4, width = 12, height = 7, dpi = 300)
#You should see signal restricted to endothelial clusters, especially in BL6J KI.
#5 Statistical comparison
  
  library(dplyr)
  
dataBalcM@meta.data %>%
    group_by(orig.ident, genotype) %>%
    summarise(mean = mean(EC_Stress1),
              sd = sd(EC_Stress1))
  
  #Then Wilcoxon or ANOVA.
 #Biological interpretation
  
  #High EC-Stress score means:
    
  #endothelial cells sense abnormal wall mechanics
  
  #but does not imply endothelial failure
  
  #It measures reactive activation, not dysfunction.
  #Optional refinement
  
  #To show ECs are activated but stable, you can compute also:
    
 #   EC_Identity <- c("Pecam1","Cdh5","Kdr","Tek","Esam","Vwf")
  
  #Stable identity + high stress = reactive but intact endothelium.

    #Expected result
    #Group
    
    #EC stress score
    #BL6J KI
    #↑↑
    #BL6J WT
    
    #low
    #BALB/c KI
    
    #slight ↑ or stable
    #BALB/c WT
    #low



#### Analyses démandées par Xavier le 16/02/2026

#1 Mitochondrial Dysfunction Score

#Gene set (mouse)

Mito_Genes <- c(
  "mt-Nd1","mt-Nd2","mt-Nd3","mt-Nd4","mt-Nd4l","mt-Nd5","mt-Nd6",
  "mt-Co1","mt-Co2","mt-Co3",
  "mt-Atp6","mt-Atp8",
  "Cox6a1","Cox6c","Cox7a2","Cox7c",
  "Ndufb1","Ndufb3","Ndufb9"
)
#Compute score

dataB6M <- AddModuleScore(
  dataB6M,
  features = list(Mito_Genes),
  name = "Mito_Score"
)
#If mitochondrial genes are downregulated in KI, then:
  
#BL6J KI → lower Mito_Score

#BL6J WT → higher

#2 Correlation with VSMC Modulation Score

#Now test whether mitochondrial loss associates with phenotypic switching.
cor.test(dataB6M$BL6J_Modulation1,# cette commande ne marche pas parceque nous ne l'avons pas défini au préalable 
         dataB6M$Mito_Score1,
         method = "spearman")
#Interpretation:
  
#Correlation
#Meaning
#Strong negative
#More modulation → less mitochondrial integrity
#Weak
#Independent processes
#Positive
#Would contradict current hypothesis
#Expected in BL6J KI:
  #→ Negative correlation

#That would mean:
  
  #Energetic collapse accompanies phenotypic modulation.

#3 Cell-type specificity check

#Now we determine whether mitochondrial decline is:
  
  #VSMC-specific

#or global across cell types

#Plot by cluster
dataB6M$celltype <- Idents(dataB6M)
table(dataB6M$celltype)
V8 <- VlnPlot(
  dataB6M,
  features = "Mito_Score1",
  group.by = "celltype",
  split.by = "genotype",
  pt.size = 0
)
ggsave(filename = "/home/afall/r-docker/Analyses_SnRNA_10_2025/Sans_Balbc/Version_nautilus/Analyses_SnRNA_01_2026/Demande_06_02_2026/Demandes_d'analyse_Xavier_suite/Vlnplot_mitochondrial_score_B6M.png", plot = V8, width = 12, height = 7, dpi = 300)

#If decrease is strongest in:
  
  #VSMCs → intrinsic medial metabolic stress

#ECs/fibroblasts also ↓ → systemic mitochondrial stress


#### Test cluster 5 d'adipocyte res 0.5 de dataB6M

#Test 1: Co-expression with VSMC core genes

#Check whether cluster 5 cells express:
  
  VSMC_Core <- c("Acta2","Tagln","Myh11","Cnn1","Myocd")
#Run:
  
  D2 <- DotPlot(dataB6M,
          features = VSMC_Core,
          group.by = "seurat_clusters") + RotatedAxis()+
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 80, hjust = 1, vjust = 1, margin = margin(t = 5)))
  ggsave(filename = "/home/afall/r-docker/Analyses_SnRNA_10_2025/Sans_Balbc/Version_nautilus/Analyses_SnRNA_01_2026/Demande_06_02_2026/Demandes_d'analyse_Xavier_suite/Dotplot_origine_adipo_cluster5_B6M_res05_markersCML.png", plot = D2, width = 12, height = 7, dpi = 300)
  
#If cluster 5:
  
  #expresses at least some Acta2/Tagln → likely VSMC-derived

#expresses none → likely non-VSMC lineage

#Test 2 : True Adipocyte Contamination Check

#True adipocytes express:
  
  Adipocyte_Markers <- c("Adipoq","Fabp4","Lpl","Leptin","Cidec")
#Run:
  Reductions(dataB6M)
  
  F5 <- FeaturePlot(
    dataB6M,
    features = Adipocyte_Markers,
    reduction = "umap"
  )
  F5
  ggsave(filename = "/home/afall/r-docker/Analyses_SnRNA_10_2025/Sans_Balbc/Version_nautilus/Analyses_SnRNA_01_2026/Demande_06_02_2026/Demandes_d'analyse_Xavier_suite/FeaturePlot_origine_adipo_cluster5_B6M_res05_markersAdipo.png", plot = F5, width = 12, height = 7, dpi = 300)
  
  D3 <- DotPlot(dataB6M,
                features = Adipocyte_Markers,
                group.by = "seurat_clusters") + RotatedAxis()+
                theme_classic() + 
    theme(axis.text.x = element_text(angle = 80, hjust = 1, vjust = 1, margin = margin(t = 5)))
  D3
  ggsave(filename = "/home/afall/r-docker/Analyses_SnRNA_10_2025/Sans_Balbc/Version_nautilus/Analyses_SnRNA_01_2026/Demande_06_02_2026/Demandes_d'analyse_Xavier_suite/Dotplot_origine_adipo_cluster5_B6M_res05_markersAdipo.png", plot = D3, width = 12, height = 7, dpi = 300)
  
#If these are absent in cluster 5:
#  → Not mature adipocytes.
  dot <- DotPlot(dataB6M, features = features) + 
    theme_classic() + #  Fond blanc, sans quadrillage ni bordures. Permet d'enlever les carreaux gris derrière mon image comme nous pouvons les voir dans dotplot1 
    theme(axis.text.x = element_text(angle = 80, hjust = 1, vjust = 1, margin = margin(t = 5)))
  dot


