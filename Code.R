load("TCGA_HGSOC_3.RData")
summary(assoc)
summary(omicas3)
dim(GeneExpr3)
summary(subtypes3)

# Asignamos un nombre a cada ómica del conjunto de datos ómicas3.
CNV = omicas3$CNV
Methyl = omicas3$Methyl
TF = omicas3$TF
dim(CNV)
dim(Methyl)
dim(TF)
comunes.Methyl=intersect(colnames(Methyl), rownames(subtypes3))
comunes.CNV=intersect(colnames(CNV), rownames(subtypes3))
comunes.TF=intersect(colnames(TF), rownames(subtypes3))

Methyl=Methyl[comunes.Methyl]
CNV= CNV[comunes.CNV]
TF= TF[comunes.TF]

par(mfrow=c(2,2))

# Boxplot para GeneExpr3
boxplot(GeneExpr3, col='bisque', outline=F, names=NULL)
title(main='Gene Expression', col.main='blue', font.main=4)

# Boxplot para CNV
boxplot(CNV, col='lightblue', outline=F, names=NULL)
title(main='Copy Number Variation', col.main='blue', font.main=4)

# Boxplot para Methyl
boxplot(Methyl, col='salmon', outline=F, names=NULL)
title(main='Methylation', col.main='blue', font.main=4)

# Boxplot para TF
boxplot(TF, col='lightgreen', outline=F, names=NULL)
title(main='Transcription Factor', col.main='blue', font.main=4)


library(mixOmics)

GeneExpr3centrado = scale(GeneExpr3, center = TRUE, scale = FALSE)
GeneExpr3centrado_trapuesto = t(GeneExpr3centrado)
mypca = mixOmics::pca(GeneExpr3centrado_trapuesto, ncomp = 10, center = FALSE, scale = FALSE)
plot(mypca)
#Vemos que hay muy pocos genes componentes principales por debajo de 3-4. Y que principalmente hay PC-1.
mypca = mixOmics::pca(GeneExpr3centrado_trapuesto, ncomp = 2, center = FALSE, scale = FALSE)

#Expression
plotVar(mypca, comp = 1:2, var.names = FALSE, pch = 20, cex = 0.5)
plotIndiv(mypca, comp = 1:2, ind.names = NULL, group = subtypes3$subtype,
          style = 'ggplot2', legend = TRUE, legend.position = 'right', 
          ellipse = TRUE, ellipse.level = 0.95, centroid = FALSE)
#### Methyl

Methyl.centrado = scale(Methyl, center = TRUE, scale = FALSE)

Methyl.centrado_trapuesto = t(Methyl.centrado)
mypca1 = mixOmics::pca(Methyl.centrado_trapuesto, ncomp = 10, center = FALSE, scale = FALSE)

plot(mypca1)
mypca1 = mixOmics::pca(Methyl.centrado_trapuesto, ncomp = 2, center = FALSE, scale = FALSE)

plotVar(mypca1, comp = 1:2, var.names = FALSE, pch = 20, cex = 0.5)
plotIndiv(mypca1, comp = 1:2, ind.names = NULL, group = subtypes3$subtype,
          style = 'ggplot2', legend = TRUE, legend.position = 'right', 
          ellipse = TRUE, ellipse.level = 0.95, centroid = FALSE)
#### TF

TF.centrado = scale(TF, center = TRUE, scale = FALSE)

TF.centrado_trapuesto = t(TF.centrado)
mypca2 = mixOmics::pca(TF.centrado_trapuesto, ncomp = 10, center = FALSE, scale = FALSE)

plot(mypca2)
mypca1 = mixOmics::pca(TF.centrado_trapuesto, ncomp = 2, center = FALSE, scale = FALSE)

plotVar(mypca2, comp = 1:2, var.names = FALSE, pch = 20, cex = 0.5)
plotIndiv(mypca2, comp = 1:2, ind.names = NULL, group = subtypes3$subtype,
          style = 'ggplot2', legend = TRUE, legend.position = 'right', 
          ellipse = TRUE, ellipse.level = 0.95, centroid = FALSE)

#### CNV
CNV.centrado = scale(CNV, center = TRUE, scale = FALSE)

CNV.centrado_trapuesto = t(CNV.centrado)
mypca3 = mixOmics::pca(CNV.centrado_trapuesto, ncomp = 10, center = FALSE, scale = FALSE)

plot(mypca3)
mypca3 = mixOmics::pca(CNV.centrado_trapuesto, ncomp = 2, center = FALSE, scale = FALSE)

plotVar(mypca3, comp = 1:2, var.names = FALSE, pch = 20, cex = 0.5)
plotIndiv(mypca3, comp = 1:2, ind.names = NULL, group = subtypes3$subtype,
          style = 'ggplot2', legend = TRUE, legend.position = 'right', 
          ellipse = TRUE, ellipse.level = 0.95, centroid = FALSE)

###POTENCIA ESTADÍSTICA

library(MultiPower)
library(FDRsampsize)

NuevoMethyl = colnames(Methyl)= subtypes3$subtype
NuevoCNV = colnames(CNV)= subtypes3$subtype
NuevoTF = colnames(TF)= subtypes3$subtype
NuevoGeneExpr3 = colnames(GeneExpr3) = subtypes3$subtype

design = list(NuevoMethyl,NuevoTF, NuevoCNV, NuevoGeneExpr3)

omicas_unidas=list(CNV=CNV, Methyl=Methyl, TF=TF, GeneExpr3= GeneExpr3)
colores= c('lightblue', 'salmon','lightgreen', 'bisque')
names(design)= names(colores)

sapply(omicas_unidas, dim)


par(mfrow=c(1,2))
type1 = rep(2, length(omicas_unidas))
PotenciaResults=MultiPower(data= omicas_unidas, groups = design , type= type1, omicPower = 0.8, averagePower = 0.6, fdr = 0.05, cost = 1, equalSize = TRUE, max.size = 600, omicCol = colores)

### DIABLO (Multiblock PLS-DA)
# Volvemos a cargar los datos correctamente
CNV = omicas3$CNV
Methyl = omicas3$Methyl
TF = omicas3$TF

comunes.Methyl=intersect(colnames(Methyl), rownames(subtypes3))
comunes.CNV=intersect(colnames(CNV), rownames(subtypes3))
comunes.TF=intersect(colnames(TF), rownames(subtypes3))

Methyl=Methyl[comunes.Methyl]
CNV= CNV[comunes.CNV]
TF= TF[comunes.TF]


#definimos X utilizando la traspuestas para que el número de filas sea el mismo.
X<- list(CNV= scale(t(CNV), center = T, scale = F),
         TF=scale(t(TF), center = T, scale = F),
         Methyl=scale(t(Methyl), center = T, scale = F))
#vectorizamos el data frame que contiene los subtipos de cáncer
Y<- unlist(as.vector(lapply(subtypes3,as.factor)), use.names = F)

#realizamos el diseño de nuestro estudio para dárselo como argumento a la función block.splsda
design.diablo<- matrix(0.1, ncol = length(X), nrow = length(X), 
                dimnames = list(names(X), names(X)))
diag(design.diablo) <- 0
design.diablo 

#Esta función tarda aproximadamente 1'5 horas, por lo que se puede omitir.
#tuneBlock=tune.block.splsda(X = X,Y = Y, ncomp=2, design = design.diablo)

ncomp = 2
list.keepX = list(CNV = c(20, 15), Methyl = c(10,5), TF = c(15, 25))

Diablo.splsda = mixOmics::block.splsda(X = X, 
                                           Y = Y, 
                                           scale = FALSE, ncomp = ncomp, 
                                           keepX = list.keepX, design = design.diablo)


plotLoadings(Diablo.splsda, comp = 1, contrib = 'max')
plotLoadings(Diablo.splsda, comp = 2, contrib = 'max')

plotVar(Diablo.splsda, style = 'graphics', legend = TRUE, cutoff = 0.5, var.names = T)


plotIndiv(Diablo.splsda, ind.names = FALSE, legend = TRUE, ellipse = T, ellipse.level = 0.95, legend.position = 'bottom')

plotDiablo(Diablo.splsda, ncomp = 1, legend.ncol = T, legend = T)
plotDiablo(Diablo.splsda, ncomp = 2, legend.ncol = T, legend=T)

circosPlot(Diablo.splsda, comp = 1, cutoff = 0.7, size.variables = 0.60)


#este script hay que copiarlo y ejecutarlo en la consola
X11()
cimDiablo(Diablo.splsda, size.legend = 1, margins = c(8,20))

#este script hay que copiarlo y ejecutarlo en la consola
set.seed(123)
X11()
network.res=network(Diablo.splsda, blocks = c(1,2,3),
        cutoff = 0.7, alpha.node = 0.5, 
        lwd.edge = 2, show.edge.labels = T,
        show.color.key=T)
##### Train Data
# Predicción del Modelo
pred <- predict(Diablo.splsda, newdata = X)
# Subtipo Original
message("Original subtype:")
head(subtypes3)

# Subtipo predicho
message("Predicted subtype:")
head(pred$AveragedPredict.class$max.dist)

message("Prediction Confusion Matrix:")

# Matriz de Confusión

confusion_matriz <- table(data.frame("prediction" = pred$AveragedPredict.class$max.dist[, 2],
                                     "trueType" = subtypes3), useNA = "ifany")

# Calcular la matriz de confusión en porcentajes
confusion_matriz_porcentaje <- prop.table(confusion_matriz, margin = 1) * 100

# Las mostramos
confusion_matriz
confusion_matriz_porcentaje

library(igraph)
x11()
myNetwork <- network(Diablo.splsda, blocks=c(1,2,3), cutoff=0.7)
write_graph(myNetwork$gR, file = "diablo.gml", format = "gml")

