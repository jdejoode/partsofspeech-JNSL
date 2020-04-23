# Author: Johan de Joode
# Journal: Journal for North-West Semitic Languages

# POS = Part of Speech
# p.pos = Proportions of Parts of Speech

library(devtools)
# Install the plotting library if you do not already have it
install_github("vqv/ggbiplot")
library(ggbiplot)
library(rgl)
library(car)

# Set the working directory to the location of this file. 
# This is specifically put in place for work with RStudio.
# Without RStudio one should setwd manually
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read the ETCBC dataset
data = readRDS(file='./data/bhsa2017_simple.rds')

# inspect the head of the table
head(data)

# The ETCBC data also tags non-consonantal articles
View(table(data[data$pdp == 'art','g_cons_utf8']))

# Let's remove non-consonantal articles
data <- data[!((data$pdp == 'art') & (data$g_cons_utf8 != 'ה')),]

# POS per book

# I am using the phrase-dependent (context-aware) POS, rather than the simple dictionary form
# because participles really functions as nouns, they should indicate nouns and not verbs
pos <- table(data$book, data$pdp)
addmargins(pos, 2)

# Change the names
# The order matters here
POS <- c(adjv = "Adj.",
         advb = "Adv.",
         art =	"Art.",
         conj =	"Conj.",
         inrg =	"Interrog. Part.",
         intj =	"Interj.",
         nega =	"Neg. Part.",
         nmpr =	"Proper Noun",
         prde =	"Demon. Pronoun",
         prep =	"Prep.",
         prin =	"Interrog. Pronoun",
         prps =	"Pers. Pronoun",
         subs =	"Noun",
         verb =	"Verb"
)
colnames(pos) <- POS
BOOKS <- c("Amos" = "Amos",
           "Canticum" = "Song",
           "Chronica_I" = "1 Chr",
           "Chronica_II" = "2 Chr",
           "Daniel" = "Dan",
           "Deuteronomium" = "Deut",
           "Ecclesiastes" = "Eccl",
           "Esra" = "Ezra",
           "Esther" = "Esth",
           "Exodus" = "Exod",
           "Ezechiel" = "Ezek",
           "Genesis" = "Gen",
           "Habakuk" = "Hab",
           "Haggai" = "Hag",
           "Hosea" = "Hos",
           "Iob" = "Job",
           "Jeremia" = "Jer",
           "Jesaia" = "Isa",
           "Joel" = "Joel",
           "Jona" = "Jonah",
           "Josua" = "Josh",
           "Judices" = "Judg",
           "Leviticus" = "Lev",
           "Maleachi" = "Mal",
           "Micha" = "Mic",
           "Nahum" = "Nah",
           "Nehemia" = "Neh",
           "Numeri" = "Num",
           "Obadia" = "Obad",
           "Proverbia" = "Prov",
           "Psalmi" = "Ps",
           "Reges_I" = "1 Kgs",
           "Reges_II" = "2 Kgs",
           "Ruth" = "Ruth",
           "Sacharia" = "Zech",
           "Samuel_I" = "1 Sam",
           "Samuel_II" = "2 Sam",
           "Threni" = "Lam",
           "Zephania" = "Zeph")
rownames(pos) <- BOOKS
rownames(pos)

# work with proportions
p.pos <- prop.table(pos, 1)
# inspect the result
round(addmargins(p.pos, 2), 4)

# save the data for inspection
write.csv(round(p.pos *10000), './results/p_pos_per_ten_thousand.csv')
write.csv(pos, './results/pos.csv')

# Visually inspect whether the distributions are normal
par(mfrow=c(4,4))  # Set the number of graphs
for (feature in colnames(p.pos)) { 
  hist(p.pos[,feature], main=NULL, ylab=NULL, xlab=NULL)
  title(feature)
}
dev.copy2pdf(file = "./results/distributions_original.pdf")

par(mfrow=c(1,1))  # Reset the number of graphs

# Transform the data to data which PCA requires
transform2 <- function(x) {
  # --------
  # Function by Dirk Speelman
  # Box-Cox power transformation that makes x as normal as possible
  #
  # remark:
  #  we initially wanted to work with
  #    transform <- function(x, pow) { return(abs(x)^pow*sign(x))}
  #  with pow equal to 0.5 for a.o. LLR measures,
  #  but then decided to systematically use Box-Cox transformations
  # --------
  
  # -- make values strictly positive (needed for powerTransform(y)) --
  if (min(x) <= 0) {
    x <- x + 0.0000001 - min(x)
  }
  # -- bc transformation
  lambda <- powerTransform(x)$lambda
  return(bcPower(x, lambda))
}

for (f in colnames(p.pos)) {
  p.pos[,f] <- transform2(p.pos[,f])
}
    
results <- numeric(length=ncol(p.pos))
for (idx in 1:ncol(p.pos)) {
  # Test each transformed POS for normality
  results[idx] <- shapiro.test(p.pos[,idx])$p.value
}
results

# Reinspect the distributions visually
par(mfrow=c(4,4))
for (feature in colnames(p.pos)) { 
  hist(p.pos[,feature], main=NULL, ylab=NULL, xlab=NULL)
  title(feature)
}
dev.copy2pdf(file = "./results/distributions_after_boxcox.pdf")
par(mfrow=c(1,1))

# Compute PCA
pca <- prcomp(p.pos, center = TRUE, scale. = TRUE)

# Create a biplot
par(mfrow=c(1,1))
ggbiplot(pca, labels = row.names(p.pos))
ggsave(file = "./results/PCA_POS.pdf")

# Create a 3D visualisation, this needs to be screencaptured in some cases
# https://cran.r-project.org/web/packages/pca3d/pca3d.pdf
library(pca3d)
pca3d(pca, show.labels = TRUE)
makeMoviePCA(dir="./results", clean=TRUE, type = "gif", movie = "POS", convert=TRUE)

# Plot the contributions of each POS to the first three dimensions
# This needs to be saved manually as you need to redimension the window
round(pca$rotation[,1:3], 3)
par(mfrow=c(1,3), mar = c(10,4,4,2) + 0.1)
for (dimen in 1:3) {
  barplot(pca$rotation[,dimen], las=2, horiz=FALSE)
  title(paste0('PC',dimen))
}

# Plot the correlations between POS to confirm our findings
library(corrplot)
par(mfrow=c(1,1))
corrplot(cor(p.pos), order = "hclust", cl.pos = "b", col = c("white", "black"), bg = "gold2")
dev.copy2pdf(file = "./results/corrplot.pdf")