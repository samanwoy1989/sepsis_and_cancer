library(grid)
library(useful)
library(magick)

# Read external images
imageA <- image_read("Fig2A.jpeg") 
imageB <- image_read("Fig2B.jpeg") 

# Create data frame
df <- data.frame(1:2)

# Create PDF
pdf("Mydocument.pdf", width = 10, height = 20)
grid.newpage() 

# Create matrix layout
pushViewport(viewport(layout = grid.layout(2,1)))

# Place elements inside grid
print(imageA, vp = vplayout(1, 1)) 
print(imageB, vp = vplayout(2, 1))
#print(df, vp = vplayout(1, 3)) 

dev.off()
