# correlation
library(corrplot)
cars <- cor(X.pool[,50:550])
corrplot(cars)
hist(X.pool[,1])

# Create a sample matrix with 0 and 1 values
mat <- matrix(sample(c(0, 1), 16, replace = TRUE), nrow = 4)

# Define a color palette with white for 0 and black for 1
palette <- c("white", "black")

# Plot the matrix using the image function and the color palette
image(X.pool[,1:10], col = palette, axes = FALSE)

# Load the ggplot2 package
library(ggplot2)
Xtx<-t(X.pool)%*%X.pool

# Assuming 'df' is your dataframe
# Calculate the correlation matrix
cor_matrix <- cor(Xtx)  # use="complete.obs" handles missing values by using complete observations

# Convert the correlation matrix into a long format for ggplot2
cor_melted <- reshape2::melt(cor_matrix)

# Draw the heatmap
ggplot(data = cor_melted, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x='', y='', title='Correlation Matrix Heatmap') +
  coord_fixed()





# Create a sample matrix with random values
set.seed(123) # Setting seed for reproducibility
matrix_data <- Xtx
library(viridis)
# Use the image function to plot the matrix
image(matrix_data, axes=FALSE, 
      col = cividis(100),
       main="Matrix Value Visualization")





matrix_data <- matrix(c(1:100), nrow=10, ncol=10)


# Convert the matrix to a long format dataframe for ggplot2
data_long <- reshape2::melt(as.data.frame(matrix_data))

# Rename the columns for clarity
colnames(data_long) <- c("X", "Y", "Value")

# Plot using ggplot2
ggplot(data_long, aes(x=X, y=Y, fill=Value)) +
  geom_tile() +
  scale_fill_gradient(low = "lightblue", high = "darkblue", name="Value") +
  labs(title="Matrix Value Visualization", x="Column", y="Row") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1))


#########Normalized
X.pool = data.matrix(read.table('./Xpool.txt', header = F))
Xtx<-t(X.pool)%*%X.pool

matrix_data <- Xtx # Create a matrix with random numbers

# Define the normalization factor (e.g., maximum absolute value, sum, etc.)
normalization_factor <- max(abs(matrix_data)) # Example: using the maximum absolute value

# Normalize the matrix
normalized_matrix <- matrix_data / normalization_factor

# View the normalized matrix
print(normalized_matrix)
library(viridis)
image(normalized_matrix, axes=FALSE, 
      col = cividis(100),
      main="Matrix Value Visualization")

cor_matrix <- cor(X.pool)
image(cor_matrix, axes=FALSE, 
      col = cividis(100),
      main="Matrix Value Visualization")