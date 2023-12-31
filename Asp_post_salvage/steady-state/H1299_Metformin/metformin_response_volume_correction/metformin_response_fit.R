
# Load data:
df <- read.csv('/Users/krdav/Desktop/Asp_post_salvage/metformin_response_volume_correction/metformin_response.csv', sep = ',')
# Generate interaction between response (peak area) and volume (volume sample dried):
df$volxR <- df$Response / df$Vol

# Generate linear model to predict Metformin concentration from the response.
# Include qubic term due to ion suppression:
linmod <- lm(Conc ~ Response + volxR + I(Response^2) + I(volxR^2), data = df)
summary(linmod)

# Plot before:
plot(y = df$Conc, x = df$Response)

# Plot prediction:
pred <- predict(object = linmod)
plot(y = df$Conc, x = pred)

# Log transform and plot again:
plot(y = log(df$Conc[df$Conc > 0]), x = log(pred[df$Conc > 0]))






### Subset based on cell line and salvage mix ###

### 143B No salvage mix
df1 <- df[df$Cell_line == '143B' & df$Salvage_mix == 'No', ]

linmod <- lm(Conc ~ Response+volxR+I(Response^2)+I(volxR^2), data = df1)
summary(linmod)

plot(y = df1$Conc, x = df1$Response)
pred <- predict(object = linmod)
plot(y = df1$Conc, x = pred)



### 143B with salvage mix
df1 <- df[df$Cell_line == '143B' & df$Salvage_mix == 'Yes', ]

linmod <- lm(Conc ~ Response+volxR+I(Response^2)+I(volxR^2), data = df1)
summary(linmod)

plot(y = df1$Conc, x = df1$Response)
pred <- predict(object = linmod)
plot(y = df1$Conc, x = pred)



### H1299 No salvage mix
df1 <- df[df$Cell_line == 'H1299' & df$Salvage_mix == 'No', ]

linmod <- lm(Conc ~ Response+volxR+I(Response^2)+I(volxR^2), data = df1)
summary(linmod)

plot(y = df1$Conc, x = df1$Response)
pred <- predict(object = linmod)
plot(y = df1$Conc, x = pred)




### H1299 with salvage mix
df1 <- df[df$Cell_line == 'H1299' & df$Salvage_mix == 'Yes', ]

linmod <- lm(Conc ~ I(Response^2)+I(volxR^2), data = df1)
summary(linmod)

plot(y = df1$Conc, x = df1$Response)
pred <- predict(object = linmod)
plot(y = df1$Conc, x = pred)










