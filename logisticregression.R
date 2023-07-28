logisticregression <- function(table, labelcolumn, label1 , label2, featurecolumn){
table[,labelcolumn] <- gsub(label1,0,table[,labelcolumn])
table[,labelcolumn] <- gsub(label2,1,table[,labelcolumn])
plot(table[,featurecolumn],table[,labelcolumn])

print(unique(table[,labelcolumn]))
plot.new()
glm_fit=glm(as.numeric(table[,labelcolumn]) ~ table[,featurecolumn], family=binomial, data = table)
lines(table[,featurecolumn], glm_fit$fitted.values)
test = predict(glm_fit,type = "response") >= 0.5
test2 <- gsub(TRUE,1,as.vector(test))
test2 <- gsub(FALSE,0,as.vector(test2))

par(pty="s")
roc.info  <- roc(as.numeric(table[,labelcolumn]),glm_fit$fitted.values,plot = TRUE)
roc.info  <- roc(as.numeric(table[,labelcolumn]),glm_fit$fitted.values,plot = TRUE, legacy.axes=TRUE,
    percent=TRUE, xlab="FALSE Positive Percentage", ylab="True Positive Percentage",
    lwd=4)

print(coords(roc.info, "best", ret = "threshold"))
fit_values <- data.frame(glm_fit$data[,labelcolumn],glm_fit$fitted.values) 
print(fit_values)
}

