---
title: "IPAH and PH-SSc vs HV and SSc-No-PH"
output: html_notebook
  ---

```r
source('~/Google Drive File Stream//My Drive/Work/Microarrays/2014 miRNA Arrays/Final Analysis/Scripts/Custom functions.R', chdir = TRUE)
dataset <- read.csv("~/Google Drive File Stream/My Drive/miRNA/dataset.csv")

library(pacman)
p_load("kableExtra","caret","tidyverse","reshape2","e1071","JamesTools","OptimalCutpoints","Boruta","ggplot2","randomForest","ROCR","rpart","party","rpart.plot","partykit","glmnet","xgboost","Ringo", "parallel", "doParallel", "DMwR", "cowplot", "pROC", "readxl")

rownames(dataset) <- dataset$X
dataset<- select(dataset, -X)

NTproBNP <- as.data.frame(read_excel("~/Google Drive File Stream/My Drive/miRNA/NTproBNP.xlsx"))
pheno <- read_excel("~/Google Drive File Stream/My Drive/miRNA/NTproBNP.xlsx", 
    sheet = "Sheet3")
```


```r
ml.Spear <- dataset[dataset$AB == "A",] %>% select(.,-PHstatus,-AB)
ml.Spear$group <- as.factor(ml.Spear$group)
setB<- dataset[dataset$AB == "B",]  %>% select(.,-PHstatus,-AB)
setB$group <- as.factor(setB$group)
setB2 <- setB
```


# Boruta

[Boruta](https://cran.r-project.org/web/packages/Boruta/Boruta.pdf) is a 'wrapper algorithm for all relevant feature selection'.


```r
#run boruta
multi.boruta<- list()
multi<- function(i){
    fit.boruta <- Boruta(factor(group)~., data=ml.Spear, maxRuns = 300, pValue = 0.01)
boruta.df <- data.frame(attStats(fit.boruta))
multi.boruta[[i]] <- rownames(boruta.df[boruta.df$decision =='Confirmed',])
}
set.seed(100)
multi.boruta<- parallel::mclapply(1:100, multi)

#See how many times each miRNA appears when boruta is run 100x
multi.boruta<- as.data.frame(table(unlist(multi.boruta)))
multi.boruta.mirs<- as.character(multi.boruta[which(multi.boruta$Freq>10),1])
multi.boruta
```

```
##           Var1 Freq
## 1    let.7d.3p  100
## 2  miR.125a.5p    5
## 3   miR.187.5p  100
## 4   miR.18b.5p   36
## 5   miR.33b.3p    2
## 6  miR.3613.3p  100
## 7     miR.451a    2
## 8  miR.4707.5p  100
## 9      miR.572  100
## 10     miR.636  100
## 11  miR.652.3p  100
## 12  miR.671.5p   97
## 13     miR.933  100
```

## Random Forest on Boruta


```r
ml.Spear$group <- as.factor((ml.Spear$group))
m<- paste(as.vector(multi.boruta.mirs, mode = "any"), collapse = "+")
RFm<- as.formula(paste("group ~ ",m,sep = ""))
Boruta.data<- ml.Spear[,c("group",multi.boruta.mirs)]
tree.Boruta<- setB[,c("group", multi.boruta.mirs)]
```

Using caret to train:

10 fold cv, with 10 repeats.

`mtry` - number of variables randomly sampled as candidates at each split (default is sqrt(no. of variables))

`ntree` - number of trees to grow

Can't optemise `ntree` and `mtry` in the same run, so optemise `ntree` first, using default `mtry`:


```r
control<- trainControl(method = 'repeatedcv',
                       number = 10,
                       repeats = 10, classProbs = TRUE, savePredictions = TRUE)
metric <- "Accuracy"
tunegrid <- expand.grid(.mtry=seq(from =1, to =4))
modellist<- list()
for (ntree in c(100, 250, 500, 750, 1000, 1250, 1500)) {
	set.seed(100)
	fit <- caret::train(RFm, data=Boruta.data, method="rf", metric=metric, tuneGrid=tunegrid, trControl=control, ntree=ntree)
	key <- toString(ntree)
	modellist[[key]] <- fit
}
# compare results
results <- resamples(modellist)
res<- summary(results)
res<- as.data.frame(res$statistics$Accuracy)
summary(results)
```

```
## 
## Call:
## summary.resamples(object = results)
## 
## Models: 100, 250, 500, 750, 1000, 1250, 1500 
## Number of resamples: 100 
## 
## Accuracy 
##           Min.   1st Qu.    Median      Mean 3rd Qu. Max. NA's
## 100  0.5714286 0.7142857 0.8571429 0.8320238       1    1    0
## 250  0.5714286 0.7142857 0.8571429 0.8373214       1    1    0
## 500  0.5714286 0.7142857 0.8571429 0.8410714       1    1    0
## 750  0.5714286 0.7142857 0.8571429 0.8426786       1    1    0
## 1000 0.5714286 0.7142857 0.8571429 0.8443452       1    1    0
## 1250 0.5714286 0.7142857 0.8571429 0.8443452       1    1    0
## 1500 0.5714286 0.7142857 0.8571429 0.8443452       1    1    0
## 
## Kappa 
##            Min.   1st Qu.    Median      Mean 3rd Qu. Max. NA's
## 100  0.00000000 0.4086538 0.6956522 0.6376523       1    1    0
## 250  0.00000000 0.4086538 0.6956522 0.6469564       1    1    0
## 500  0.00000000 0.4086538 0.6956522 0.6548932       1    1    0
## 750  0.00000000 0.4086538 0.6956522 0.6578694       1    1    0
## 1000 0.08695652 0.4166667 0.6956522 0.6635938       1    1    0
## 1250 0.08695652 0.4166667 0.6956522 0.6635938       1    1    0
## 1500 0.08695652 0.4166667 0.6956522 0.6635938       1    1    0
```

Now optmise `mtry`:


```r
ntree = as.numeric(rownames(res)[which(res$Mean == max(res$Mean))])[1]
set.seed(100)
tunegrid <- expand.grid(.mtry=seq(from = 1, to=4, by = 0.5))
modellist<- list()
control<- trainControl(method = 'repeatedcv',
                       number = 10,
                       repeats = 10, classProbs = TRUE, savePredictions = TRUE)
	set.seed(100)
fit <- caret::train(RFm, data=Boruta.data, method="rf", metric="Accuracy", tuneGrid=tunegrid, trControl=control, ntree=ntree)
fit
```

```
## Random Forest 
## 
## 71 samples
## 10 predictors
##  2 classes: 'HV', 'PH' 
## 
## No pre-processing
## Resampling: Cross-Validated (10 fold, repeated 10 times) 
## Summary of sample sizes: 64, 63, 64, 64, 64, 65, ... 
## Resampling results across tuning parameters:
## 
##   mtry  Accuracy   Kappa    
##   1.0   0.8400595  0.6546686
##   1.5   0.8193452  0.6146054
##   2.0   0.8160714  0.6057677
##   2.5   0.8236310  0.6243189
##   3.0   0.8050000  0.5864137
##   3.5   0.7912500  0.5607784
##   4.0   0.8004762  0.5739003
## 
## Accuracy was used to select the optimal model using the largest value.
## The final value used for the model was mtry = 1.
```



```r
paste("ntree used:", ntree)
```

```
## [1] "ntree used: 1000"
```

```r
tunegrid <- expand.grid(.mtry=fit$bestTune$mtry)
control<- trainControl(method = 'repeatedcv',
                       number = 10,
                       repeats = 10,
                       savePredictions = TRUE,
                       classProbs = TRUE)
fit.Boruta <- caret::train(RFm, data=Boruta.data, method="rf", metric="Accuracy", tuneGrid=tunegrid, trControl=control, ntree=ntree)
RFBoruta.train <- predict(fit.Boruta, tree.Boruta)
BorutaInterim<- confusionMatrix(as.factor(RFBoruta.train), as.factor(tree.Boruta$group), positive = "PH")
BorutaInterim
```

```
## Confusion Matrix and Statistics
## 
##           Reference
## Prediction HV PH
##         HV 10  3
##         PH  4 19
##                                           
##                Accuracy : 0.8056          
##                  95% CI : (0.6398, 0.9181)
##     No Information Rate : 0.6111          
##     P-Value [Acc > NIR] : 0.01065         
##                                           
##                   Kappa : 0.5855          
##                                           
##  Mcnemar's Test P-Value : 1.00000         
##                                           
##             Sensitivity : 0.8636          
##             Specificity : 0.7143          
##          Pos Pred Value : 0.8261          
##          Neg Pred Value : 0.7692          
##              Prevalence : 0.6111          
##          Detection Rate : 0.5278          
##    Detection Prevalence : 0.6389          
##       Balanced Accuracy : 0.7890          
##                                           
##        'Positive' Class : PH              
## 
```

# Rpart 

The `rpart` function from the [`rpart`](https://cran.r-project.org/web/packages/rpart/rpart.pdf) package can be utilised to grow a regression tree. This tree is built by splitting the data on the single variable which best splits the group in 2. Once the data is separated, this process is applied to each sub-group separately recursively until the subgroups reach a minimum size (defined as 3 below), or no improvement can be made. 

Rpart uses a variable selection algorithm called recursive feature elimination (REF, also known as backward selection). 

`minsplit` - min no of observations that must exist in a node in order for a split to be attempted
`minbucket` - min no of observations in any terminal node


```r
		#train control
tc <- trainControl(method="repeatedcv", number=10, repeats = 10, classProbs=TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE)
#maxcompete

set.seed(20)
fit.caret.rpart <- caret::train(group ~ ., data=ml.Spear, method='rpart', metric="ROC", trControl=tc, control=rpart.control(minsplit=2, minbucket=3, surrogatestyle = 1, maxcompete = 0)) 
fit.caret.rpart$bestTune
```

```
##   cp
## 1  0
```


```r
prp(fit.caret.rpart$finalModel, main="PH from HV Rpart model", extra=2, varlen=0)
```

![](ensembleshort_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

```r
plot(as.party(fit.caret.rpart$finalModel), main="PH from HV Rpart model", drop_terminal=F)
```

![](ensembleshort_files/figure-html/unnamed-chunk-9-2.png)<!-- -->

### B validation


```r
rpartmirs<-labels(fit.caret.rpart$finalModel)[-1]
rpartmirs<- gsub("<.*","", rpartmirs)
rpartmirs<- unique(gsub(">.*","", rpartmirs))

predrpart2 <- predict(fit.caret.rpart, setB2)

fit.preds.table.Rpart<- cbind(rownames(setB2),predrpart2,as.character(setB2$group))
fit.preds.table.Rpart<- as.data.frame(fit.preds.table.Rpart)
colnames(fit.preds.table.Rpart) <- c("sample","fit.preds","group")

RpartInterim<- confusionMatrix(predrpart2, setB2$group, positive = "PH")
RpartInterim
```

```
## Confusion Matrix and Statistics
## 
##           Reference
## Prediction HV PH
##         HV  9  2
##         PH  5 20
##                                           
##                Accuracy : 0.8056          
##                  95% CI : (0.6398, 0.9181)
##     No Information Rate : 0.6111          
##     P-Value [Acc > NIR] : 0.01065         
##                                           
##                   Kappa : 0.5743          
##                                           
##  Mcnemar's Test P-Value : 0.44969         
##                                           
##             Sensitivity : 0.9091          
##             Specificity : 0.6429          
##          Pos Pred Value : 0.8000          
##          Neg Pred Value : 0.8182          
##              Prevalence : 0.6111          
##          Detection Rate : 0.5556          
##    Detection Prevalence : 0.6944          
##       Balanced Accuracy : 0.7760          
##                                           
##        'Positive' Class : PH              
## 
```


# LASSO

LASSO (Least Absolute Shrinkage and Selection Operator) is a feature selection method designed to reduce over-fitting. It automatically selects the significant variables by shrinking the coefficients of predictors deemed unimportant to zero.

`Glmnet` ([Vignette](http://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html)) is a package that uses penalised maximum likelihood to fit a generalised linear model. 


```r
fit.glmcv <- cv.glmnet(x=as.matrix(ml.Spear[-1]), y=as.factor(ml.Spear$group), alpha=1, family='binomial', nfolds=10)
#summary(fit.glmcv)
other.glmcv <- cv.glmnet(x=as.matrix(ml.Spear[-1]), y=as.factor(ml.Spear$group), alpha=1, family='binomial', nfolds=10, type.measure = "class")
model.lambda <- fit.glmcv$lambda.min

plot(fit.glmcv, cex.axis=1, cex.lab=1,cex.main=1)
```

![](ensembleshort_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

```r
plot(other.glmcv, cex.axis=1, cex.lab=1, cex.main=1)
```

![](ensembleshort_files/figure-html/unnamed-chunk-11-2.png)<!-- -->

In the plot above, the red indicates the cross validation curve, and the the error bars show upper and lower standard deviation curves. 

The two dotted lines show `lambda.min` - the value of $\lambda$ that gives minimum mean cross-validated error. The other $\lambda$ is `lambda.1se`, which gives the most regularised model such that the error is within 1 standard error of the minimum. 



```r
#refit model for new lambda
paste("lambda value:", fit.glmcv$lambda.min)
```

```
## [1] "lambda value: 0.0416631184443247"
```

```r
tuneLASSO<- expand.grid(.alpha = 1, .lambda = fit.glmcv$lambda.min)
LASSO.min<- caret::train(group ~ ., data = ml.Spear, method = "glmnet", trControl = control, family = 'binomial', tuneGrid = tuneLASSO)
lasso.model<- coef(LASSO.min$finalModel, LASSO.min$bestTune$lambda) %>% as.matrix() %>% as.data.frame() %>% rownames_to_column %>% filter(abs(`1`) >0)

ml.Spear.LASSO <- cbind("group" = ml.Spear$group, ml.Spear[,colnames(ml.Spear) %in% lasso.model$rowname])

paste("Number of miRs in LASSO model:",length(lasso.model$rowname)-1)
```

```
## [1] "Number of miRs in LASSO model: 13"
```

```r
kable(lasso.model, caption = "miRs retained by LASSO") %>% kable_styling(full_width = TRUE)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>miRs retained by LASSO</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> rowname </th>
   <th style="text-align:right;"> 1 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> (Intercept) </td>
   <td style="text-align:right;"> -15.0237033 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> let.7d.3p </td>
   <td style="text-align:right;"> 0.4949133 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.1306.3p </td>
   <td style="text-align:right;"> 0.0077301 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.148a.3p </td>
   <td style="text-align:right;"> 0.5795444 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.187.5p </td>
   <td style="text-align:right;"> 0.1561147 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.34a.5p </td>
   <td style="text-align:right;"> 2.2138448 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.451a </td>
   <td style="text-align:right;"> -0.1104036 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.4707.5p </td>
   <td style="text-align:right;"> 2.1832830 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.484 </td>
   <td style="text-align:right;"> 0.9258305 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.494 </td>
   <td style="text-align:right;"> -0.0658810 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.548am.5p </td>
   <td style="text-align:right;"> 0.2649530 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.572 </td>
   <td style="text-align:right;"> 0.4641874 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.636 </td>
   <td style="text-align:right;"> -1.3620149 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.671.5p </td>
   <td style="text-align:right;"> 0.0123772 </td>
  </tr>
</tbody>
</table>

## Validation on B


```r
val.LASSO<- setB[, colnames(setB) %in% colnames(ml.Spear)]
LASSO.pred.min <- predict(LASSO.min, newdata=val.LASSO[-1])
predicted.classes<-as.character(LASSO.pred.min)
fit.preds.table.LASSO<- as.data.frame(cbind(as.character(rownames(setB)),as.character(predicted.classes),as.character(setB2$group)))
colnames(fit.preds.table.LASSO) <- c("sample","LASSO","group")
LASSOInterim<- confusionMatrix(LASSO.pred.min,val.LASSO$group, positive = "PH")
LASSOInterim
```

```
## Confusion Matrix and Statistics
## 
##           Reference
## Prediction HV PH
##         HV  9  5
##         PH  5 17
##                                          
##                Accuracy : 0.7222         
##                  95% CI : (0.5481, 0.858)
##     No Information Rate : 0.6111         
##     P-Value [Acc > NIR] : 0.1143         
##                                          
##                   Kappa : 0.4156         
##                                          
##  Mcnemar's Test P-Value : 1.0000         
##                                          
##             Sensitivity : 0.7727         
##             Specificity : 0.6429         
##          Pos Pred Value : 0.7727         
##          Neg Pred Value : 0.6429         
##              Prevalence : 0.6111         
##          Detection Rate : 0.4722         
##    Detection Prevalence : 0.6111         
##       Balanced Accuracy : 0.7078         
##                                          
##        'Positive' Class : PH             
## 
```



# XGBoost

[XGBoost](https://cran.r-project.org/web/packages/xgboost/xgboost.pdf) (**Ex**treme **G**radient **B**oosting) is an optimised distributed gradient boosting library that performs better than gradient boosting (GBM) framework alone. 

## Parameters for tuning

`nrounds` - maximum number of iterations (similar to no. of trees grown). 

`eta` - [range: (0,1)] - learning rate. After every round, it shrinks the feature weights to reach the best optimum. Lower eta = slower computation 

`gamma` [range: $(0,\infty)$] - controls regularisation (prevents over-fitting)

`max_depth` [range: $(0,\infty)$] - controls the tree depth. The larger the depth, the more complex the model and the higher the chances of over-fitting. Larger data sets require deep trees

`min_child_weight` [range: $(0,\infty)$] - In classification, if the leaf node has a minimum sum of instance weight lower than min_child_weight, the tree splitting stops. I.e. it stops potential future interactions to reduce over-fitting

`subsample` [range: (0,1)] - Controls the number of samples supplied to a tree

`colsample_bytree` [range: (0,1)] - controls the no. of features (variables) supplied to a tree

| Parameter          | Default |Range       | Values attempted                | Final model |
|--------------------|---------|------------|---------------------------------|-------------|
|`nrounds`          | 100     |             | 100 - 10 000                    | 200 |
|`eta`              | 0.3     |(0,1)         |0.01, 0.025, 0.05, 0.1, 0.2, 0.3 | 0.025 |
|`gamma`            | 0       | $(0,\infty)$ | 0,0.05, 0.1, 0.5, 0.7, 0.9, 1   | 0.05 |
| `max_depth`       | 6       | $(0,\infty)$ | 1, 2, 3, 4, 5, 6                   | 1 |
|`colsample_bytree` | 1       | (0,1)       | 0.4, 0.6, 0.8, 1.0              | 0.6 |
|`subsample`        | 1       | (0,1)       |0.5, 0.75, 1.0                   | 0.5 |
|`min_child_weight` | 1       | $(0,\infty)$ | 1, 2, 3, 4                         | 1 |

The default metric to determine the best settings in train is accuracy, with Kappa also being calculated for a classification model.

## Default 


```r
library(xgboost)
train<- ml.Spear[-1]
train<- data.matrix(train)
validation<- setB[-1]
validation<- data.matrix(validation)

labels<- ml.Spear$group 

#Set up XGBoost with default hyperparameters
grid_default <- expand.grid(
  nrounds = 100,
  max_depth = 6,
  eta = 0.3,
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 1
)

train_control <- caret::trainControl(
  method = "none",
  verboseIter = FALSE, # no training log
  allowParallel = TRUE 
)
```


```r
xgb_base<- caret::train( 
  x = train,
  y = as.factor(labels),
  trControl = train_control,
  tuneGrid = grid_default,
  method = "xgbTree",
  verbose = TRUE,
  scale_pos_weight = (72/43) #(Total no. samples / no. positive in groupA)
  )

xgbpredict_base<- predict(xgb_base, validation)
confusionMatrix(xgbpredict_base, as.factor(setB$group), positive = "PH")
```

```
## Confusion Matrix and Statistics
## 
##           Reference
## Prediction HV PH
##         HV  9  5
##         PH  5 17
##                                          
##                Accuracy : 0.7222         
##                  95% CI : (0.5481, 0.858)
##     No Information Rate : 0.6111         
##     P-Value [Acc > NIR] : 0.1143         
##                                          
##                   Kappa : 0.4156         
##                                          
##  Mcnemar's Test P-Value : 1.0000         
##                                          
##             Sensitivity : 0.7727         
##             Specificity : 0.6429         
##          Pos Pred Value : 0.7727         
##          Neg Pred Value : 0.6429         
##              Prevalence : 0.6111         
##          Detection Rate : 0.4722         
##    Detection Prevalence : 0.6111         
##       Balanced Accuracy : 0.7078         
##                                          
##        'Positive' Class : PH             
## 
```


```r
nrounds<- seq(from = 100, to =1000, by = 50)
#tune using caret
tune_grid <- expand.grid(
  #nrounds = seq(from = 100, to = 1000, by = 20),
  nrounds = nrounds,
  eta = 0.05,
  max_depth = c(1, 2, 3),
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 1
)

set.seed(25)
tune_control <- caret::trainControl(
  method = "repeatedcv", # cross-validation
  number = 10, # with n folds
  repeats = 10, #the no. of complete sets of folds to compute
  #index = createFolds(tr_treated$Id_clean), # fix the folds
  verboseIter = FALSE, # no training log
  allowParallel = TRUE 
)

set.seed(25)
xgb_tune <- caret::train(
  x = train,
  y = as.factor(labels),
  trControl = tune_control,
  tuneGrid = tune_grid,
  method = "xgbTree",
  verbose = TRUE
  #,   scale_pos_weight = (35/21)- 1
)
kable(xgb_tune$bestTune, caption = "1st tuning, best parameters")%>% kable_styling(full_width = TRUE)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>1st tuning, best parameters</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> nrounds </th>
   <th style="text-align:right;"> max_depth </th>
   <th style="text-align:right;"> eta </th>
   <th style="text-align:right;"> gamma </th>
   <th style="text-align:right;"> colsample_bytree </th>
   <th style="text-align:right;"> min_child_weight </th>
   <th style="text-align:right;"> subsample </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 17 </td>
   <td style="text-align:right;"> 900 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.05 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
</tbody>
</table>


```r
tune_grid2 <- expand.grid(
  nrounds = nrounds,
  eta = 0.05,
  max_depth = c(xgb_tune$bestTune$max_depth -1, xgb_tune$bestTune$max_depth, xgb_tune$bestTune$max_depth +1),
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = c(0,1,2,3,4),
  subsample = 1
)

set.seed(25)
xgb_tune2 <- caret::train(
  x = train,
  y = as.factor(labels),
  trControl = tune_control,
  tuneGrid = tune_grid2,
  method = "xgbTree",
  verbose = TRUE
  #,  scale_pos_weight = (35/21) -1
)

kable(xgb_tune2$bestTune, caption = "2nd tuning, best parameters")%>% kable_styling(full_width = TRUE)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>2nd tuning, best parameters</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> nrounds </th>
   <th style="text-align:right;"> max_depth </th>
   <th style="text-align:right;"> eta </th>
   <th style="text-align:right;"> gamma </th>
   <th style="text-align:right;"> colsample_bytree </th>
   <th style="text-align:right;"> min_child_weight </th>
   <th style="text-align:right;"> subsample </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 142 </td>
   <td style="text-align:right;"> 500 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.05 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
</tbody>
</table>


```r
tune_grid3 <- expand.grid(
  nrounds = nrounds,
  eta = 0.05,
  max_depth = xgb_tune2$bestTune$max_depth,
  gamma = 0,
  colsample_bytree = c(0.4, 0.6, 0.8, 1.0),
  min_child_weight = xgb_tune2$bestTune$min_child_weight,
  subsample = c(0.5, 0.75, 1.0)
)

set.seed(25)
xgb_tune3 <- caret::train(
  x = train,
  y = as.factor(labels),
  trControl = tune_control,
  tuneGrid = tune_grid3,
  method = "xgbTree",
  verbose = TRUE
  #,   scale_pos_weight = (35/21) -1
)

kable(xgb_tune3$bestTune, caption = "3rd tuning, best parameters")%>% kable_styling(full_width = TRUE)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>3rd tuning, best parameters</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> nrounds </th>
   <th style="text-align:right;"> max_depth </th>
   <th style="text-align:right;"> eta </th>
   <th style="text-align:right;"> gamma </th>
   <th style="text-align:right;"> colsample_bytree </th>
   <th style="text-align:right;"> min_child_weight </th>
   <th style="text-align:right;"> subsample </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 12 </td>
   <td style="text-align:right;"> 650 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.05 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.4 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.5 </td>
  </tr>
</tbody>
</table>


```r
tune_grid4 <- expand.grid(
  nrounds = nrounds,
  eta = 0.05,
  max_depth = xgb_tune3$bestTune$max_depth,
  gamma = c(0,0.05, 0.1, 0.5, 0.7, 0.9, 1),
  colsample_bytree = xgb_tune3$bestTune$colsample_bytree,
  min_child_weight = xgb_tune2$bestTune$min_child_weight,
  subsample = xgb_tune3$bestTune$subsample
)

set.seed(25)
xgb_tune4 <- caret::train(
  x = train,
  y = as.factor(labels),
  trControl = tune_control,
  tuneGrid = tune_grid4,
  method = "xgbTree",
  verbose = TRUE
  #,  scale_pos_weight = (35/21) -1
)

kable(xgb_tune4$bestTune, caption = "4th tuning, best parameters")%>% kable_styling(full_width = TRUE)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>4th tuning, best parameters</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> nrounds </th>
   <th style="text-align:right;"> max_depth </th>
   <th style="text-align:right;"> eta </th>
   <th style="text-align:right;"> gamma </th>
   <th style="text-align:right;"> colsample_bytree </th>
   <th style="text-align:right;"> min_child_weight </th>
   <th style="text-align:right;"> subsample </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 28 </td>
   <td style="text-align:right;"> 500 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.05 </td>
   <td style="text-align:right;"> 0.05 </td>
   <td style="text-align:right;"> 0.4 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.5 </td>
  </tr>
</tbody>
</table>


```r
tune_grid5 <- expand.grid(
  nrounds = seq(from = 100, to = 10000, by = 50),
  eta = c(0.01,0.025,0.05,0.1),
  max_depth = xgb_tune3$bestTune$max_depth,
  gamma = xgb_tune4$bestTune$gamma,
  colsample_bytree = xgb_tune3$bestTune$colsample_bytree,
  min_child_weight = xgb_tune2$bestTune$min_child_weight,
  subsample = xgb_tune3$bestTune$subsample
)

set.seed(25)
xgb_tune5 <- caret::train(
  x = train,
  y = as.factor(labels),
  trControl = tune_control,
  tuneGrid = tune_grid5,
  method = "xgbTree",
  verbose = TRUE
  #,  scale_pos_weight = (35/21) -1
)

kable(xgb_tune5$bestTune, caption = "5th tuning, best parameters")%>% kable_styling(full_width = TRUE)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>5th tuning, best parameters</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> nrounds </th>
   <th style="text-align:right;"> max_depth </th>
   <th style="text-align:right;"> eta </th>
   <th style="text-align:right;"> gamma </th>
   <th style="text-align:right;"> colsample_bytree </th>
   <th style="text-align:right;"> min_child_weight </th>
   <th style="text-align:right;"> subsample </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 284 </td>
   <td style="text-align:right;"> 4300 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.025 </td>
   <td style="text-align:right;"> 0.05 </td>
   <td style="text-align:right;"> 0.4 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.5 </td>
  </tr>
</tbody>
</table>




```r
#look at features from final model
xgblabels <- ifelse(labels == "PH", 1,0)
xgb_final<- xgboost(data = train, label = xgblabels, nrounds = xgb_tune5$bestTune$nrounds, objective = "binary:logistic", max_depth = xgb_tune5$bestTune$max_depth, eta = xgb_tune5$bestTune$eta, min_child_weight = xgb_tune5$bestTune$min_child_weight, verbose = 0)
importance_final<- xgb.importance(feature_names = colnames(train), model = xgb_final)
xgb_mirs<- importance_final[importance_final$Gain>0.05,]

kable(xgb_mirs, caption = "Important Features (Gain > 0.05)")%>% kable_styling(full_width = TRUE)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>Important Features (Gain &gt; 0.05)</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Feature </th>
   <th style="text-align:right;"> Gain </th>
   <th style="text-align:right;"> Cover </th>
   <th style="text-align:right;"> Frequency </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> miR.636 </td>
   <td style="text-align:right;"> 0.1645373 </td>
   <td style="text-align:right;"> 0.1020921 </td>
   <td style="text-align:right;"> 0.0912052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.187.5p </td>
   <td style="text-align:right;"> 0.1396281 </td>
   <td style="text-align:right;"> 0.0728265 </td>
   <td style="text-align:right;"> 0.0464169 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> let.7d.3p </td>
   <td style="text-align:right;"> 0.1206077 </td>
   <td style="text-align:right;"> 0.1114459 </td>
   <td style="text-align:right;"> 0.0960912 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.3613.3p </td>
   <td style="text-align:right;"> 0.1008889 </td>
   <td style="text-align:right;"> 0.0771543 </td>
   <td style="text-align:right;"> 0.0602606 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.4707.5p </td>
   <td style="text-align:right;"> 0.0825629 </td>
   <td style="text-align:right;"> 0.0396843 </td>
   <td style="text-align:right;"> 0.0203583 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.18b.5p </td>
   <td style="text-align:right;"> 0.0760481 </td>
   <td style="text-align:right;"> 0.0960875 </td>
   <td style="text-align:right;"> 0.0944625 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.572 </td>
   <td style="text-align:right;"> 0.0706146 </td>
   <td style="text-align:right;"> 0.0634856 </td>
   <td style="text-align:right;"> 0.0570033 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.125a.5p </td>
   <td style="text-align:right;"> 0.0616937 </td>
   <td style="text-align:right;"> 0.0762966 </td>
   <td style="text-align:right;"> 0.0814332 </td>
  </tr>
</tbody>
</table>

```r
gg<- xgb.ggplot.importance(importance_matrix = xgb_mirs)
gg+ggplot2::theme(legend.position = "none", text = element_text(size = 22), axis.text = element_text(size = 18))
```

![](ensembleshort_files/figure-html/unnamed-chunk-21-1.png)<!-- -->

```r
xgbpredictfinal<- predict(xgb_tune5, validation)
xgbpredictfinal_probs<- predict(xgb_tune5, validation, type = "prob")
confusionMatrix(xgbpredictfinal, as.factor(setB$group), positive = "PH")
```

```
## Confusion Matrix and Statistics
## 
##           Reference
## Prediction HV PH
##         HV  9  4
##         PH  5 18
##                                          
##                Accuracy : 0.75           
##                  95% CI : (0.578, 0.8788)
##     No Information Rate : 0.6111         
##     P-Value [Acc > NIR] : 0.05907        
##                                          
##                   Kappa : 0.4671         
##                                          
##  Mcnemar's Test P-Value : 1.00000        
##                                          
##             Sensitivity : 0.8182         
##             Specificity : 0.6429         
##          Pos Pred Value : 0.7826         
##          Neg Pred Value : 0.6923         
##              Prevalence : 0.6111         
##          Detection Rate : 0.5000         
##    Detection Prevalence : 0.6389         
##       Balanced Accuracy : 0.7305         
##                                          
##        'Positive' Class : PH             
## 
```

# XGBoost on selcted miRs


```r
setAshort<- ml.Spear[,c("group",xgb_mirs$Feature)]
setBshort<- setB[,c("group",xgb_mirs$Feature)]
train.short<- data.matrix(setAshort[-1])
validation.short<- data.matrix(setBshort[-1])

labels<- setAshort$group 
setBshort$group <- setB$group

#Set up XGBoost with default hyperparameters
grid_default <- expand.grid(
  nrounds = 100,
  max_depth = 6,
  eta = 0.3,
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 1
)

train_control <- caret::trainControl(
  method = "none",
  verboseIter = FALSE, # no training log
  allowParallel = TRUE 
)

xgb_base_short<- caret::train( 
  x = train.short,
  y = as.factor(labels),
  trControl = train_control,
  tuneGrid = grid_default,
  method = "xgbTree",
  verbose = TRUE
  )

xgbpredict_base_short<- predict(xgb_base_short, validation.short)
confusionMatrix(xgbpredict_base_short, as.factor(setBshort$group), positive = "PH")
```

```
## Confusion Matrix and Statistics
## 
##           Reference
## Prediction HV PH
##         HV  7  6
##         PH  7 16
##                                           
##                Accuracy : 0.6389          
##                  95% CI : (0.4622, 0.7918)
##     No Information Rate : 0.6111          
##     P-Value [Acc > NIR] : 0.4372          
##                                           
##                   Kappa : 0.2303          
##                                           
##  Mcnemar's Test P-Value : 1.0000          
##                                           
##             Sensitivity : 0.7273          
##             Specificity : 0.5000          
##          Pos Pred Value : 0.6957          
##          Neg Pred Value : 0.5385          
##              Prevalence : 0.6111          
##          Detection Rate : 0.4444          
##    Detection Prevalence : 0.6389          
##       Balanced Accuracy : 0.6136          
##                                           
##        'Positive' Class : PH              
## 
```

## First Tune

Check model accuracy for different tree depths, for nrounds between 100 and 1000


```r
nrounds<- seq(from = 100, to =1000, by = 50)
#tune using caret
tune_grid <- expand.grid(
  #nrounds = seq(from = 100, to = 1000, by = 50),
  nrounds = nrounds,
  eta = c(0.025, 0.05, 0.1, 0.2, 0.3),
  max_depth = c(1, 2, 3, 4, 5, 6),
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 1
)

set.seed(25)
tune_control <- caret::trainControl(
  method = "repeatedcv", # cross-validation
  number = 10, # with n folds
  repeats = 10, #the no. of complete sets of folds to compute
  #index = createFolds(tr_treated$Id_clean), # fix the folds
  verboseIter = FALSE, # no training log
  allowParallel = TRUE,
  savePredictions = TRUE
)

set.seed(25)
xgb_tune_short <- caret::train(
  x = train.short,
  y = as.factor(labels),
  trControl = tune_control,
  tuneGrid = tune_grid,
  method = "xgbTree",
  verbose = TRUE
)

kable(xgb_tune_short$bestTune, caption = "1st tuning, best parameters")%>% kable_styling(full_width = TRUE)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>1st tuning, best parameters</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> nrounds </th>
   <th style="text-align:right;"> max_depth </th>
   <th style="text-align:right;"> eta </th>
   <th style="text-align:right;"> gamma </th>
   <th style="text-align:right;"> colsample_bytree </th>
   <th style="text-align:right;"> min_child_weight </th>
   <th style="text-align:right;"> subsample </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 10 </td>
   <td style="text-align:right;"> 550 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.025 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
</tbody>
</table>

## Second Tune - Maximum Depth and Minimum Child Weight


```r
tune_grid2_short <- expand.grid(
  nrounds = nrounds,
  eta = xgb_tune_short$bestTune$eta,
  max_depth = c(xgb_tune_short$bestTune$max_depth -1, xgb_tune_short$bestTune$max_depth, xgb_tune_short$bestTune$max_depth +1),
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = c(1,2,3,4),
  subsample = 1
)

set.seed(25)
xgb_tune2_short <- caret::train(
  x = train.short,
  y = as.factor(labels),
  trControl = tune_control,
  tuneGrid = tune_grid2_short,
  method = "xgbTree",
  verbose = TRUE
)

kable(xgb_tune2_short$bestTune, caption = "2nd tuning, best parameters")%>% kable_styling(full_width = TRUE)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>2nd tuning, best parameters</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> nrounds </th>
   <th style="text-align:right;"> max_depth </th>
   <th style="text-align:right;"> eta </th>
   <th style="text-align:right;"> gamma </th>
   <th style="text-align:right;"> colsample_bytree </th>
   <th style="text-align:right;"> min_child_weight </th>
   <th style="text-align:right;"> subsample </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 86 </td>
   <td style="text-align:right;"> 550 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.025 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
</tbody>
</table>

Minimum Sum of Instance Weight refers to the `min_child_weight` parameter

## Third tune - Column and Row sampling


```r
tune_grid3_short <- expand.grid(
  nrounds = nrounds,
  eta = xgb_tune_short$bestTune$eta,
  max_depth = xgb_tune2_short$bestTune$max_depth,
  gamma = 0,
  colsample_bytree = c(0.4, 0.6, 0.8, 1.0),
  min_child_weight = xgb_tune2_short$bestTune$min_child_weight,
  subsample = c(0.5, 0.75, 1.0)
)

set.seed(25)
xgb_tune3_short <- caret::train(
  x = train.short,
  y = as.factor(labels),
  trControl = tune_control,
  tuneGrid = tune_grid3_short,
  method = "xgbTree",
  verbose = TRUE
)

kable(xgb_tune3_short$bestTune, caption = "3rd tuning, best parameters")%>% kable_styling(full_width = TRUE)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>3rd tuning, best parameters</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> nrounds </th>
   <th style="text-align:right;"> max_depth </th>
   <th style="text-align:right;"> eta </th>
   <th style="text-align:right;"> gamma </th>
   <th style="text-align:right;"> colsample_bytree </th>
   <th style="text-align:right;"> min_child_weight </th>
   <th style="text-align:right;"> subsample </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 13 </td>
   <td style="text-align:right;"> 700 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.025 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.4 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.5 </td>
  </tr>
</tbody>
</table>

## Fourth tune - gamma


```r
tune_grid4_short <- expand.grid(
  nrounds = nrounds,
  eta = xgb_tune_short$bestTune$eta,
  max_depth = xgb_tune3_short$bestTune$max_depth,
  gamma = c(0,0.05, 0.1, 0.5, 0.7, 0.9, 1),
  colsample_bytree = xgb_tune3_short$bestTune$colsample_bytree,
  min_child_weight = xgb_tune2_short$bestTune$min_child_weight,
  subsample = xgb_tune3_short$bestTune$subsample
)

set.seed(25)
xgb_tune4_short <- caret::train(
  x = train.short,
  y = as.factor(labels),
  trControl = tune_control,
  tuneGrid = tune_grid4_short,
  method = "xgbTree",
  verbose = TRUE
)

kable(xgb_tune4_short$bestTune, caption = "4th tuning, best parameters")%>% kable_styling(full_width = TRUE)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>4th tuning, best parameters</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> nrounds </th>
   <th style="text-align:right;"> max_depth </th>
   <th style="text-align:right;"> eta </th>
   <th style="text-align:right;"> gamma </th>
   <th style="text-align:right;"> colsample_bytree </th>
   <th style="text-align:right;"> min_child_weight </th>
   <th style="text-align:right;"> subsample </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 119 </td>
   <td style="text-align:right;"> 300 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.025 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.4 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.5 </td>
  </tr>
</tbody>
</table>

## Fifth tune - Reducing Learning Rate


```r
tune_grid5_short <- expand.grid(
  nrounds = seq(from = 100, to = 10000, by = 50),
  eta = c(0.01,0.025,0.05,0.1),
  max_depth = xgb_tune3_short$bestTune$max_depth,
  gamma = xgb_tune4_short$bestTune$gamma,
  colsample_bytree = xgb_tune3_short$bestTune$colsample_bytree,
  min_child_weight = xgb_tune2_short$bestTune$min_child_weight,
  subsample = xgb_tune3_short$bestTune$subsample
)

set.seed(25)
xgb_tune5_short <- caret::train(
  x = train.short,
  y = as.factor(labels),
  trControl = tune_control,
  tuneGrid = tune_grid5_short,
  method = "xgbTree",
  verbose = TRUE
)

kable(xgb_tune5_short$bestTune, caption = "5th tuning, best parameters")%>% kable_styling(full_width = TRUE)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>5th tuning, best parameters</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> nrounds </th>
   <th style="text-align:right;"> max_depth </th>
   <th style="text-align:right;"> eta </th>
   <th style="text-align:right;"> gamma </th>
   <th style="text-align:right;"> colsample_bytree </th>
   <th style="text-align:right;"> min_child_weight </th>
   <th style="text-align:right;"> subsample </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 202 </td>
   <td style="text-align:right;"> 200 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.025 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.4 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.5 </td>
  </tr>
</tbody>
</table>

## Final model


```r
final_grid_short <- expand.grid(
  nrounds = xgb_tune5_short$bestTune$nrounds,
  eta = xgb_tune5_short$bestTune$eta,
  max_depth = xgb_tune5_short$bestTune$max_depth,
  gamma = xgb_tune5_short$bestTune$gamma,
  colsample_bytree = xgb_tune5_short$bestTune$colsample_bytree,
  min_child_weight = xgb_tune5_short$bestTune$min_child_weight,
  subsample = xgb_tune5_short$bestTune$subsample
)

set.seed(25)
tune_control_f <- caret::trainControl(
  method = "repeatedcv", # cross-validation
  number = 10, # with n folds
  repeats = 10, #the no. of complete sets of folds to compute
  #index = createFolds(tr_treated$Id_clean), # fix the folds
  verboseIter = FALSE, # no training log
  allowParallel = TRUE,
  savePredictions = TRUE,
  classProbs = TRUE
)
set.seed(25)
xgb_tune_final_short <- caret::train(
  x = train.short,
  y = as.factor(labels),
  trControl = tune_control_f,
  tuneGrid = final_grid_short,
  method = "xgbTree",
  verbose = TRUE
)

xgbpredictfinalshort<- predict(xgb_tune_final_short, validation.short)
XGBInterim<- confusionMatrix(xgbpredictfinalshort, as.factor(setBshort$group), positive = "PH")
XGBInterim
```

```
## Confusion Matrix and Statistics
## 
##           Reference
## Prediction HV PH
##         HV 10  2
##         PH  4 20
##                                           
##                Accuracy : 0.8333          
##                  95% CI : (0.6719, 0.9363)
##     No Information Rate : 0.6111          
##     P-Value [Acc > NIR] : 0.003604        
##                                           
##                   Kappa : 0.64            
##                                           
##  Mcnemar's Test P-Value : 0.683091        
##                                           
##             Sensitivity : 0.9091          
##             Specificity : 0.7143          
##          Pos Pred Value : 0.8333          
##          Neg Pred Value : 0.8333          
##              Prevalence : 0.6111          
##          Detection Rate : 0.5556          
##    Detection Prevalence : 0.6667          
##       Balanced Accuracy : 0.8117          
##                                           
##        'Positive' Class : PH              
## 
```


# Univariate analysis 


```r
#mir in any feature selection method
all.mirs<- unique(c(lasso.model$rowname[-1],multi.boruta.mirs, as.character(rpartmirs), xgb_mirs$Feature))
#miRs in at least 2 selection methods
mirnames<- duplicated(c(lasso.model$rowname[-1],multi.boruta.mirs, rpartmirs, xgb_mirs$Feature))
mirnames<- c(lasso.model$rowname[-1],multi.boruta.mirs, rpartmirs, xgb_mirs$Feature)[mirnames] %>% unique(.)
```

miRNA in any feature selection method: 

`all.mirs`

miRNAs in at least 2 feature selection methods:

`mirnames`

## ROC

ROC curves for all microRNAs selected by either LASSO, Rpart, Boruta or XGBoost. NB calculated across all HV and PAH , not traning set


```r
#ROC curves for IPAH
library(pROC)
roc.res.i <- list()
     for(i in all.mirs){
          #Save ROC data
         roc.res.i[[i]] <- roc(dataset$group, dataset[,i], plot=F)
         plot.roc(roc.res.i[[i]], lwd = 5, cex.axis = 1, cex.lab = 1, cex.main = 1, main=paste("ROC for",i))
     }
```

![](ensembleshort_files/figure-html/unnamed-chunk-30-1.png)<!-- -->![](ensembleshort_files/figure-html/unnamed-chunk-30-2.png)<!-- -->![](ensembleshort_files/figure-html/unnamed-chunk-30-3.png)<!-- -->![](ensembleshort_files/figure-html/unnamed-chunk-30-4.png)<!-- -->![](ensembleshort_files/figure-html/unnamed-chunk-30-5.png)<!-- -->![](ensembleshort_files/figure-html/unnamed-chunk-30-6.png)<!-- -->![](ensembleshort_files/figure-html/unnamed-chunk-30-7.png)<!-- -->![](ensembleshort_files/figure-html/unnamed-chunk-30-8.png)<!-- -->![](ensembleshort_files/figure-html/unnamed-chunk-30-9.png)<!-- -->![](ensembleshort_files/figure-html/unnamed-chunk-30-10.png)<!-- -->![](ensembleshort_files/figure-html/unnamed-chunk-30-11.png)<!-- -->![](ensembleshort_files/figure-html/unnamed-chunk-30-12.png)<!-- -->![](ensembleshort_files/figure-html/unnamed-chunk-30-13.png)<!-- -->![](ensembleshort_files/figure-html/unnamed-chunk-30-14.png)<!-- -->![](ensembleshort_files/figure-html/unnamed-chunk-30-15.png)<!-- -->![](ensembleshort_files/figure-html/unnamed-chunk-30-16.png)<!-- -->![](ensembleshort_files/figure-html/unnamed-chunk-30-17.png)<!-- -->![](ensembleshort_files/figure-html/unnamed-chunk-30-18.png)<!-- -->![](ensembleshort_files/figure-html/unnamed-chunk-30-19.png)<!-- -->![](ensembleshort_files/figure-html/unnamed-chunk-30-20.png)<!-- -->


```r
cutoffs <- find.Cutpoints.Loop(all.mirs, data=ml.Spear[-1], pheno=ml.Spear[,1], healthytag="PH", method='MaxSpSe')
kable(cutoffs, caption = "Cutpoints for miRNAs selected by the above methods, calculated on the training set") %>% kable_styling(full_width = TRUE)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>Cutpoints for miRNAs selected by the above methods, calculated on the training set</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> mir </th>
   <th style="text-align:right;"> cutpoint </th>
   <th style="text-align:left;"> direction </th>
   <th style="text-align:right;"> sensitivity </th>
   <th style="text-align:right;"> specificity </th>
   <th style="text-align:left;"> auc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> let.7d.3p </td>
   <td style="text-align:right;"> 3.085032 </td>
   <td style="text-align:left;"> &gt; </td>
   <td style="text-align:right;"> 0.6551724 </td>
   <td style="text-align:right;"> 0.6428571 </td>
   <td style="text-align:left;"> 0.663 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.1306.3p </td>
   <td style="text-align:right;"> 2.788380 </td>
   <td style="text-align:left;"> &gt; </td>
   <td style="text-align:right;"> 0.6206897 </td>
   <td style="text-align:right;"> 0.6190476 </td>
   <td style="text-align:left;"> 0.65 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.148a.3p </td>
   <td style="text-align:right;"> 2.686263 </td>
   <td style="text-align:left;"> &gt; </td>
   <td style="text-align:right;"> 0.6551724 </td>
   <td style="text-align:right;"> 0.6428571 </td>
   <td style="text-align:left;"> 0.668 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.187.5p </td>
   <td style="text-align:right;"> 2.555543 </td>
   <td style="text-align:left;"> &gt; </td>
   <td style="text-align:right;"> 0.7241379 </td>
   <td style="text-align:right;"> 0.7142857 </td>
   <td style="text-align:left;"> 0.782 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.34a.5p </td>
   <td style="text-align:right;"> 2.792594 </td>
   <td style="text-align:left;"> &gt; </td>
   <td style="text-align:right;"> 0.5862069 </td>
   <td style="text-align:right;"> 0.5952381 </td>
   <td style="text-align:left;"> 0.621 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.451a </td>
   <td style="text-align:right;"> 11.518547 </td>
   <td style="text-align:left;"> &lt; </td>
   <td style="text-align:right;"> 0.6206897 </td>
   <td style="text-align:right;"> 0.6428571 </td>
   <td style="text-align:left;"> 0.687 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.4707.5p </td>
   <td style="text-align:right;"> 2.642439 </td>
   <td style="text-align:left;"> &gt; </td>
   <td style="text-align:right;"> 0.6551724 </td>
   <td style="text-align:right;"> 0.6428571 </td>
   <td style="text-align:left;"> 0.66 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.484 </td>
   <td style="text-align:right;"> 3.057025 </td>
   <td style="text-align:left;"> &gt; </td>
   <td style="text-align:right;"> 0.5517241 </td>
   <td style="text-align:right;"> 0.5476190 </td>
   <td style="text-align:left;"> 0.584 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.494 </td>
   <td style="text-align:right;"> 3.870084 </td>
   <td style="text-align:left;"> &lt; </td>
   <td style="text-align:right;"> 0.6206897 </td>
   <td style="text-align:right;"> 0.5238095 </td>
   <td style="text-align:left;"> 0.571 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.548am.5p </td>
   <td style="text-align:right;"> 2.880562 </td>
   <td style="text-align:left;"> &gt; </td>
   <td style="text-align:right;"> 0.5517241 </td>
   <td style="text-align:right;"> 0.5714286 </td>
   <td style="text-align:left;"> 0.541 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.572 </td>
   <td style="text-align:right;"> 4.507390 </td>
   <td style="text-align:left;"> &gt; </td>
   <td style="text-align:right;"> 0.6206897 </td>
   <td style="text-align:right;"> 0.5952381 </td>
   <td style="text-align:left;"> 0.678 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.636 </td>
   <td style="text-align:right;"> 3.408426 </td>
   <td style="text-align:left;"> &lt; </td>
   <td style="text-align:right;"> 0.7241379 </td>
   <td style="text-align:right;"> 0.7380952 </td>
   <td style="text-align:left;"> 0.8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.671.5p </td>
   <td style="text-align:right;"> 3.092697 </td>
   <td style="text-align:left;"> &gt; </td>
   <td style="text-align:right;"> 0.6206897 </td>
   <td style="text-align:right;"> 0.6190476 </td>
   <td style="text-align:left;"> 0.67 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.18b.5p </td>
   <td style="text-align:right;"> 2.627057 </td>
   <td style="text-align:left;"> &lt; </td>
   <td style="text-align:right;"> 0.6896552 </td>
   <td style="text-align:right;"> 0.6190476 </td>
   <td style="text-align:left;"> 0.567 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.3613.3p </td>
   <td style="text-align:right;"> 3.705675 </td>
   <td style="text-align:left;"> &lt; </td>
   <td style="text-align:right;"> 0.8275862 </td>
   <td style="text-align:right;"> 0.6666667 </td>
   <td style="text-align:left;"> 0.746 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.652.3p </td>
   <td style="text-align:right;"> 2.537779 </td>
   <td style="text-align:left;"> &gt; </td>
   <td style="text-align:right;"> 0.6896552 </td>
   <td style="text-align:right;"> 0.6904762 </td>
   <td style="text-align:left;"> 0.706 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.933 </td>
   <td style="text-align:right;"> 2.987387 </td>
   <td style="text-align:left;"> &lt; </td>
   <td style="text-align:right;"> 0.7586207 </td>
   <td style="text-align:right;"> 0.6904762 </td>
   <td style="text-align:left;"> 0.737 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.151a.3p </td>
   <td style="text-align:right;"> 2.639110 </td>
   <td style="text-align:left;"> &gt; </td>
   <td style="text-align:right;"> 0.6551724 </td>
   <td style="text-align:right;"> 0.6428571 </td>
   <td style="text-align:left;"> 0.667 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.1246 </td>
   <td style="text-align:right;"> 3.706243 </td>
   <td style="text-align:left;"> &gt; </td>
   <td style="text-align:right;"> 0.5172414 </td>
   <td style="text-align:right;"> 0.5000000 </td>
   <td style="text-align:left;"> 0.522 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.125a.5p </td>
   <td style="text-align:right;"> 2.687269 </td>
   <td style="text-align:left;"> &gt; </td>
   <td style="text-align:right;"> 0.6896552 </td>
   <td style="text-align:right;"> 0.6904762 </td>
   <td style="text-align:left;"> 0.679 </td>
  </tr>
</tbody>
</table>

miR-187-5p


```r
m187cut <- cutoffs %>% filter(mir == "miR.187.5p")
m187tab<- table(ifelse(setB$miR.187.5p > m187cut$cutpoint, "PH", "HV"), setB$group)
(m187tab[1]+m187tab[4])/sum(m187tab)
```

```
## [1] 0.7777778
```

miR-636


```r
m636cut <- cutoffs %>% filter(mir == "miR.636")
m636tab<- table(ifelse(setB$miR.636 < m636cut$cutpoint, "PH","HV"), setB$group)
(m636tab[1]+m636tab[4])/sum(m636tab)
```

```
## [1] 0.6944444
```

## Wilcox Tests


```r
compared<- dataset[,c("group", all.mirs)]

sigs<- list()
  for(i in 2:(length(compared))) {
    sigs[[i]] <- wilcox.test(
      compared[compared$group == "HV",i],
      compared[compared$group == "PH",i], alternative = "two.sided"
    )
  }
  names(sigs)<- colnames(compared)
sigs <-  sigs[-1]

df<- data.frame(mir.name=character(),p.value=double() ,stringsAsFactors = FALSE)

for(i in 1:length(sigs)) {
  unlist(sigs[i])
 df[i,] <- c(names(sigs[i]),sigs[[i]]$p.value)
}

df$p.adj.BH<- p.adjust(df$p.value, method = "BH")
rownames(df)<- df$mir.name

kable(df[-1], caption = "Wilcox tests, adjusted p-value BH method, miRs selcted by feature selection") %>%  kable_styling(full_width = TRUE)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>Wilcox tests, adjusted p-value BH method, miRs selcted by feature selection</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> p.value </th>
   <th style="text-align:right;"> p.adj.BH </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> let.7d.3p </td>
   <td style="text-align:left;"> 0.421515639091875 </td>
   <td style="text-align:right;"> 0.4437007 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.1306.3p </td>
   <td style="text-align:left;"> 0.0170342695582462 </td>
   <td style="text-align:right;"> 0.0283904 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.148a.3p </td>
   <td style="text-align:left;"> 0.0242887381217238 </td>
   <td style="text-align:right;"> 0.0373673 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.187.5p </td>
   <td style="text-align:left;"> 9.59281571025942e-08 </td>
   <td style="text-align:right;"> 0.0000019 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.34a.5p </td>
   <td style="text-align:left;"> 0.115801851075843 </td>
   <td style="text-align:right;"> 0.1415116 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.451a </td>
   <td style="text-align:left;"> 0.00635005507654273 </td>
   <td style="text-align:right;"> 0.0127001 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.4707.5p </td>
   <td style="text-align:left;"> 0.0022645387494186 </td>
   <td style="text-align:right;"> 0.0057824 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.484 </td>
   <td style="text-align:left;"> 0.0545909875121561 </td>
   <td style="text-align:right;"> 0.0779871 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.494 </td>
   <td style="text-align:left;"> 0.351917353704318 </td>
   <td style="text-align:right;"> 0.3910193 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.548am.5p </td>
   <td style="text-align:left;"> 0.946806105179864 </td>
   <td style="text-align:right;"> 0.9468061 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.572 </td>
   <td style="text-align:left;"> 0.00231295892090655 </td>
   <td style="text-align:right;"> 0.0057824 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.636 </td>
   <td style="text-align:left;"> 2.16919244720858e-06 </td>
   <td style="text-align:right;"> 0.0000217 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.671.5p </td>
   <td style="text-align:left;"> 0.000261797038495134 </td>
   <td style="text-align:right;"> 0.0013090 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.18b.5p </td>
   <td style="text-align:left;"> 0.0928287163659173 </td>
   <td style="text-align:right;"> 0.1237716 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.3613.3p </td>
   <td style="text-align:left;"> 2.1595648726422e-05 </td>
   <td style="text-align:right;"> 0.0001440 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.652.3p </td>
   <td style="text-align:left;"> 0.00178991207196812 </td>
   <td style="text-align:right;"> 0.0057824 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.933 </td>
   <td style="text-align:left;"> 0.00164133147224172 </td>
   <td style="text-align:right;"> 0.0057824 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.151a.3p </td>
   <td style="text-align:left;"> 0.00599323424402427 </td>
   <td style="text-align:right;"> 0.0127001 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.1246 </td>
   <td style="text-align:left;"> 0.120284827480817 </td>
   <td style="text-align:right;"> 0.1415116 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.125a.5p </td>
   <td style="text-align:left;"> 0.0135653887117804 </td>
   <td style="text-align:right;"> 0.0246643 </td>
  </tr>
</tbody>
</table>

# Model MicroRNAs

Summary statistics for PH group:


```r
#pick mirs selected by any of the 3 methods
selectedmirs <- c("group",all.mirs)

#reduce dataframe to contain only mirs in selectedmirs
selecmirs<- ml.Spear[, colnames(ml.Spear) %in% selectedmirs]
#melt dataframe so ggplot can be used
df<- melt(selecmirs, id.var = "group") %>% mutate (group = factor(group, ))

#print summaries of PH and HV groups
H<- selecmirs[selecmirs$group == "HV",]
PH<- selecmirs[selecmirs$group == "PH",]
summary(PH)
```

```
##  group     let.7d.3p        miR.1246      miR.125a.5p     miR.1306.3p   
##  HV: 0   Min.   :2.726   Min.   :2.595   Min.   :2.608   Min.   :2.585  
##  PH:42   1st Qu.:3.035   1st Qu.:3.079   1st Qu.:2.675   1st Qu.:2.727  
##          Median :3.231   Median :3.704   Median :2.712   Median :2.858  
##          Mean   :3.440   Mean   :3.938   Mean   :2.879   Mean   :2.985  
##          3rd Qu.:3.521   3rd Qu.:4.564   3rd Qu.:2.895   3rd Qu.:3.085  
##          Max.   :7.135   Max.   :6.123   Max.   :5.046   Max.   :4.329  
##   miR.148a.3p     miR.151a.3p      miR.187.5p      miR.18b.5p   
##  Min.   :2.581   Min.   :2.577   Min.   :2.473   Min.   :2.576  
##  1st Qu.:2.662   1st Qu.:2.632   1st Qu.:2.545   1st Qu.:2.601  
##  Median :2.710   Median :2.653   Median :2.695   Median :2.618  
##  Mean   :2.831   Mean   :2.810   Mean   :2.829   Mean   :2.635  
##  3rd Qu.:2.903   3rd Qu.:2.829   3rd Qu.:2.833   3rd Qu.:2.644  
##  Max.   :3.771   Max.   :4.406   Max.   :5.340   Max.   :2.787  
##    miR.34a.5p     miR.3613.3p       miR.451a       miR.4707.5p   
##  Min.   :2.682   Min.   :2.882   Min.   : 6.925   Min.   :2.499  
##  1st Qu.:2.754   1st Qu.:3.102   1st Qu.: 9.349   1st Qu.:2.571  
##  Median :2.808   Median :3.517   Median :10.379   Median :2.737  
##  Mean   :2.895   Mean   :3.828   Mean   :10.509   Mean   :2.869  
##  3rd Qu.:2.951   3rd Qu.:4.082   3rd Qu.:11.861   3rd Qu.:3.073  
##  Max.   :3.478   Max.   :8.708   Max.   :14.237   Max.   :3.648  
##     miR.484         miR.494       miR.548am.5p      miR.572     
##  Min.   :2.707   Min.   :2.642   Min.   :2.688   Min.   :2.541  
##  1st Qu.:2.995   1st Qu.:2.945   1st Qu.:2.807   1st Qu.:4.199  
##  Median :3.079   Median :3.777   Median :2.884   Median :4.743  
##  Mean   :3.113   Mean   :3.993   Mean   :2.997   Mean   :4.673  
##  3rd Qu.:3.185   3rd Qu.:4.544   3rd Qu.:3.038   3rd Qu.:5.231  
##  Max.   :3.821   Max.   :8.154   Max.   :4.141   Max.   :7.474  
##     miR.636        miR.652.3p      miR.671.5p       miR.933     
##  Min.   :2.721   Min.   :2.498   Min.   :2.492   Min.   :2.601  
##  1st Qu.:3.006   1st Qu.:2.533   1st Qu.:2.919   1st Qu.:2.832  
##  Median :3.221   Median :2.564   Median :3.399   Median :2.938  
##  Mean   :3.271   Mean   :2.629   Mean   :4.557   Mean   :2.966  
##  3rd Qu.:3.407   3rd Qu.:2.631   3rd Qu.:6.167   3rd Qu.:3.016  
##  Max.   :4.954   Max.   :3.596   Max.   :9.917   Max.   :3.670
```

Summary statistics for healthy volunteers:


```r
summary(H)
```

```
##  group     let.7d.3p        miR.1246       miR.125a.5p     miR.1306.3p   
##  HV:29   Min.   :2.752   Min.   : 2.649   Min.   :2.626   Min.   :2.581  
##  PH: 0   1st Qu.:2.964   1st Qu.: 3.060   1st Qu.:2.668   1st Qu.:2.676  
##          Median :3.058   Median : 3.706   Median :2.678   Median :2.758  
##          Mean   :3.132   Mean   : 4.067   Mean   :2.714   Mean   :2.806  
##          3rd Qu.:3.136   3rd Qu.: 4.475   3rd Qu.:2.702   3rd Qu.:2.825  
##          Max.   :4.537   Max.   :13.529   Max.   :3.657   Max.   :3.545  
##   miR.148a.3p     miR.151a.3p      miR.187.5p      miR.18b.5p   
##  Min.   :2.612   Min.   :2.593   Min.   :2.453   Min.   :2.586  
##  1st Qu.:2.650   1st Qu.:2.620   1st Qu.:2.496   1st Qu.:2.607  
##  Median :2.681   Median :2.628   Median :2.534   Median :2.632  
##  Mean   :2.691   Mean   :2.657   Mean   :2.558   Mean   :2.631  
##  3rd Qu.:2.696   3rd Qu.:2.649   3rd Qu.:2.576   3rd Qu.:2.651  
##  Max.   :2.998   Max.   :3.215   Max.   :3.024   Max.   :2.684  
##    miR.34a.5p     miR.3613.3p       miR.451a       miR.4707.5p   
##  Min.   :2.649   Min.   :3.455   Min.   : 7.365   Min.   :2.482  
##  1st Qu.:2.751   1st Qu.:3.752   1st Qu.:10.870   1st Qu.:2.583  
##  Median :2.791   Median :4.042   Median :12.053   Median :2.631  
##  Mean   :2.798   Mean   :4.394   Mean   :11.645   Mean   :2.638  
##  3rd Qu.:2.835   3rd Qu.:4.970   3rd Qu.:12.649   3rd Qu.:2.680  
##  Max.   :3.087   Max.   :6.831   Max.   :14.218   Max.   :2.836  
##     miR.484         miR.494       miR.548am.5p      miR.572     
##  Min.   :2.765   Min.   :2.637   Min.   :2.701   Min.   :2.620  
##  1st Qu.:2.917   1st Qu.:3.421   1st Qu.:2.818   1st Qu.:3.798  
##  Median :3.047   Median :3.942   Median :2.865   Median :4.372  
##  Mean   :3.067   Mean   :4.426   Mean   :2.908   Mean   :4.203  
##  3rd Qu.:3.131   3rd Qu.:5.366   3rd Qu.:2.933   3rd Qu.:4.681  
##  Max.   :3.849   Max.   :8.353   Max.   :3.563   Max.   :5.936  
##     miR.636        miR.652.3p      miR.671.5p       miR.933     
##  Min.   :2.737   Min.   :2.506   Min.   :2.538   Min.   :2.796  
##  1st Qu.:3.349   1st Qu.:2.523   1st Qu.:2.638   1st Qu.:2.987  
##  Median :3.825   Median :2.534   Median :2.852   Median :3.056  
##  Mean   :3.950   Mean   :2.543   Mean   :3.699   Mean   :3.143  
##  3rd Qu.:4.382   3rd Qu.:2.542   3rd Qu.:3.756   3rd Qu.:3.245  
##  Max.   :5.509   Max.   :2.737   Max.   :9.023   Max.   :3.800
```

```r
ggplot(data = df, aes(x= variable, y=value)) + geom_boxplot(aes(fill=group)) + theme(text = element_text(size=12), axis.text.x=element_text(angle=90, vjust=0.5, size = 12),axis.text.y=element_text(size = 12)) + labs(x="", y="Expression", size = 12)
```

![](ensembleshort_files/figure-html/unnamed-chunk-36-1.png)<!-- -->



```r
#melt dataframe into SSc, SSc-no-ph, IPAH, HV
selectedmirs <- c("PHstatus",all.mirs)

#reduce dataset to contain only mirs in selectedmirs - all subjects not just set A
selecmirs2<- dataset[,which(colnames(dataset) %in% selectedmirs)]
dfSSc<- melt(selecmirs2, id.var = "PHstatus")
dfSSc$PHstatus<- factor(dfSSc$PHstatus, levels = c("HV","No_PH_CTD","CTD_PAH","IPAH"))
dfSScPH<- dfSSc[dfSSc$PHstatus == 'CTD_PAH' | dfSSc$PHstatus == 'IPAH', ]
ggplot(data = dfSScPH, aes(x= variable, y=value)) + geom_boxplot(aes(fill=PHstatus)) + theme(axis.text.x=element_text(angle=90, vjust=0.5, size = 12),axis.text.y=element_text(size = 12)) + labs(x="", y="Expression")
```

![](ensembleshort_files/figure-html/unnamed-chunk-37-1.png)<!-- -->


```r
ggplot(data = dfSSc, aes(x= variable, y=value)) + geom_boxplot(aes(fill=PHstatus)) + theme(axis.text.x=element_text(angle=90, vjust=0.5, size = 12),axis.text.y=element_text(size = 12)) + labs(x="", y="Expression")
```

![](ensembleshort_files/figure-html/unnamed-chunk-38-1.png)<!-- -->


## miRs in more than one feature selection method


```r
#melt dataframe into SSc, SSc-no-ph, IPAH, HV
mormirs <- c("group",unique(mirnames))

#reduce datset to contain only mirs in selectedmirs
moremirs<- dataset[,mormirs]
moremirs<- melt(moremirs, id.var = "group")

ggplot(data = moremirs, aes(x= variable, y=value)) + geom_boxplot(aes(fill=group)) + theme(axis.text.x=element_text(angle=90, vjust=0.5, size = 12),axis.text.y=element_text(size = 12),axis.title.y=element_text(size=12)) + labs(x="", y="Expression") 
```

![](ensembleshort_files/figure-html/unnamed-chunk-39-1.png)<!-- -->

```r
moremirsAB<- dataset[,which(colnames(dataset) %in% c(unique(mirnames), "AB"))]
AB<- moremirsAB %>% melt(id.var = "AB")

ggplot(data = AB, aes(x= variable, y=value)) + geom_boxplot(aes(fill=AB)) + theme(axis.text.x=element_text(angle=90, vjust=0.5, size = 12),axis.text.y=element_text(size = 12)) + labs(x="", y="Expression")
```

![](ensembleshort_files/figure-html/unnamed-chunk-39-2.png)<!-- -->

```r
final2<- moremirs[moremirs$variable =="miR-636" | moremirs$variable =="miR-187-5p" ,]
#ggplot(data = final3, aes(x= variable, y=value)) + geom_boxplot(aes(fill=group),size=1.2) + theme(axis.text.x=element_text(angle=60, vjust=0.5, size = 16),axis.text.y=element_text(size = 12),axis.title.y=element_text(size=12), legend.text = element_text(size = 12)) + labs(x="", y="Expression")
```

# CV ROCs

## Boruta


```r
boruta.split<- split(fit.Boruta$pred, fit.Boruta$pred$Resample)
parallel::mclapply(1:100, function(d) {
  pROC::auc(pROC::roc(predictor = boruta.split[[d]]$PH, response = boruta.split[[d]]$obs))[1]
}) %>% unlist() %>% hist(main = paste("Boruta CV mean AUC:", round(mean(.),2)))
```

![](ensembleshort_files/figure-html/unnamed-chunk-40-1.png)<!-- -->

## Rpart


```r
rpart.auc<- list()
rpart.split<- split(fit.caret.rpart$pred, fit.caret.rpart$pred$Resample)
parallel::mclapply(1:100, function(d) {
  pROC::auc(pROC::roc(predictor = rpart.split[[d]]$PH, response = rpart.split[[d]]$obs))[1]
}) %>% unlist() %>% hist(main = paste("Rpart CV mean AUC:", round(mean(.),2)))
```

![](ensembleshort_files/figure-html/unnamed-chunk-41-1.png)<!-- -->

## LASSO (caret model - lambda min)


```r
LASSO.split.min<- split(LASSO.min$pred, LASSO.min$pred$Resample)
parallel::mclapply(1:100, function(d) {
  pROC::auc(pROC::roc(predictor = LASSO.split.min[[d]]$PH, response = LASSO.split.min[[d]]$obs))[1]
}) %>% unlist() %>% hist(main = paste("LASSO CV mean AUC:", round(mean(.),2)))
```

![](ensembleshort_files/figure-html/unnamed-chunk-42-1.png)<!-- -->

## XGBoost


```r
tosplit<- xgb_tune_final_short
XGB.split<- split(tosplit$pred, tosplit$pred$Resample)
parallel::mclapply(1:100, function(d) {
  pROC::auc(pROC::roc(predictor = XGB.split[[d]]$PH, response = XGB.split[[d]]$obs))[1]
}) %>% unlist() %>% hist(main = paste("XGBoost CV mean AUC:", round(mean(.),2)))
```

![](ensembleshort_files/figure-html/unnamed-chunk-43-1.png)<!-- -->

# ROC on interim set

## Boruta


```r
Boruta.perf<- predict(fit.Boruta, tree.Boruta[-1], type= "prob")
Boruta.pred<- prediction(Boruta.perf$PH, tree.Boruta$group)
Boruta.perfs<- ROCR::performance(Boruta.pred,"tpr","fpr")
Boruta.sens<- ROCR::performance(Boruta.pred,"sens","spec")
Boruta.auc <- pROC::auc(tree.Boruta$group, Boruta.perf[,2])
```

```
## Setting levels: control = HV, case = PH
```

```
## Setting direction: controls < cases
```

```r
Boruta.aucCI<- pROC::ci.auc(tree.Boruta$group, Boruta.perf[,2])
```

```
## Setting levels: control = HV, case = PH
## Setting direction: controls < cases
```

```r
par(mfrow=c(1,2))
plot(Boruta.perfs, main = paste("AUC:", round(Boruta.auc,2)))
plot(Boruta.sens)
```

![](ensembleshort_files/figure-html/unnamed-chunk-44-1.png)<!-- -->

```r
Boruta.aucCI
```

```
## 95% CI: 0.6902-0.9981 (DeLong)
```

## Rpart


```r
rpart.perf<- predict(fit.caret.rpart, setB2[-1], type= "prob")
rpart.pred<- prediction(rpart.perf$PH, setB2$group)
rpart.perfs<- ROCR::performance(rpart.pred,"tpr","fpr")
rpart.sens<- ROCR::performance(rpart.pred,"sens","spec")
rpart.auc<- pROC::auc(setB2$group, rpart.perf[,2])
```

```
## Setting levels: control = HV, case = PH
```

```
## Setting direction: controls < cases
```

```r
rpart.aucCI<- pROC::ci.auc(setB2$group, rpart.perf[,2])
```

```
## Setting levels: control = HV, case = PH
## Setting direction: controls < cases
```

```r
par(mfrow=c(1,2))
plot(rpart.perfs, main = paste("AUC:", round(rpart.auc,2)))
plot(rpart.sens)
```

![](ensembleshort_files/figure-html/unnamed-chunk-45-1.png)<!-- -->

```r
rpart.aucCI
```

```
## 95% CI: 0.6321-0.9458 (DeLong)
```

## LASSO 


```r
LASSO.probs.min <- predict(LASSO.min, newdata=val.LASSO[-1], type = "prob")
LASSO.pred.min<- prediction(LASSO.probs.min$PH, val.LASSO$group)
LASSO.perfs.min<- ROCR::performance(LASSO.pred.min,"tpr","fpr")
LASSO.sens.min<- ROCR::performance(LASSO.pred.min,"sens","spec")
LASSO.auc.min<- pROC::auc(val.LASSO$group, LASSO.probs.min[,2])
```

```
## Setting levels: control = HV, case = PH
```

```
## Setting direction: controls < cases
```

```r
LASSO.aucCI<- pROC::ci.auc(val.LASSO$group, LASSO.probs.min[,2])
```

```
## Setting levels: control = HV, case = PH
## Setting direction: controls < cases
```

```r
par(mfrow=c(1,2))
plot(LASSO.perfs.min, main = paste("AUC:", round(LASSO.auc.min,2)))
plot(LASSO.sens.min)
```

![](ensembleshort_files/figure-html/unnamed-chunk-46-1.png)<!-- -->

```r
LASSO.aucCI
```

```
## 95% CI: 0.6343-0.9436 (DeLong)
```

## XGBoost


```r
xgbpreds<- predict(xgb_tune_final_short, validation.short, type = "prob")
xgb.pred<- prediction(xgbpreds$PH, as.factor(setBshort$group))
xgb.perfs<- ROCR::performance(xgb.pred,"tpr","fpr")
xgb.sens<- ROCR::performance(xgb.pred,"sens","spec")
xgb.auc<- pROC::auc(as.factor(setBshort$group), xgbpreds[,2])
```

```
## Setting levels: control = HV, case = PH
```

```
## Setting direction: controls < cases
```

```r
xgb.aucCI<- pROC::ci.auc(as.factor(setBshort$group), xgbpreds[,2])
```

```
## Setting levels: control = HV, case = PH
## Setting direction: controls < cases
```

```r
par(mfrow=c(1,2))
plot(xgb.perfs, main = paste("AUC:", round(xgb.auc,2)))
plot(xgb.sens)
```

![](ensembleshort_files/figure-html/unnamed-chunk-47-1.png)<!-- -->

```r
xgb.aucCI
```

```
## 95% CI: 0.6576-0.9852 (DeLong)
```


```r
plot(xgb.perfs, col = "blue", lwd = 2.5, cex = 2)
plot(LASSO.perfs.min, col = "light blue", add = TRUE, lwd = 2.5)
plot(rpart.perfs, col = "dark blue", add = TRUE, lwd = 2.5)
plot(Boruta.perfs, col = "turquoise", add = TRUE, lwd = 2.5)
abline(0,1, lty= 2, lwd = 2)
legend("bottomright", legend = c("XGBoost", "LASSO", "Rpart", "Random Forest"), col = c("blue", "light blue", "dark blue", "turquoise"), lty = 1, inset = 0.1, lwd = 2.5)
```

![](ensembleshort_files/figure-html/unnamed-chunk-48-1.png)<!-- -->

# Selected miRs


```r
kable(multi.boruta.mirs, col.names = paste("Boruta miRs: (",length(multi.boruta.mirs),")"))
```

<table>
 <thead>
  <tr>
   <th style="text-align:left;"> Boruta miRs: ( 10 ) </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> let.7d.3p </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.187.5p </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.18b.5p </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.3613.3p </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.4707.5p </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.572 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.636 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.652.3p </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.671.5p </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.933 </td>
  </tr>
</tbody>
</table>

```r
kable(as.character(rpartmirs), col.names = paste("Rpart miRs: (",length(rpartmirs),")"))
```

<table>
 <thead>
  <tr>
   <th style="text-align:left;"> Rpart miRs: ( 4 ) </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> miR.187.5p </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.636 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.151a.3p </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.1246 </td>
  </tr>
</tbody>
</table>

```r
kable(lasso.model$rowname[-1], col.names = paste("LASSO miRs: (", length(lasso.model$rowname[-1]),")"))
```

<table>
 <thead>
  <tr>
   <th style="text-align:left;"> LASSO miRs: ( 13 ) </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> let.7d.3p </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.1306.3p </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.148a.3p </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.187.5p </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.34a.5p </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.451a </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.4707.5p </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.484 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.494 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.548am.5p </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.572 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.636 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.671.5p </td>
  </tr>
</tbody>
</table>

```r
kable(xgb_mirs$Feature, col.names = paste("XGBoost miRs: (",length(xgb_mirs$Feature),")"))
```

<table>
 <thead>
  <tr>
   <th style="text-align:left;"> XGBoost miRs: ( 8 ) </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> miR.636 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.187.5p </td>
  </tr>
  <tr>
   <td style="text-align:left;"> let.7d.3p </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.3613.3p </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.4707.5p </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.18b.5p </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.572 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miR.125a.5p </td>
  </tr>
</tbody>
</table>

```r
paste("Unique miRs:", length(unique(c(multi.boruta.mirs, rpartmirs, lasso.model$rowname[-1], xgb_mirs$Feature))))
```

```
## [1] "Unique miRs: 20"
```


```r
library(VennDiagram)
v1<- venn.diagram(x=list(A= as.vector(lasso.model$rowname),B=as.vector(rpartmirs),C=as.vector(multi.boruta.mirs), D=as.vector(xgb_mirs$Feature)), category = c("LASSO miRs","Rpart miRs","Random forest miRs","XGBoost miRs"), fill = c("skyblue","pink1","mediumorchid","dark blue"), lty = "blank", fontfamily = "sans", cat.fontfamily = "sans", filename=NULL, simplify=TRUE, cex = 1.3, cat.cex = 1.3)
grid.newpage()
grid.draw(v1)
```

![](ensembleshort_files/figure-html/unnamed-chunk-50-1.png)<!-- -->

# Ensemble approach


```r
totalprobs<- data.frame(patientID = rownames(Boruta.perf), group = setB$group, RandomForest = Boruta.perf$PH, rpart = rpart.perf$PH, LASSO = LASSO.probs.min$PH, XGB = xgbpreds$PH)
totalprobs <- mutate(totalprobs, Mean = rowMeans(totalprobs[3:6]))
totalprobs$EnsemblePreds <- ifelse(totalprobs$Mean > 0.5, "PH", "HV")
Ensemble.pred<-prediction(totalprobs$Mean, totalprobs$group) 
Ensemble.perfs<- ROCR::performance(Ensemble.pred, "tpr", "fpr")
kable(totalprobs)  %>% kable_styling(full_width = TRUE)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> patientID </th>
   <th style="text-align:left;"> group </th>
   <th style="text-align:right;"> RandomForest </th>
   <th style="text-align:right;"> rpart </th>
   <th style="text-align:right;"> LASSO </th>
   <th style="text-align:right;"> XGB </th>
   <th style="text-align:right;"> Mean </th>
   <th style="text-align:left;"> EnsemblePreds </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 773_1_3.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.861 </td>
   <td style="text-align:right;"> 0.8666667 </td>
   <td style="text-align:right;"> 0.8610238 </td>
   <td style="text-align:right;"> 0.9266345 </td>
   <td style="text-align:right;"> 0.8788313 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 988_2_1.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.625 </td>
   <td style="text-align:right;"> 0.8666667 </td>
   <td style="text-align:right;"> 0.6270985 </td>
   <td style="text-align:right;"> 0.8060075 </td>
   <td style="text-align:right;"> 0.7311932 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 987_2_2.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.867 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.9681849 </td>
   <td style="text-align:right;"> 0.7909521 </td>
   <td style="text-align:right;"> 0.8873035 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 988_2_2.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.642 </td>
   <td style="text-align:right;"> 0.8666667 </td>
   <td style="text-align:right;"> 0.7828652 </td>
   <td style="text-align:right;"> 0.8272272 </td>
   <td style="text-align:right;"> 0.7796898 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 988_2_3.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.469 </td>
   <td style="text-align:right;"> 0.8666667 </td>
   <td style="text-align:right;"> 0.6646667 </td>
   <td style="text-align:right;"> 0.4601494 </td>
   <td style="text-align:right;"> 0.6151207 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 991_1_1.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.892 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.9760982 </td>
   <td style="text-align:right;"> 0.9143730 </td>
   <td style="text-align:right;"> 0.9263870 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 994_2_4.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.112 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.4433199 </td>
   <td style="text-align:right;"> 0.1719547 </td>
   <td style="text-align:right;"> 0.1931823 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 756_1_4.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.827 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.9729316 </td>
   <td style="text-align:right;"> 0.8753925 </td>
   <td style="text-align:right;"> 0.8996003 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 780_2_1.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.588 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.5069344 </td>
   <td style="text-align:right;"> 0.3914426 </td>
   <td style="text-align:right;"> 0.6023635 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 756_1_1.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.805 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.8110336 </td>
   <td style="text-align:right;"> 0.8562466 </td>
   <td style="text-align:right;"> 0.8488393 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 990_1_4.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.908 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.8860309 </td>
   <td style="text-align:right;"> 0.9222433 </td>
   <td style="text-align:right;"> 0.9098378 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 990_1_3.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.671 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.4652240 </td>
   <td style="text-align:right;"> 0.7577598 </td>
   <td style="text-align:right;"> 0.7042652 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 759_2_4.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.918 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.7153382 </td>
   <td style="text-align:right;"> 0.7805382 </td>
   <td style="text-align:right;"> 0.8342383 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 760_1_1.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.888 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.9594336 </td>
   <td style="text-align:right;"> 0.7051297 </td>
   <td style="text-align:right;"> 0.8689100 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 990_1_1.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.799 </td>
   <td style="text-align:right;"> 0.8666667 </td>
   <td style="text-align:right;"> 0.7554185 </td>
   <td style="text-align:right;"> 0.8007538 </td>
   <td style="text-align:right;"> 0.8054598 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 756_2_2.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.621 </td>
   <td style="text-align:right;"> 0.8666667 </td>
   <td style="text-align:right;"> 0.7759468 </td>
   <td style="text-align:right;"> 0.6305977 </td>
   <td style="text-align:right;"> 0.7235528 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 991_2_3.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.858 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.8211853 </td>
   <td style="text-align:right;"> 0.8018932 </td>
   <td style="text-align:right;"> 0.8510389 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 989_2_2.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.713 </td>
   <td style="text-align:right;"> 0.8666667 </td>
   <td style="text-align:right;"> 0.7855566 </td>
   <td style="text-align:right;"> 0.7922395 </td>
   <td style="text-align:right;"> 0.7893657 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 760_2_1.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.710 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.6660852 </td>
   <td style="text-align:right;"> 0.7215491 </td>
   <td style="text-align:right;"> 0.5357722 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 770_1_2.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.877 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.8527485 </td>
   <td style="text-align:right;"> 0.9214918 </td>
   <td style="text-align:right;"> 0.8935793 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 779_1_2.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.713 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.5425231 </td>
   <td style="text-align:right;"> 0.7151947 </td>
   <td style="text-align:right;"> 0.7234487 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 779_1_4.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.631 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.2968497 </td>
   <td style="text-align:right;"> 0.5210212 </td>
   <td style="text-align:right;"> 0.5929869 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 986_1_4.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.540 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.3144224 </td>
   <td style="text-align:right;"> 0.5690590 </td>
   <td style="text-align:right;"> 0.3672340 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 994_1_3.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.240 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.0996630 </td>
   <td style="text-align:right;"> 0.2303311 </td>
   <td style="text-align:right;"> 0.1538622 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 994_1_4.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.251 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.1631164 </td>
   <td style="text-align:right;"> 0.2609853 </td>
   <td style="text-align:right;"> 0.1801391 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 980_1_1.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.327 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.2625634 </td>
   <td style="text-align:right;"> 0.2860165 </td>
   <td style="text-align:right;"> 0.2302586 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 980_1_2.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.251 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.3797412 </td>
   <td style="text-align:right;"> 0.4075463 </td>
   <td style="text-align:right;"> 0.2709355 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 982_1_3.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.489 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.4489525 </td>
   <td style="text-align:right;"> 0.4375283 </td>
   <td style="text-align:right;"> 0.3552338 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 758_1_2.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.298 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.3030106 </td>
   <td style="text-align:right;"> 0.2726333 </td>
   <td style="text-align:right;"> 0.2297746 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 982_1_1.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.420 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.5840516 </td>
   <td style="text-align:right;"> 0.4051780 </td>
   <td style="text-align:right;"> 0.3636710 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 780_1_1.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.396 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.4244571 </td>
   <td style="text-align:right;"> 0.4439179 </td>
   <td style="text-align:right;"> 0.5468630 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 773_2_1.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.942 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.9280066 </td>
   <td style="text-align:right;"> 0.8514320 </td>
   <td style="text-align:right;"> 0.9111289 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 767_2_3.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.492 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.2876439 </td>
   <td style="text-align:right;"> 0.5135296 </td>
   <td style="text-align:right;"> 0.3346570 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 770_2_2.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.893 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.9282741 </td>
   <td style="text-align:right;"> 0.8327190 </td>
   <td style="text-align:right;"> 0.8942675 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 767_2_4.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.255 </td>
   <td style="text-align:right;"> 0.8000000 </td>
   <td style="text-align:right;"> 0.2361473 </td>
   <td style="text-align:right;"> 0.2724249 </td>
   <td style="text-align:right;"> 0.3908930 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 774_1_3.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.401 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.3690218 </td>
   <td style="text-align:right;"> 0.5977038 </td>
   <td style="text-align:right;"> 0.5727006 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
</tbody>
</table>


```r
pROC::auc(totalprobs$group, totalprobs$Mean)
```

```
## Setting levels: control = HV, case = PH
```

```
## Setting direction: controls < cases
```

```
## Area under the curve: 0.8474
```

```r
pROC::ci.auc(totalprobs$group, totalprobs$Mean)
```

```
## Setting levels: control = HV, case = PH
## Setting direction: controls < cases
```

```
## 95% CI: 0.695-0.9998 (DeLong)
```

```r
table(totalprobs$EnsemblePreds, totalprobs$group)
```

```
##     
##      HV PH
##   HV  9  2
##   PH  5 20
```


# Variable Importance 


```r
RFvars<-varImp(fit.Boruta)[[1]]%>% rownames_to_column() %>% select(miR = rowname, RFIMP = Overall)
rpartvars<- varImp(fit.caret.rpart)[[1]] %>% filter(Overall > 0) %>% rownames_to_column() %>% select(miR = rowname, rpartIMP = Overall)
LASSOvars<- varImp(LASSO.min)[[1]]  %>% filter(Overall > 0)%>% rownames_to_column() %>% select(miR = rowname, LASSOIMP = Overall)
XGBvars<- varImp(xgb_tune_final_short)[[1]]%>% rownames_to_column() %>% select(miR = rowname, XGBIMP = Overall)
n<- max(dim(RFvars)[1], dim(rpartvars)[1], dim(LASSOvars)[1], dim(XGBvars)[1])
allimps <- dplyr::full_join(x = RFvars, y=  rpartvars) %>% full_join(x=., y=LASSOvars) %>% full_join(x=., XGBvars)
```

```
## Joining, by = "miR"
## Joining, by = "miR"
## Joining, by = "miR"
```

```r
plot(varImp(fit.Boruta))
```

![](ensembleshort_files/figure-html/unnamed-chunk-53-1.png)<!-- -->

```r
plot(varImp(fit.caret.rpart))
```

![](ensembleshort_files/figure-html/unnamed-chunk-53-2.png)<!-- -->

```r
plot(varImp(LASSO.min))
```

![](ensembleshort_files/figure-html/unnamed-chunk-53-3.png)<!-- -->

```r
plot(varImp(xgb_tune_final_short))
```

![](ensembleshort_files/figure-html/unnamed-chunk-53-4.png)<!-- -->


```r
reshaped<- melt(allimps, id.vars = 1)
ggplot(reshaped, aes(x=miR, y = value)) + geom_bar(aes(fill=variable),stat="identity",position="dodge")  + theme(axis.text.x=element_text(angle=90, vjust=0.5, size = 12), panel.background = element_blank(), panel.grid = element_line(colour="grey")) + labs(x="", y="Importance", size=12) + scale_fill_brewer(palette = "RdYlBu") 
```

```
## Warning: Removed 45 rows containing missing values (geom_bar).
```

![](ensembleshort_files/figure-html/unnamed-chunk-54-1.png)<!-- -->

```r
ggplot(reshaped, aes(x=miR, y = value)) + geom_bar(aes(fill=variable),stat="identity")  + theme(axis.text.x=element_text(angle=90, vjust=0.5, size = 12), panel.background = element_blank(), panel.grid = element_line(colour="grey")) + labs(x="", y="Importance", size=12) +  scale_fill_manual(name = "", labels = c("Random Forest", "rpart", "LASSO", "XGBoost"), values = c("dark red", "tomato2", "midnight blue", "cornflowerblue"))
```

```
## Warning: Removed 45 rows containing missing values (position_stack).
```

![](ensembleshort_files/figure-html/unnamed-chunk-54-2.png)<!-- -->

# Mean centered data (Figures 2 and 4)

Figure includes training and validation set.


```r
#melt dataframe into SSc, SSc-no-ph, IPAH, HV
mormirs <- c("group",mirnames)

#reduce dataset to contain only mirs in selectedmirs
moremirs<- dataset[,mormirs]

#Mean centre data:
centre_colmeans <- function(x) {
    xcentre = colMeans(x)
    x - rep(xcentre, rep.int(nrow(x), ncol(x)))
}

meanc<- centre_colmeans(moremirs[-1])
meanc$group <- moremirs$group
meanc.melt<- melt(meanc, id.var = "group")

ggplot(data = meanc.melt, aes(x= variable, y=value)) + geom_boxplot(aes(fill=group)) + theme(axis.text.x=element_text(angle=90, vjust=0.5, size = 12),axis.text.y=element_text(size = 12),axis.title.y=element_text(size=12), panel.background = element_rect(fill = 'white'), axis.line = element_line(colour = 'black')) + labs(x="", y="Mean centered expression") + scale_colour_manual(values = c("Dark Blue", "Light Blue")) + scale_fill_manual(values = c("Midnight blue", "cornflowerblue")) 
```

![](ensembleshort_files/figure-html/unnamed-chunk-55-1.png)<!-- -->


```r
mirheat<- data.frame(Boruta = all.mirs, Rpart = all.mirs, LASSO = all.mirs, XGBoost = all.mirs, stringsAsFactors = FALSE)
rownames(mirheat)<- all.mirs

mirheat$Boruta <- ifelse(mirheat$Boruta %in% multi.boruta.mirs, "miR selected", "miR not selected")
mirheat$Rpart<- ifelse(mirheat$Rpart %in% as.character(rpartmirs),"miR selected", "miR not selected")
mirheat$LASSO <- ifelse(mirheat$LASSO %in% lasso.model$rowname[-1], "miR selected", "miR not selected")
mirheat$XGBoost <- ifelse(mirheat$XGBoost %in% xgb_mirs$Feature,"miR selected", "miR not selected")

library(pheatmap)

colors <- colorRampPalette(brewer.pal(9, "Blues"))(255)
my_colour = list(
    XGBoost = c(`miR selected` = "midnightblue", `miR not selected`= "cornsilk1"),
    LASSO = c(`miR selected` = "midnightblue", `miR not selected`= "cornsilk1"),
    Rpart = c(`miR selected` = "midnightblue", `miR not selected`= "cornsilk1"),
    Boruta = c(`miR selected` = "midnightblue", `miR not selected`= "cornsilk1")
)

forheatmap<- dataset[,colnames(dataset) %in% all.mirs]
sampleDists<- cor(forheatmap, method = "spearman")
sampleDists<- abs(sampleDists)
sampleDistsM<- as.matrix(sampleDists)
pheatmap(sampleDistsM, annotation_row = mirheat, col = colors, annotation_colors = my_colour, clustering_method = "average")
```

![](ensembleshort_files/figure-html/unnamed-chunk-56-1.png)<!-- -->


# With NTproBNP

Classification of patients based on:

Normal <125 pg/ml in under 75's

Normal <450 pg/ml in over 75's


```r
library(readxl)
 NTproBNP <- as.data.frame(read_excel("~/Google Drive File Stream/My Drive/miRNA/NTproBNP.xlsx"))
rownames(NTproBNP)<- NTproBNP$filename
table(NTproBNP[,6:7])
```

```
##       Predicted
## Status HV PH
##     HV 21  4
##     PH  2 42
```


```r
dataset$filename <- rownames(dataset)
withBNP<- dplyr::left_join(dataset, NTproBNP, by = "filename")
rownames(withBNP)<- withBNP$filename
withBNP<- dplyr::select(withBNP, -biobankid, -uid, -filename, -Status, -Predicted, -PHstatus, -Normal.range) 
withBNPa<- withBNP[withBNP$AB == "A",] %>% select(-AB) %>% filter(NTproBNP > 0)
withBNPb<- withBNP[withBNP$AB == "B",]%>% select(-AB) %>% filter(NTproBNP > 0)
```


# Logistic regression NTproBNP


```r
withBNPa$group <- as.factor(withBNPa$group)
proBNP<- glm(group ~ NTproBNP, data = withBNPa, family = "binomial")
summary(proBNP)
```

```
## 
## Call:
## glm(formula = group ~ NTproBNP, family = "binomial", data = withBNPa)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -2.8404  -0.6935   0.0031   0.3460   1.7358  
## 
## Coefficients:
##              Estimate Std. Error z value Pr(>|z|)  
## (Intercept) -1.489534   0.605324  -2.461   0.0139 *
## NTproBNP     0.004907   0.001966   2.496   0.0126 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 58.704  on 43  degrees of freedom
## Residual deviance: 31.976  on 42  degrees of freedom
## AIC: 35.976
## 
## Number of Fisher Scoring iterations: 8
```

```r
withBNPb$group<- as.factor(withBNPb$group)
NTproBNPprob<- predict(proBNP, newdata = withBNPb, type = "response")
NTproBNPpred<- rep("HV", dim(withBNPb)[1])
NTproBNPpred[NTproBNPprob > 0.5] = "PH"
NTproBNPpred<- as.factor(NTproBNPpred)
confusionMatrix(NTproBNPpred,withBNPb$group, positive = "PH")
```

```
## Confusion Matrix and Statistics
## 
##           Reference
## Prediction HV PH
##         HV  8  4
##         PH  0 12
##                                           
##                Accuracy : 0.8333          
##                  95% CI : (0.6262, 0.9526)
##     No Information Rate : 0.6667          
##     P-Value [Acc > NIR] : 0.05935         
##                                           
##                   Kappa : 0.6667          
##                                           
##  Mcnemar's Test P-Value : 0.13361         
##                                           
##             Sensitivity : 0.7500          
##             Specificity : 1.0000          
##          Pos Pred Value : 1.0000          
##          Neg Pred Value : 0.6667          
##              Prevalence : 0.6667          
##          Detection Rate : 0.5000          
##    Detection Prevalence : 0.5000          
##       Balanced Accuracy : 0.8750          
##                                           
##        'Positive' Class : PH              
## 
```

```r
pred.BNP<- prediction(NTproBNPprob, as.factor(withBNPb$group), label.ordering = c("HV","PH"))
perfs.BNP<- ROCR::performance(pred.BNP,"tpr","fpr")
sens.BNP<- ROCR::performance(pred.BNP,"sens","spec")
auc.BNP<- pROC::auc(as.factor(withBNPb$group), NTproBNPprob)
aucCI.BNP<- pROC::ci.auc(as.factor(withBNPb$group), NTproBNPprob)
aucCI.BNP
```

```
## 95% CI: 0.8404-1 (DeLong)
```

```r
par(mfrow=c(1,2))
plot(perfs.BNP, main = paste("AUC:", round(auc.BNP,2)))
plot(sens.BNP)
```

![](ensembleshort_files/figure-html/unnamed-chunk-59-1.png)<!-- -->


# Re run models with NTproBNP

Add NTproBNP to existing models, using continuous NTproBNP

## Boruta


```r
m2<- paste(m, "+", "NTproBNP")
RFm2<- as.formula(paste("group ~ ",m2,sep = ""))
Boruta.data2<- withBNPa[,c("group",multi.boruta.mirs, "NTproBNP")]
tree.Boruta2<- withBNPb[,c("group", multi.boruta.mirs, "NTproBNP")]
```


```r
tunegrid <- expand.grid(.mtry=1)
control<- trainControl(method = 'repeatedcv',
                       number = 10,
                       repeats = 10,
                       savePredictions = TRUE,
                       classProbs = TRUE)
fit.Boruta2 <- caret::train(RFm2, data=Boruta.data2, method="rf", metric="Accuracy", tuneGrid=tunegrid, trControl=control, ntree=ntree)
RFBoruta.train2 <- predict(fit.Boruta2, tree.Boruta2)
confusionMatrix(as.factor(RFBoruta.train2), tree.Boruta2$group, positive = "PH")
```

```
## Confusion Matrix and Statistics
## 
##           Reference
## Prediction HV PH
##         HV  6  0
##         PH  2 16
##                                         
##                Accuracy : 0.9167        
##                  95% CI : (0.73, 0.9897)
##     No Information Rate : 0.6667        
##     P-Value [Acc > NIR] : 0.004871      
##                                         
##                   Kappa : 0.8           
##                                         
##  Mcnemar's Test P-Value : 0.479500      
##                                         
##             Sensitivity : 1.0000        
##             Specificity : 0.7500        
##          Pos Pred Value : 0.8889        
##          Neg Pred Value : 1.0000        
##              Prevalence : 0.6667        
##          Detection Rate : 0.6667        
##    Detection Prevalence : 0.7500        
##       Balanced Accuracy : 0.8750        
##                                         
##        'Positive' Class : PH            
## 
```

## Rpart



```r
rpartBNPdata<- withBNPa[,colnames(withBNPa) %in% c("group", rpartmirs, "NTproBNP")]
tc <- trainControl(method="repeatedcv", number=10, repeats = 10, classProbs=TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE)

	set.seed(100)
fit.caret.rpart2 <- caret::train(group ~ ., data=rpartBNPdata, method='rpart', metric="ROC", trControl=tc, control=rpart.control(minsplit=2, minbucket=3, surrogatestyle = 1, maxcompete = 0)) 
fit.caret.rpart2$bestTune
```

```
##          cp
## 2 0.4117647
```


```r
prp(fit.caret.rpart2$finalModel, main="PH from HV Rpart model", extra=2, varlen=0)
```

![](ensembleshort_files/figure-html/unnamed-chunk-63-1.png)<!-- -->

```r
plot(as.party(fit.caret.rpart2$finalModel), main="PH from HV Rpart model", drop_terminal=F)
```

![](ensembleshort_files/figure-html/unnamed-chunk-63-2.png)<!-- -->


```r
predrpart2 <- predict(fit.caret.rpart2, withBNPb)
predrpart2<- as.character(predrpart2)
confusionMatrix(as.factor(predrpart2), withBNPb$group, positive = "PH")
```

```
## Confusion Matrix and Statistics
## 
##           Reference
## Prediction HV PH
##         HV  6  2
##         PH  2 14
##                                           
##                Accuracy : 0.8333          
##                  95% CI : (0.6262, 0.9526)
##     No Information Rate : 0.6667          
##     P-Value [Acc > NIR] : 0.05935         
##                                           
##                   Kappa : 0.625           
##                                           
##  Mcnemar's Test P-Value : 1.00000         
##                                           
##             Sensitivity : 0.8750          
##             Specificity : 0.7500          
##          Pos Pred Value : 0.8750          
##          Neg Pred Value : 0.7500          
##              Prevalence : 0.6667          
##          Detection Rate : 0.5833          
##    Detection Prevalence : 0.6667          
##       Balanced Accuracy : 0.8125          
##                                           
##        'Positive' Class : PH              
## 
```

## LASSO


```r
newLASSO <- withBNPa[, colnames(withBNPa) %in% c(colnames(ml.Spear.LASSO), "NTproBNP")]
tuneLASSO<- expand.grid(.alpha = 1, .lambda = fit.glmcv$lambda.min)
LASSO.min2<- caret::train(group ~ ., data = newLASSO, method = "glmnet", trControl = control, family = 'binomial', tuneGrid = tuneLASSO)
lasso.model2<- coef(LASSO.min2$finalModel, LASSO.min2$bestTune$lambda) %>% as.matrix() %>% as.data.frame() %>% rownames_to_column %>% filter(abs(`1`) >0)
```


```r
val.LASSO2<- withBNPb[, colnames(withBNPb) %in% colnames(newLASSO)]
LASSO.pred.min2 <- predict(LASSO.min2, newdata=val.LASSO2[-1])
predicted.classes<-as.character(LASSO.pred.min2)
confusionMatrix(LASSO.pred.min2,val.LASSO2$group, positive = "PH")
```

```
## Confusion Matrix and Statistics
## 
##           Reference
## Prediction HV PH
##         HV  5  3
##         PH  3 13
##                                           
##                Accuracy : 0.75            
##                  95% CI : (0.5329, 0.9023)
##     No Information Rate : 0.6667          
##     P-Value [Acc > NIR] : 0.2632          
##                                           
##                   Kappa : 0.4375          
##                                           
##  Mcnemar's Test P-Value : 1.0000          
##                                           
##             Sensitivity : 0.8125          
##             Specificity : 0.6250          
##          Pos Pred Value : 0.8125          
##          Neg Pred Value : 0.6250          
##              Prevalence : 0.6667          
##          Detection Rate : 0.5417          
##    Detection Prevalence : 0.6667          
##       Balanced Accuracy : 0.7188          
##                                           
##        'Positive' Class : PH              
## 
```



## XGBoost



```r
train2<- withBNPa[, colnames(withBNPa) %in% c(colnames(train.short), "NTproBNP")]
train2<- data.matrix(train2)
validation2<- withBNPb[, colnames(withBNPb) %in% c(colnames(train.short), "NTproBNP")]
validation2<- data.matrix(validation2)

labels<- withBNPa$group
labelsB <- withBNPb$group

set.seed(25)
xgb_tune_final_short2 <- caret::train(
  x = train2,
  y = as.factor(labels),
  trControl = tune_control,
  tuneGrid = final_grid_short,
  method = "xgbTree",
  verbose = TRUE
)

xgbpredictfinalshort2<- predict(xgb_tune_final_short2, validation2)
confusionMatrix(xgbpredictfinalshort2, as.factor(labelsB), positive = "PH")
```

```
## Confusion Matrix and Statistics
## 
##           Reference
## Prediction HV PH
##         HV  6  1
##         PH  2 15
##                                           
##                Accuracy : 0.875           
##                  95% CI : (0.6764, 0.9734)
##     No Information Rate : 0.6667          
##     P-Value [Acc > NIR] : 0.0199          
##                                           
##                   Kappa : 0.7097          
##                                           
##  Mcnemar's Test P-Value : 1.0000          
##                                           
##             Sensitivity : 0.9375          
##             Specificity : 0.7500          
##          Pos Pred Value : 0.8824          
##          Neg Pred Value : 0.8571          
##              Prevalence : 0.6667          
##          Detection Rate : 0.6250          
##    Detection Prevalence : 0.7083          
##       Balanced Accuracy : 0.8438          
##                                           
##        'Positive' Class : PH              
## 
```


# ROC on interim set for models including NTproBNP

## Boruta


```r
Boruta.perf2<- predict(fit.Boruta2, tree.Boruta2[-1], type= "prob")
Boruta.pred2<- prediction(Boruta.perf2$PH, tree.Boruta2$group)
Boruta.perfs2<- ROCR::performance(Boruta.pred2,"tpr","fpr")
Boruta.sens2<- ROCR::performance(Boruta.pred2,"sens","spec")
Boruta.auc2<- pROC::auc(tree.Boruta2$group, Boruta.perf2[,2])
Boruta.aucCI2<- pROC::ci.auc(tree.Boruta2$group, Boruta.perf2[,2])
par(mfrow=c(1,2))
plot(Boruta.perfs2, main = paste("AUC:", round(Boruta.auc2,2)))
plot(Boruta.sens2)
```

![](ensembleshort_files/figure-html/unnamed-chunk-68-1.png)<!-- -->

```r
Boruta.aucCI2
```

```
## 95% CI: 0.9211-1 (DeLong)
```

## Rpart


```r
rpart.perf2<- predict(fit.caret.rpart2, withBNPb[-1], type= "prob")
rpart.pred2<- prediction(rpart.perf2$PH, withBNPb$group)
rpart.perfs2<- ROCR::performance(rpart.pred2,"tpr","fpr")
rpart.sens2<- ROCR::performance(rpart.pred2,"sens","spec")
rpart.auc2<- pROC::auc(withBNPb$group, rpart.perf2[,2])
rpart.aucCI2<- pROC::ci.auc(withBNPb$group, rpart.perf2[,2])
par(mfrow=c(1,2))
plot(rpart.perfs2, main = paste("AUC:", round(rpart.auc2,2)))
plot(rpart.sens2)
```

![](ensembleshort_files/figure-html/unnamed-chunk-69-1.png)<!-- -->

```r
rpart.aucCI2
```

```
## 95% CI: 0.6316-0.9934 (DeLong)
```

## LASSO 


```r
LASSO.probs.min2 <- predict(LASSO.min2, newdata=withBNPb[-1], type = "prob")
LASSO.pred.min2<- prediction(LASSO.probs.min2$PH, withBNPb$group)
LASSO.perfs.min2<- ROCR::performance(LASSO.pred.min2,"tpr","fpr")
LASSO.sens.min2<- ROCR::performance(LASSO.pred.min2,"sens","spec")
LASSO.auc.min2<- pROC::auc(withBNPb$group, LASSO.probs.min2[,2])
LASSO.aucCI2<- pROC::ci.auc(withBNPb$group, LASSO.probs.min2[,2])
par(mfrow=c(1,2))
plot(LASSO.perfs.min2, main = paste("AUC:", round(LASSO.auc.min2,2)))
plot(LASSO.sens.min2)
```

![](ensembleshort_files/figure-html/unnamed-chunk-70-1.png)<!-- -->

```r
LASSO.aucCI2
```

```
## 95% CI: 0.8296-1 (DeLong)
```

## XGboost


```r
xgbpreds2<- predict(xgb_tune_final_short2, validation2, type = "prob")
xgb.pred2<- prediction(xgbpreds2$PH, as.factor(labelsB))
xgb.perfs2<- ROCR::performance(xgb.pred2,"tpr","fpr")
xgb.sens2<- ROCR::performance(xgb.pred2,"sens","spec")
xgb.auc2<- pROC::auc(as.factor(labelsB), xgbpreds2[,2])
xgb.aucCI2<- pROC::ci.auc(as.factor(labelsB), xgbpreds2[,2])
par(mfrow=c(1,2))
plot(xgb.perfs2, main = paste("AUC:", round(xgb.auc2,2)))
plot(xgb.sens2)
```

![](ensembleshort_files/figure-html/unnamed-chunk-71-1.png)<!-- -->

## Ensemble


```r
totalprobs2<- data.frame(patientID = rownames(Boruta.perf2), group = tree.Boruta2$group, RandomForest = Boruta.perf2$PH, rpart = rpart.perf2$PH, LASSO = LASSO.probs.min2$PH, XGB = xgbpreds2$PH)
totalprobs2 <- mutate(totalprobs2, Mean = rowMeans(totalprobs2[3:6]))
totalprobs2$EnsemblePreds2 <- ifelse(totalprobs2$Mean > 0.5, "PH", "HV")
Ensemble.pred2<-prediction(totalprobs2$Mean, totalprobs2$group) 
Ensemble.perfs2<- ROCR::performance(Ensemble.pred2, "tpr", "fpr")
kable(totalprobs2)  %>% kable_styling(full_width = TRUE)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> patientID </th>
   <th style="text-align:left;"> group </th>
   <th style="text-align:right;"> RandomForest </th>
   <th style="text-align:right;"> rpart </th>
   <th style="text-align:right;"> LASSO </th>
   <th style="text-align:right;"> XGB </th>
   <th style="text-align:right;"> Mean </th>
   <th style="text-align:left;"> EnsemblePreds2 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 988_2_1.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.741 </td>
   <td style="text-align:right;"> 0.9285714 </td>
   <td style="text-align:right;"> 0.5754949 </td>
   <td style="text-align:right;"> 0.8038322 </td>
   <td style="text-align:right;"> 0.7622246 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 987_2_2.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.856 </td>
   <td style="text-align:right;"> 0.9285714 </td>
   <td style="text-align:right;"> 0.9999971 </td>
   <td style="text-align:right;"> 0.9065028 </td>
   <td style="text-align:right;"> 0.9227678 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 988_2_2.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.532 </td>
   <td style="text-align:right;"> 0.0625000 </td>
   <td style="text-align:right;"> 0.5965555 </td>
   <td style="text-align:right;"> 0.5060492 </td>
   <td style="text-align:right;"> 0.4242762 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 991_1_1.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.908 </td>
   <td style="text-align:right;"> 0.9285714 </td>
   <td style="text-align:right;"> 0.9998723 </td>
   <td style="text-align:right;"> 0.9368197 </td>
   <td style="text-align:right;"> 0.9433159 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 756_1_4.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.913 </td>
   <td style="text-align:right;"> 0.9285714 </td>
   <td style="text-align:right;"> 0.9636261 </td>
   <td style="text-align:right;"> 0.9433010 </td>
   <td style="text-align:right;"> 0.9371246 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 756_1_1.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.892 </td>
   <td style="text-align:right;"> 0.9285714 </td>
   <td style="text-align:right;"> 0.9472897 </td>
   <td style="text-align:right;"> 0.8939501 </td>
   <td style="text-align:right;"> 0.9154528 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 990_1_4.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.950 </td>
   <td style="text-align:right;"> 0.9285714 </td>
   <td style="text-align:right;"> 0.8271773 </td>
   <td style="text-align:right;"> 0.9509395 </td>
   <td style="text-align:right;"> 0.9141721 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 990_1_3.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.744 </td>
   <td style="text-align:right;"> 0.9285714 </td>
   <td style="text-align:right;"> 0.4253599 </td>
   <td style="text-align:right;"> 0.7725604 </td>
   <td style="text-align:right;"> 0.7176229 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 759_2_4.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.941 </td>
   <td style="text-align:right;"> 0.9285714 </td>
   <td style="text-align:right;"> 0.9949804 </td>
   <td style="text-align:right;"> 0.8888467 </td>
   <td style="text-align:right;"> 0.9383496 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 760_1_1.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.917 </td>
   <td style="text-align:right;"> 0.9285714 </td>
   <td style="text-align:right;"> 0.9741941 </td>
   <td style="text-align:right;"> 0.8625593 </td>
   <td style="text-align:right;"> 0.9205812 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 990_1_1.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.848 </td>
   <td style="text-align:right;"> 0.9285714 </td>
   <td style="text-align:right;"> 0.9999998 </td>
   <td style="text-align:right;"> 0.8455595 </td>
   <td style="text-align:right;"> 0.9055327 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 756_2_2.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.746 </td>
   <td style="text-align:right;"> 0.9285714 </td>
   <td style="text-align:right;"> 0.9991012 </td>
   <td style="text-align:right;"> 0.8210467 </td>
   <td style="text-align:right;"> 0.8736798 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 991_2_3.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.741 </td>
   <td style="text-align:right;"> 0.0625000 </td>
   <td style="text-align:right;"> 0.6952485 </td>
   <td style="text-align:right;"> 0.5080644 </td>
   <td style="text-align:right;"> 0.5017032 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 989_2_2.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.787 </td>
   <td style="text-align:right;"> 0.9285714 </td>
   <td style="text-align:right;"> 0.7774324 </td>
   <td style="text-align:right;"> 0.9088891 </td>
   <td style="text-align:right;"> 0.8504732 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 760_2_1.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.694 </td>
   <td style="text-align:right;"> 0.9285714 </td>
   <td style="text-align:right;"> 0.6341162 </td>
   <td style="text-align:right;"> 0.7872966 </td>
   <td style="text-align:right;"> 0.7609961 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 770_1_2.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.899 </td>
   <td style="text-align:right;"> 0.9285714 </td>
   <td style="text-align:right;"> 0.9441111 </td>
   <td style="text-align:right;"> 0.9493203 </td>
   <td style="text-align:right;"> 0.9302507 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 779_1_2.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.767 </td>
   <td style="text-align:right;"> 0.9285714 </td>
   <td style="text-align:right;"> 0.4301157 </td>
   <td style="text-align:right;"> 0.7762326 </td>
   <td style="text-align:right;"> 0.7254799 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 986_1_4.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.529 </td>
   <td style="text-align:right;"> 0.0625000 </td>
   <td style="text-align:right;"> 0.4565787 </td>
   <td style="text-align:right;"> 0.4813725 </td>
   <td style="text-align:right;"> 0.3823628 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 994_1_4.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.246 </td>
   <td style="text-align:right;"> 0.0625000 </td>
   <td style="text-align:right;"> 0.1335094 </td>
   <td style="text-align:right;"> 0.1271603 </td>
   <td style="text-align:right;"> 0.1422924 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 980_1_1.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.256 </td>
   <td style="text-align:right;"> 0.0625000 </td>
   <td style="text-align:right;"> 0.2953754 </td>
   <td style="text-align:right;"> 0.2057881 </td>
   <td style="text-align:right;"> 0.2049159 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 980_1_2.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.142 </td>
   <td style="text-align:right;"> 0.0625000 </td>
   <td style="text-align:right;"> 0.2926489 </td>
   <td style="text-align:right;"> 0.2262291 </td>
   <td style="text-align:right;"> 0.1808445 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 982_1_3.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.313 </td>
   <td style="text-align:right;"> 0.0625000 </td>
   <td style="text-align:right;"> 0.2995181 </td>
   <td style="text-align:right;"> 0.2219070 </td>
   <td style="text-align:right;"> 0.2242313 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 758_1_2.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.219 </td>
   <td style="text-align:right;"> 0.0625000 </td>
   <td style="text-align:right;"> 0.2667769 </td>
   <td style="text-align:right;"> 0.1904820 </td>
   <td style="text-align:right;"> 0.1846897 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 982_1_1.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.371 </td>
   <td style="text-align:right;"> 0.9285714 </td>
   <td style="text-align:right;"> 0.5217201 </td>
   <td style="text-align:right;"> 0.4548452 </td>
   <td style="text-align:right;"> 0.5690342 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
</tbody>
</table>


```r
pROC::auc(totalprobs2$group, totalprobs2$Mean)
```

```
## Area under the curve: 0.9375
```

```r
pROC::ci.auc(totalprobs2$group, totalprobs2$Mean)
```

```
## 95% CI: 0.8432-1 (DeLong)
```

```r
table(totalprobs2$EnsemblePreds2, totalprobs2$group)
```

```
##     
##      HV PH
##   HV  6  1
##   PH  2 15
```

## Figure 3

Dashed line models include NTproBNP


```r
layout(mat = matrix(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,7,7), nrow = 2, byrow = TRUE))
#layout.show(n = 5)
plot(xgb.perfs, col = "dark blue", lty = 1, main = "XGBoost", cex.lab = 1.3, cex.main = 1.5)
plot(xgb.perfs2, col = "dark blue", lty = 2, add = TRUE, cex.lab = 2, cex.main = 1.5)
plot(LASSO.perfs.min, col = "dark blue", lty = 1, main = "LASSO", cex.lab = 1.3, cex.main = 1.5)
plot(LASSO.perfs.min2, col = "dark blue", lty = 2, add = TRUE, cex.lab = 1.3, cex.main = 1.5)
plot(rpart.perfs, col = "dark blue", lty = 1, main = "Rpart", cex.lab = 1.3, cex.main = 1.5)
plot(rpart.perfs2, col = "dark blue", lty = 2, add = TRUE, cex.lab = 1.3, cex.main = 1.5)
plot(Boruta.perfs, col = "dark blue", lty = 1, main = "Random Forest", cex.lab = 1.3, cex.main = 1.5)
plot(Boruta.perfs2, col = "dark blue", lty = 2, add = TRUE, cex.lab = 1.3, cex.main = 1.5)
plot(Ensemble.perfs, col = "dark blue", lty = 1, main = "Ensemble", cex.lab = 1.3, cex.main = 1.5)
plot(Ensemble.perfs2, col = "dark blue", lty = 2, add = TRUE, cex.lab = 1.3, cex.main = 1.5)
plot(perfs.BNP, main = "NTproBNP", cex.lab = 1.3, cex.main = 1.5, col = "dark blue", lty = 2)
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", legend = c("base model (n = 32)", "with NTproBNP (n = 24)"), lty = c(1,2), cex = 1.3, bty = "n") 
```

![](ensembleshort_files/figure-html/unnamed-chunk-74-1.png)<!-- -->


```r
layout(mat = matrix(c(1,1,2,2,3,3,4,4,5,5,6,6), nrow = 2, byrow = TRUE))
#layout.show(n = 5)
plot(xgb.perfs, col = "dark blue", lty = 1, main = "XGBoost", cex.lab = 1.3, cex.main = 1.5, lwd = 3)
plot(LASSO.perfs.min, col = "dark blue", lty = 1, main = "LASSO", cex.lab = 1.3, cex.main = 1.5, lwd = 3)
plot(rpart.perfs, col = "dark blue", lty = 1, main = "Rpart", cex.lab = 1.3, cex.main = 1.5, lwd = 3)
plot(Boruta.perfs, col = "dark blue", lty = 1, main = "Random Forest", cex.lab = 1.3, cex.main = 1.5, lwd = 3)
plot(Ensemble.perfs, col = "dark blue", lty = 1, main = "Ensemble", cex.lab = 1.3, cex.main = 1.5, lwd = 3)
plot(perfs.BNP, main = "NTproBNP", cex.lab = 1.3, cex.main = 1.5, col = "dark blue", lty = 1, lwd = 3)
```

![](ensembleshort_files/figure-html/unnamed-chunk-75-1.png)<!-- -->



# Cross Validations on training set


```r
separatesplittoconfusion<- function(newlist, splits, d) {
  newlist[[d]]<- confusionMatrix(splits[[d]]$pred, splits[[d]]$obs, positive = "PH")
}
confusionaccuracy<- function(x){
  (x$table[1]+x$table[4])/(x$table[1]+x$table[2]+x$table[3]+x$table[4])
}
confusionsensitivity<- function(x){
   x$table[1]/(x$table[1]+x$table[2])
}
confusionspecificity<- function(x){
 (x$table[4]/(x$table[3]+x$table[4]))
}
confusionPPV<- function(x){
 (x$table[1]/(x$table[1]+x$table[3]))
}
confusionNPV<- function(x){
  (x$table[4]/(x$table[2]+x$table[4]))
}
```

## Random Forest


```r
borutaconfusions<- list()
borutaconfusions<- lapply(1:100, separatesplittoconfusion, newlist = borutaconfusions, splits = boruta.split)

cvborutaaccuracy<- unlist(lapply(borutaconfusions, confusionaccuracy))
cvborutasensitivity<- unlist(lapply(borutaconfusions, confusionsensitivity))
cvborutaspecificity<- unlist(lapply(borutaconfusions, confusionspecificity))
cvborutaPPV <- unlist(lapply(borutaconfusions, confusionPPV))
cvborutaNPV<- unlist(lapply(borutaconfusions, confusionNPV))
summary(cvborutaaccuracy)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.4286  0.8125  0.8571  0.8461  0.8750  1.0000
```

```r
summary(cvborutasensitivity)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.0000  0.6667  0.6667  0.7250  1.0000  1.0000
```

```r
summary(cvborutaspecificity)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.5000  0.8000  1.0000  0.9305  1.0000  1.0000
```

```r
summary(cvborutaPPV)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
##  0.3333  0.7500  1.0000  0.9056  1.0000  1.0000       2
```

```r
summary(cvborutaNPV)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.5000  0.8000  0.8000  0.8485  1.0000  1.0000
```

## rpart


```r
rpartconfusions<- list()
rpartconfusions<- lapply(1:100, separatesplittoconfusion, newlist = rpartconfusions, splits = rpart.split)
cvrpartaccuracy<- unlist(lapply(rpartconfusions, confusionaccuracy))
cvrpartsensitivity<- unlist(lapply(rpartconfusions, confusionsensitivity))
cvrpartspecificity<- unlist(lapply(rpartconfusions, confusionspecificity))
cvrpartPPV <- unlist(lapply(rpartconfusions, confusionPPV))
cvrpartNPV<- unlist(lapply(rpartconfusions, confusionNPV))
summary(cvrpartaccuracy)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.2500  0.4762  0.5714  0.5774  0.6667  0.8571
```

```r
summary(cvrpartsensitivity)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.0000  0.3333  0.5000  0.4939  0.6667  1.0000
```

```r
summary(cvrpartspecificity)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.2000  0.4917  0.6667  0.6358  0.8083  1.0000
```

```r
summary(cvrpartPPV)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
##  0.0000  0.4000  0.5000  0.5115  0.6250  1.0000       3
```

```r
summary(cvrpartNPV)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.3333  0.5714  0.6667  0.6474  0.7143  1.0000
```

## LASSO


```r
LASSOconfusions<- list()
LASSOconfusions<- lapply(1:100, separatesplittoconfusion, newlist = LASSOconfusions, splits = LASSO.split.min)
cvLASSOaccuracy<- unlist(lapply(LASSOconfusions, confusionaccuracy))
cvLASSOsensitivity<- unlist(lapply(LASSOconfusions, confusionsensitivity))
cvLASSOspecificity<- unlist(lapply(LASSOconfusions, confusionspecificity))
cvLASSOPPV <- unlist(lapply(LASSOconfusions, confusionPPV))
cvLASSONPV<- unlist(lapply(LASSOconfusions, confusionNPV))
summary(cvLASSOaccuracy)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.2500  0.7024  0.7500  0.7596  0.8571  1.0000
```

```r
summary(cvLASSOsensitivity)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.0000  0.3333  0.6667  0.6433  1.0000  1.0000
```

```r
summary(cvLASSOspecificity)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.4000  0.7500  0.8000  0.8375  1.0000  1.0000
```

```r
summary(cvLASSOPPV)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
##  0.0000  0.6000  0.7500  0.7633  1.0000  1.0000       2
```

```r
summary(cvLASSONPV)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.4000  0.6667  0.8000  0.7951  1.0000  1.0000
```

## XGBoost


```r
XGBconfusions<- list()
XGBconfusions<- lapply(1:100, separatesplittoconfusion, newlist = XGBconfusions, splits = XGB.split)
cvXGBaccuracy<- unlist(lapply(XGBconfusions, confusionaccuracy))
cvXGBsensitivity<- unlist(lapply(XGBconfusions, confusionsensitivity))
cvXGBspecificity<- unlist(lapply(XGBconfusions, confusionspecificity))
cvXGBPPV <- unlist(lapply(XGBconfusions, confusionPPV))
cvXGBNPV<- unlist(lapply(XGBconfusions, confusionNPV))
summary(cvXGBaccuracy)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.4286  0.7143  0.8571  0.8296  1.0000  1.0000
```

```r
summary(cvXGBsensitivity)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.0000  0.6667  0.6667  0.7533  1.0000  1.0000
```

```r
summary(cvXGBspecificity)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.2500  0.7500  1.0000  0.8845  1.0000  1.0000
```

```r
summary(cvXGBPPV)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
##  0.3333  0.6667  1.0000  0.8473  1.0000  1.0000       2
```

```r
summary(cvXGBNPV)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.5000  0.7500  0.8000  0.8571  1.0000  1.0000
```


# ROC tests


```r
roc.out_NTproBNP<- roc(as.factor(withBNPb$group), NTproBNPprob)
roc.out_XGB_NTproBNP<- roc(as.factor(labelsB), xgbpreds2[,2])
roc.out_RF_NTproBNP <- roc(tree.Boruta2$group, Boruta.perf2[,2])
roc.out_rpart_NTproBNP <- roc(withBNPb$group, rpart.perf2[,2])
roc.out_LASSO_NTproBNP <- roc(withBNPb$group, LASSO.probs.min2[,2])

roc.out_rpart <- roc(setB2$group, rpart.perf[,2])
roc.out_LASSO <- roc(val.LASSO$group, LASSO.probs.min[,2])
roc.out_XGB <- roc(as.factor(setBshort$group), xgbpreds[,2])
roc.out_RF <- roc(tree.Boruta$group, Boruta.perf[,2])

roc.test(roc.out_RF, roc.out_RF_NTproBNP, paired = FALSE)
```

```
## 
## 	DeLong's test for two ROC curves
## 
## data:  roc.out_RF and roc.out_RF_NTproBNP
## D = -1.5516, df = 42.469, p-value = 0.1282
## alternative hypothesis: true difference in AUC is not equal to 0
## sample estimates:
## AUC of roc1 AUC of roc2 
##   0.8441558   0.9726562
```

```r
roc.test(roc.out_rpart, roc.out_rpart_NTproBNP, paired = FALSE)
```

```
## 
## 	DeLong's test for two ROC curves
## 
## data:  roc.out_rpart and roc.out_rpart_NTproBNP
## D = -0.19268, df = 51.469, p-value = 0.848
## alternative hypothesis: true difference in AUC is not equal to 0
## sample estimates:
## AUC of roc1 AUC of roc2 
##    0.788961    0.812500
```

```r
roc.test(roc.out_LASSO, roc.out_LASSO_NTproBNP, paired = FALSE)
```

```
## 
## 	DeLong's test for two ROC curves
## 
## data:  roc.out_LASSO and roc.out_LASSO_NTproBNP
## D = -1.4974, df = 55.603, p-value = 0.1399
## alternative hypothesis: true difference in AUC is not equal to 0
## sample estimates:
## AUC of roc1 AUC of roc2 
##   0.7889610   0.9296875
```

```r
roc.test(roc.out_XGB, roc.out_XGB_NTproBNP, paired = FALSE)
```

```
## 
## 	DeLong's test for two ROC curves
## 
## data:  roc.out_XGB and roc.out_XGB_NTproBNP
## D = -1.3993, df = 50.743, p-value = 0.1678
## alternative hypothesis: true difference in AUC is not equal to 0
## sample estimates:
## AUC of roc1 AUC of roc2 
##   0.8214286   0.9531250
```

```r
roc.test(roc.out_NTproBNP, roc.out_RF_NTproBNP, paired = FALSE)
```

```
## Warning in roc.test.roc(roc.out_NTproBNP, roc.out_RF_NTproBNP, paired = FALSE):
## The ROC curves seem to be paired. Consider performing a paired roc.test.
```

```
## 
## 	DeLong's test for two ROC curves
## 
## data:  roc.out_NTproBNP and roc.out_RF_NTproBNP
## D = -0.62677, df = 34.993, p-value = 0.5349
## alternative hypothesis: true difference in AUC is not equal to 0
## sample estimates:
## AUC of roc1 AUC of roc2 
##   0.9375000   0.9726562
```

```r
roc.test(roc.out_NTproBNP, roc.out_rpart_NTproBNP, paired = FALSE)
```

```
## Warning in roc.test.roc(roc.out_NTproBNP, roc.out_rpart_NTproBNP, paired =
## FALSE): The ROC curves seem to be paired. Consider performing a paired roc.test.
```

```
## 
## 	DeLong's test for two ROC curves
## 
## data:  roc.out_NTproBNP and roc.out_rpart_NTproBNP
## D = 1.1932, df = 35.241, p-value = 0.2408
## alternative hypothesis: true difference in AUC is not equal to 0
## sample estimates:
## AUC of roc1 AUC of roc2 
##      0.9375      0.8125
```

```r
roc.test(roc.out_NTproBNP, roc.out_LASSO_NTproBNP, paired = FALSE)
```

```
## Warning in roc.test.roc(roc.out_NTproBNP, roc.out_LASSO_NTproBNP, paired =
## FALSE): The ROC curves seem to be paired. Consider performing a paired roc.test.
```

```
## 
## 	DeLong's test for two ROC curves
## 
## data:  roc.out_NTproBNP and roc.out_LASSO_NTproBNP
## D = 0.10982, df = 45.96, p-value = 0.913
## alternative hypothesis: true difference in AUC is not equal to 0
## sample estimates:
## AUC of roc1 AUC of roc2 
##   0.9375000   0.9296875
```

```r
roc.test(roc.out_NTproBNP, roc.out_XGB_NTproBNP, paired = FALSE)
```

```
## Warning in roc.test.roc(roc.out_NTproBNP, roc.out_XGB_NTproBNP, paired = FALSE):
## The ROC curves seem to be paired. Consider performing a paired roc.test.
```

```
## 
## 	DeLong's test for two ROC curves
## 
## data:  roc.out_NTproBNP and roc.out_XGB_NTproBNP
## D = -0.23747, df = 45.185, p-value = 0.8134
## alternative hypothesis: true difference in AUC is not equal to 0
## sample estimates:
## AUC of roc1 AUC of roc2 
##    0.937500    0.953125
```




# Validation gender split AUC


```r
male<- filter(pheno, Gender == "Male" & AorB == "B") %>% select(filename) %>% as.list(.)
female<- filter(pheno, Gender == "Female" & AorB == "B") %>% select(filename) 

Male<- setB[rownames(setB) %in% male[[1]],]
Female<- setB[rownames(setB) %in% female[[1]],]
```

### Random Forest


```r
RFBoruta.trainM <- predict(fit.Boruta, Male)
BorutaInterimM<- confusionMatrix(as.factor(RFBoruta.trainM), as.factor(Male$group), positive = "PH")
BorutaInterimM
```

```
## Confusion Matrix and Statistics
## 
##           Reference
## Prediction HV PH
##         HV  3  1
##         PH  1  7
##                                           
##                Accuracy : 0.8333          
##                  95% CI : (0.5159, 0.9791)
##     No Information Rate : 0.6667          
##     P-Value [Acc > NIR] : 0.1811          
##                                           
##                   Kappa : 0.625           
##                                           
##  Mcnemar's Test P-Value : 1.0000          
##                                           
##             Sensitivity : 0.8750          
##             Specificity : 0.7500          
##          Pos Pred Value : 0.8750          
##          Neg Pred Value : 0.7500          
##              Prevalence : 0.6667          
##          Detection Rate : 0.5833          
##    Detection Prevalence : 0.6667          
##       Balanced Accuracy : 0.8125          
##                                           
##        'Positive' Class : PH              
## 
```

```r
RF.perfM<- predict(fit.Boruta, Male[-1], type= "prob")
RF.aucM<- pROC::auc(Male$group, RF.perfM[,2])
RF.aucCIM<- pROC::ci.auc(Male$group, RF.perfM[,2])
RF.aucM
```

```
## Area under the curve: 0.9375
```

```r
RF.aucCIM
```

```
## 95% CI: 0.7911-1 (DeLong)
```


```r
RFBoruta.trainF <- predict(fit.Boruta, Female)
BorutaInterimF<- confusionMatrix(as.factor(RFBoruta.trainF), as.factor(Female$group), positive = "PH")
BorutaInterimF
```

```
## Confusion Matrix and Statistics
## 
##           Reference
## Prediction HV PH
##         HV  7  2
##         PH  3 11
##                                          
##                Accuracy : 0.7826         
##                  95% CI : (0.563, 0.9254)
##     No Information Rate : 0.5652         
##     P-Value [Acc > NIR] : 0.02627        
##                                          
##                   Kappa : 0.5525         
##                                          
##  Mcnemar's Test P-Value : 1.00000        
##                                          
##             Sensitivity : 0.8462         
##             Specificity : 0.7000         
##          Pos Pred Value : 0.7857         
##          Neg Pred Value : 0.7778         
##              Prevalence : 0.5652         
##          Detection Rate : 0.4783         
##    Detection Prevalence : 0.6087         
##       Balanced Accuracy : 0.7731         
##                                          
##        'Positive' Class : PH             
## 
```

```r
RF.perfF<- predict(fit.Boruta, Female[-1], type= "prob")
RF.aucF<- pROC::auc(Female$group, RF.perfF[,2])
```

```
## Setting levels: control = HV, case = PH
```

```
## Setting direction: controls < cases
```

```r
RF.aucCIF<- pROC::ci.auc(Female$group, RF.perfF[,2])
```

```
## Setting levels: control = HV, case = PH
## Setting direction: controls < cases
```

```r
RF.aucF
```

```
## Area under the curve: 0.8154
```

```r
RF.aucCIF
```

```
## 95% CI: 0.5994-1 (DeLong)
```

### rpart


```r
predrpartM <- predict(fit.caret.rpart, Male)
confusionMatrix(Male$group, predrpartM)
```

```
## Confusion Matrix and Statistics
## 
##           Reference
## Prediction HV PH
##         HV  2  2
##         PH  1  7
##                                           
##                Accuracy : 0.75            
##                  95% CI : (0.4281, 0.9451)
##     No Information Rate : 0.75            
##     P-Value [Acc > NIR] : 0.6488          
##                                           
##                   Kappa : 0.4             
##                                           
##  Mcnemar's Test P-Value : 1.0000          
##                                           
##             Sensitivity : 0.6667          
##             Specificity : 0.7778          
##          Pos Pred Value : 0.5000          
##          Neg Pred Value : 0.8750          
##              Prevalence : 0.2500          
##          Detection Rate : 0.1667          
##    Detection Prevalence : 0.3333          
##       Balanced Accuracy : 0.7222          
##                                           
##        'Positive' Class : HV              
## 
```

```r
rpart.perfM<- predict(fit.caret.rpart, Male[-1], type= "prob")
rpart.aucM<- pROC::auc(Male$group, rpart.perfM[,2])
rpart.aucCIM<- pROC::ci.auc(Male$group, rpart.perfM[,2])
rpart.aucM
```

```
## Area under the curve: 0.625
```

```r
rpart.aucCIM
```

```
## 95% CI: 0.2479-1 (DeLong)
```

```r
predrpartF <- predict(fit.caret.rpart, Female)
confusionMatrix(Female$group, predrpartF, positive = "PH")
```

```
## Confusion Matrix and Statistics
## 
##           Reference
## Prediction HV PH
##         HV  7  3
##         PH  1 12
##                                           
##                Accuracy : 0.8261          
##                  95% CI : (0.6122, 0.9505)
##     No Information Rate : 0.6522          
##     P-Value [Acc > NIR] : 0.05753         
##                                           
##                   Kappa : 0.6378          
##                                           
##  Mcnemar's Test P-Value : 0.61708         
##                                           
##             Sensitivity : 0.8000          
##             Specificity : 0.8750          
##          Pos Pred Value : 0.9231          
##          Neg Pred Value : 0.7000          
##              Prevalence : 0.6522          
##          Detection Rate : 0.5217          
##    Detection Prevalence : 0.5652          
##       Balanced Accuracy : 0.8375          
##                                           
##        'Positive' Class : PH              
## 
```

```r
rpart.perfF<- predict(fit.caret.rpart, Female[-1], type= "prob")
rpart.aucF<- pROC::auc(Female$group, rpart.perfF[,2])
rpart.aucCIF<- pROC::ci.auc(Female$group, rpart.perfF[,2])
rpart.aucF
```

```
## Area under the curve: 0.85
```

```r
rpart.aucCIF
```

```
## 95% CI: 0.6877-1 (DeLong)
```

### LASSO

```r
LASSO.pred.minM <- predict(LASSO.min, newdata=Male[-1])
LASSOInterimM<- confusionMatrix(LASSO.pred.minM, Male$group, positive = "PH")
LASSOInterimM
```

```
## Confusion Matrix and Statistics
## 
##           Reference
## Prediction HV PH
##         HV  3  3
##         PH  1  5
##                                           
##                Accuracy : 0.6667          
##                  95% CI : (0.3489, 0.9008)
##     No Information Rate : 0.6667          
##     P-Value [Acc > NIR] : 0.6315          
##                                           
##                   Kappa : 0.3333          
##                                           
##  Mcnemar's Test P-Value : 0.6171          
##                                           
##             Sensitivity : 0.6250          
##             Specificity : 0.7500          
##          Pos Pred Value : 0.8333          
##          Neg Pred Value : 0.5000          
##              Prevalence : 0.6667          
##          Detection Rate : 0.4167          
##    Detection Prevalence : 0.5000          
##       Balanced Accuracy : 0.6875          
##                                           
##        'Positive' Class : PH              
## 
```

```r
LASSO.perfM<- predict(LASSO.min, Male[-1], type= "prob")
LASSO.aucM<- pROC::auc(Male$group, LASSO.perfM[,2])
LASSO.aucCIM<- pROC::ci.auc(Male$group, LASSO.perfM[,2])
LASSO.aucM
```

```
## Area under the curve: 0.7188
```

```r
LASSO.aucCIM
```

```
## 95% CI: 0.393-1 (DeLong)
```


```r
LASSO.pred.minF <- predict(LASSO.min, newdata=Female[-1])
LASSOInterimF<- confusionMatrix(LASSO.pred.minF, Female$group, positive = "PH")
LASSOInterimF
```

```
## Confusion Matrix and Statistics
## 
##           Reference
## Prediction HV PH
##         HV  6  1
##         PH  4 12
##                                          
##                Accuracy : 0.7826         
##                  95% CI : (0.563, 0.9254)
##     No Information Rate : 0.5652         
##     P-Value [Acc > NIR] : 0.02627        
##                                          
##                   Kappa : 0.5418         
##                                          
##  Mcnemar's Test P-Value : 0.37109        
##                                          
##             Sensitivity : 0.9231         
##             Specificity : 0.6000         
##          Pos Pred Value : 0.7500         
##          Neg Pred Value : 0.8571         
##              Prevalence : 0.5652         
##          Detection Rate : 0.5217         
##    Detection Prevalence : 0.6957         
##       Balanced Accuracy : 0.7615         
##                                          
##        'Positive' Class : PH             
## 
```

```r
LASSO.perfF<- predict(LASSO.min, Female[-1], type= "prob")
LASSO.aucF<- pROC::auc(Female$group, LASSO.perfF[,2])
LASSO.aucCIF<- pROC::ci.auc(Female$group, LASSO.perfF[,2])
LASSO.aucF
```

```
## Area under the curve: 0.8538
```

```r
LASSO.aucCIF
```

```
## 95% CI: 0.6809-1 (DeLong)
```

### XGBoost


```r
xgbpredictfinalshortM<- predict(xgb_tune_final_short, Male[,c("group",xgb_mirs$Feature)])
XGBInterimM<- confusionMatrix(xgbpredictfinalshortM, as.factor(Male$group), positive = "PH")
XGBInterimM
```

```
## Confusion Matrix and Statistics
## 
##           Reference
## Prediction HV PH
##         HV  4  0
##         PH  0  8
##                                      
##                Accuracy : 1          
##                  95% CI : (0.7354, 1)
##     No Information Rate : 0.6667     
##     P-Value [Acc > NIR] : 0.007707   
##                                      
##                   Kappa : 1          
##                                      
##  Mcnemar's Test P-Value : NA         
##                                      
##             Sensitivity : 1.0000     
##             Specificity : 1.0000     
##          Pos Pred Value : 1.0000     
##          Neg Pred Value : 1.0000     
##              Prevalence : 0.6667     
##          Detection Rate : 0.6667     
##    Detection Prevalence : 0.6667     
##       Balanced Accuracy : 1.0000     
##                                      
##        'Positive' Class : PH         
## 
```

```r
XGB.perfM<- predict(xgb_tune_final_short, Male[,c("group",xgb_mirs$Feature)], type= "prob")
XGB.aucM<- pROC::auc(Male$group, XGB.perfM[,2])
XGB.aucCIM<- pROC::ci.auc(Male$group, XGB.perfM[,2])
```

```
## Warning in ci.auc.roc(roc = roc, ...): ci.auc() of a ROC curve with AUC == 1 is
## always 1-1 and can be misleading.
```

```r
XGB.aucM
```

```
## Area under the curve: 1
```

```r
XGB.aucCIM
```

```
## 95% CI: 1-1 (DeLong)
```


```r
xgbpredictfinalshortF<- predict(xgb_tune_final_short, Female[,c("group",xgb_mirs$Feature)])
XGBInterimF<- confusionMatrix(xgbpredictfinalshortF, as.factor(Female$group), positive = "PH")
XGBInterimF
```

```
## Confusion Matrix and Statistics
## 
##           Reference
## Prediction HV PH
##         HV  6  2
##         PH  4 11
##                                           
##                Accuracy : 0.7391          
##                  95% CI : (0.5159, 0.8977)
##     No Information Rate : 0.5652          
##     P-Value [Acc > NIR] : 0.06809         
##                                           
##                   Kappa : 0.4567          
##                                           
##  Mcnemar's Test P-Value : 0.68309         
##                                           
##             Sensitivity : 0.8462          
##             Specificity : 0.6000          
##          Pos Pred Value : 0.7333          
##          Neg Pred Value : 0.7500          
##              Prevalence : 0.5652          
##          Detection Rate : 0.4783          
##    Detection Prevalence : 0.6522          
##       Balanced Accuracy : 0.7231          
##                                           
##        'Positive' Class : PH              
## 
```

```r
XGB.perfF<- predict(xgb_tune_final_short, Female[,c("group",xgb_mirs$Feature)], type= "prob")
XGB.aucF<- pROC::auc(Female$group, XGB.perfF[,2])
XGB.aucCIF<- pROC::ci.auc(Female$group, XGB.perfF[,2])
XGB.aucF
```

```
## Area under the curve: 0.7615
```

```r
XGB.aucCIF
```

```
## 95% CI: 0.5481-0.975 (DeLong)
```




# Model predictions for all samples, training and validation sets


```r
all.Boruta.perf<- predict(fit.Boruta, select(dataset, -group, -PHstatus, -AB, -filename), type= "prob")
all.rpart.perf<- predict(fit.caret.rpart, select(dataset, -group, -PHstatus, -AB, -filename), type= "prob")
all.LASSO.probs.min <- predict(LASSO.min, newdata=select(dataset, -group, -PHstatus, -AB, -filename), type = "prob")
all.xgbpreds<- predict(xgb_tune_final_short, select(dataset, xgb_mirs$Feature, -group, -PHstatus, -AB, -filename), type = "prob")

allprobs<- data.frame(patientID = dataset$filename, group = dataset$group, RandomForest = all.Boruta.perf$PH, rpart = all.rpart.perf$PH, LASSO = all.LASSO.probs.min$PH, XGB = all.xgbpreds$PH)
allprobs <- mutate(allprobs, Mean = rowMeans(allprobs[3:6]))
allprobs$EnsemblePreds <- ifelse(allprobs$Mean > 0.5, "PH", "HV")
#write.csv(allprobs, "~/Google Drive File Stream/My Drive/miRNA/modelpredictions.csv")
kable(allprobs)  %>% kable_styling(full_width = TRUE)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> patientID </th>
   <th style="text-align:left;"> group </th>
   <th style="text-align:right;"> RandomForest </th>
   <th style="text-align:right;"> rpart </th>
   <th style="text-align:right;"> LASSO </th>
   <th style="text-align:right;"> XGB </th>
   <th style="text-align:right;"> Mean </th>
   <th style="text-align:left;"> EnsemblePreds </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 773_1_3.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.861 </td>
   <td style="text-align:right;"> 0.8666667 </td>
   <td style="text-align:right;"> 0.8610238 </td>
   <td style="text-align:right;"> 0.9266345 </td>
   <td style="text-align:right;"> 0.8788313 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 761_2_4.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.975 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.8864699 </td>
   <td style="text-align:right;"> 0.8823745 </td>
   <td style="text-align:right;"> 0.9167303 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 988_2_1.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.625 </td>
   <td style="text-align:right;"> 0.8666667 </td>
   <td style="text-align:right;"> 0.6270985 </td>
   <td style="text-align:right;"> 0.8060075 </td>
   <td style="text-align:right;"> 0.7311932 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 987_2_1.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.861 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.8465702 </td>
   <td style="text-align:right;"> 0.7811339 </td>
   <td style="text-align:right;"> 0.8529452 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 991_1_4.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.204 </td>
   <td style="text-align:right;"> 0.8666667 </td>
   <td style="text-align:right;"> 0.6789439 </td>
   <td style="text-align:right;"> 0.3211375 </td>
   <td style="text-align:right;"> 0.5176870 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 762_1_1.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.952 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.6710947 </td>
   <td style="text-align:right;"> 0.8816717 </td>
   <td style="text-align:right;"> 0.8569608 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 991_1_2.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.865 </td>
   <td style="text-align:right;"> 0.8666667 </td>
   <td style="text-align:right;"> 0.8497567 </td>
   <td style="text-align:right;"> 0.7917249 </td>
   <td style="text-align:right;"> 0.8432871 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 765_1_1.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.926 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.9355463 </td>
   <td style="text-align:right;"> 0.7913492 </td>
   <td style="text-align:right;"> 0.8939931 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 987_2_2.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.867 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.9681849 </td>
   <td style="text-align:right;"> 0.7909521 </td>
   <td style="text-align:right;"> 0.8873035 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 988_2_2.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.642 </td>
   <td style="text-align:right;"> 0.8666667 </td>
   <td style="text-align:right;"> 0.7828652 </td>
   <td style="text-align:right;"> 0.8272272 </td>
   <td style="text-align:right;"> 0.7796898 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 988_2_3.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.469 </td>
   <td style="text-align:right;"> 0.8666667 </td>
   <td style="text-align:right;"> 0.6646667 </td>
   <td style="text-align:right;"> 0.4601494 </td>
   <td style="text-align:right;"> 0.6151207 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 991_1_1.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.892 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.9760982 </td>
   <td style="text-align:right;"> 0.9143730 </td>
   <td style="text-align:right;"> 0.9263870 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 991_2_4.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.949 </td>
   <td style="text-align:right;"> 0.8666667 </td>
   <td style="text-align:right;"> 0.7041091 </td>
   <td style="text-align:right;"> 0.9094258 </td>
   <td style="text-align:right;"> 0.8573004 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 990_2_4.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.909 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.5705039 </td>
   <td style="text-align:right;"> 0.8755173 </td>
   <td style="text-align:right;"> 0.8195245 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 988_2_4.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.301 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.6179052 </td>
   <td style="text-align:right;"> 0.6803067 </td>
   <td style="text-align:right;"> 0.3998030 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 988_1_1.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.266 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.4856227 </td>
   <td style="text-align:right;"> 0.4597225 </td>
   <td style="text-align:right;"> 0.3141999 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 986_1_1.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.866 </td>
   <td style="text-align:right;"> 0.8666667 </td>
   <td style="text-align:right;"> 0.8318733 </td>
   <td style="text-align:right;"> 0.7475660 </td>
   <td style="text-align:right;"> 0.8280265 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 988_1_2.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.908 </td>
   <td style="text-align:right;"> 0.8666667 </td>
   <td style="text-align:right;"> 0.9557704 </td>
   <td style="text-align:right;"> 0.7759580 </td>
   <td style="text-align:right;"> 0.8765988 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 981_1_3.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.145 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.4265124 </td>
   <td style="text-align:right;"> 0.3136015 </td>
   <td style="text-align:right;"> 0.2212785 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 986_1_3.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.877 </td>
   <td style="text-align:right;"> 0.8666667 </td>
   <td style="text-align:right;"> 0.6068253 </td>
   <td style="text-align:right;"> 0.6847930 </td>
   <td style="text-align:right;"> 0.7588212 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 994_2_4.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.112 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.4433199 </td>
   <td style="text-align:right;"> 0.1719547 </td>
   <td style="text-align:right;"> 0.1931823 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 988_1_3.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.244 </td>
   <td style="text-align:right;"> 0.8666667 </td>
   <td style="text-align:right;"> 0.6279423 </td>
   <td style="text-align:right;"> 0.5507385 </td>
   <td style="text-align:right;"> 0.5723369 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 756_1_4.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.827 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.9729316 </td>
   <td style="text-align:right;"> 0.8753925 </td>
   <td style="text-align:right;"> 0.8996003 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 766_1_3.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.075 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.0754960 </td>
   <td style="text-align:right;"> 0.1578765 </td>
   <td style="text-align:right;"> 0.0884568 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 773_1_4.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.115 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.2784155 </td>
   <td style="text-align:right;"> 0.2450149 </td>
   <td style="text-align:right;"> 0.1709712 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 780_2_1.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.588 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.5069344 </td>
   <td style="text-align:right;"> 0.3914426 </td>
   <td style="text-align:right;"> 0.6023635 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 993_2_4.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.753 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.2654226 </td>
   <td style="text-align:right;"> 0.4854265 </td>
   <td style="text-align:right;"> 0.3873259 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 756_1_1.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.805 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.8110336 </td>
   <td style="text-align:right;"> 0.8562466 </td>
   <td style="text-align:right;"> 0.8488393 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 987_1_2.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.915 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.5369077 </td>
   <td style="text-align:right;"> 0.8458510 </td>
   <td style="text-align:right;"> 0.8052089 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 990_1_4.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.908 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.8860309 </td>
   <td style="text-align:right;"> 0.9222433 </td>
   <td style="text-align:right;"> 0.9098378 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 991_1_3.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.804 </td>
   <td style="text-align:right;"> 0.8666667 </td>
   <td style="text-align:right;"> 0.7985939 </td>
   <td style="text-align:right;"> 0.4786575 </td>
   <td style="text-align:right;"> 0.7369795 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 990_1_3.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.671 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.4652240 </td>
   <td style="text-align:right;"> 0.7577598 </td>
   <td style="text-align:right;"> 0.7042652 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 986_2_1.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.877 </td>
   <td style="text-align:right;"> 0.8666667 </td>
   <td style="text-align:right;"> 0.7292580 </td>
   <td style="text-align:right;"> 0.7152265 </td>
   <td style="text-align:right;"> 0.7970378 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 759_2_4.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.918 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.7153382 </td>
   <td style="text-align:right;"> 0.7805382 </td>
   <td style="text-align:right;"> 0.8342383 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 760_1_1.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.888 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.9594336 </td>
   <td style="text-align:right;"> 0.7051297 </td>
   <td style="text-align:right;"> 0.8689100 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 990_1_1.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.799 </td>
   <td style="text-align:right;"> 0.8666667 </td>
   <td style="text-align:right;"> 0.7554185 </td>
   <td style="text-align:right;"> 0.8007538 </td>
   <td style="text-align:right;"> 0.8054598 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 756_2_2.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.621 </td>
   <td style="text-align:right;"> 0.8666667 </td>
   <td style="text-align:right;"> 0.7759468 </td>
   <td style="text-align:right;"> 0.6305977 </td>
   <td style="text-align:right;"> 0.7235528 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 989_2_4.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.859 </td>
   <td style="text-align:right;"> 0.8000000 </td>
   <td style="text-align:right;"> 0.4648375 </td>
   <td style="text-align:right;"> 0.6327240 </td>
   <td style="text-align:right;"> 0.6891404 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 986_1_2.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.906 </td>
   <td style="text-align:right;"> 0.8666667 </td>
   <td style="text-align:right;"> 0.8104660 </td>
   <td style="text-align:right;"> 0.7922019 </td>
   <td style="text-align:right;"> 0.8438336 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 991_2_3.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.858 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.8211853 </td>
   <td style="text-align:right;"> 0.8018932 </td>
   <td style="text-align:right;"> 0.8510389 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 989_2_2.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.713 </td>
   <td style="text-align:right;"> 0.8666667 </td>
   <td style="text-align:right;"> 0.7855566 </td>
   <td style="text-align:right;"> 0.7922395 </td>
   <td style="text-align:right;"> 0.7893657 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 989_2_1.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.916 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.6768916 </td>
   <td style="text-align:right;"> 0.8776831 </td>
   <td style="text-align:right;"> 0.8484129 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 760_1_3.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.865 </td>
   <td style="text-align:right;"> 0.8000000 </td>
   <td style="text-align:right;"> 0.6817949 </td>
   <td style="text-align:right;"> 0.7203315 </td>
   <td style="text-align:right;"> 0.7667816 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 760_1_4.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.865 </td>
   <td style="text-align:right;"> 0.8000000 </td>
   <td style="text-align:right;"> 0.9294178 </td>
   <td style="text-align:right;"> 0.6711552 </td>
   <td style="text-align:right;"> 0.8163932 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 989_1_3.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.786 </td>
   <td style="text-align:right;"> 0.8666667 </td>
   <td style="text-align:right;"> 0.5553606 </td>
   <td style="text-align:right;"> 0.6654003 </td>
   <td style="text-align:right;"> 0.7183569 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 756_2_3.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.844 </td>
   <td style="text-align:right;"> 0.8666667 </td>
   <td style="text-align:right;"> 0.6277201 </td>
   <td style="text-align:right;"> 0.6890227 </td>
   <td style="text-align:right;"> 0.7568524 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 760_2_1.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.710 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.6660852 </td>
   <td style="text-align:right;"> 0.7215491 </td>
   <td style="text-align:right;"> 0.5357722 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 994_2_2.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.047 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.0970803 </td>
   <td style="text-align:right;"> 0.1024557 </td>
   <td style="text-align:right;"> 0.0729976 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 989_1_2.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.885 </td>
   <td style="text-align:right;"> 0.8666667 </td>
   <td style="text-align:right;"> 0.8363723 </td>
   <td style="text-align:right;"> 0.7992688 </td>
   <td style="text-align:right;"> 0.8468269 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 770_1_2.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.877 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.8527485 </td>
   <td style="text-align:right;"> 0.9214918 </td>
   <td style="text-align:right;"> 0.8935793 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 770_1_3.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.981 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.9623633 </td>
   <td style="text-align:right;"> 0.9516933 </td>
   <td style="text-align:right;"> 0.9545334 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 761_1_1.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.964 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.9188469 </td>
   <td style="text-align:right;"> 0.8966982 </td>
   <td style="text-align:right;"> 0.9256555 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 779_1_2.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.713 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.5425231 </td>
   <td style="text-align:right;"> 0.7151947 </td>
   <td style="text-align:right;"> 0.7234487 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 779_1_4.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.631 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.2968497 </td>
   <td style="text-align:right;"> 0.5210212 </td>
   <td style="text-align:right;"> 0.5929869 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 762_2_4.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.978 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.9306949 </td>
   <td style="text-align:right;"> 0.9167427 </td>
   <td style="text-align:right;"> 0.9371286 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 761_1_2.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.980 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.8585149 </td>
   <td style="text-align:right;"> 0.9279700 </td>
   <td style="text-align:right;"> 0.9223905 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 761_1_3.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.950 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.8871988 </td>
   <td style="text-align:right;"> 0.9021694 </td>
   <td style="text-align:right;"> 0.9156113 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 761_1_4.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.987 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.9567790 </td>
   <td style="text-align:right;"> 0.9192775 </td>
   <td style="text-align:right;"> 0.9465334 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 761_2_1.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.976 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.7876303 </td>
   <td style="text-align:right;"> 0.8600428 </td>
   <td style="text-align:right;"> 0.8866875 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 988_1_4.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.896 </td>
   <td style="text-align:right;"> 0.8666667 </td>
   <td style="text-align:right;"> 0.7124172 </td>
   <td style="text-align:right;"> 0.7363076 </td>
   <td style="text-align:right;"> 0.8028479 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 761_2_3.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.922 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.7181012 </td>
   <td style="text-align:right;"> 0.8198347 </td>
   <td style="text-align:right;"> 0.8457532 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 986_1_4.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.540 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.3144224 </td>
   <td style="text-align:right;"> 0.5690590 </td>
   <td style="text-align:right;"> 0.3672340 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 757_2_3.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.063 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.1894561 </td>
   <td style="text-align:right;"> 0.1470265 </td>
   <td style="text-align:right;"> 0.1112343 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 994_1_3.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.240 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.0996630 </td>
   <td style="text-align:right;"> 0.2303311 </td>
   <td style="text-align:right;"> 0.1538622 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 994_1_4.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.251 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.1631164 </td>
   <td style="text-align:right;"> 0.2609853 </td>
   <td style="text-align:right;"> 0.1801391 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 980_1_1.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.327 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.2625634 </td>
   <td style="text-align:right;"> 0.2860165 </td>
   <td style="text-align:right;"> 0.2302586 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 757_2_4.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.147 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.2152335 </td>
   <td style="text-align:right;"> 0.2754282 </td>
   <td style="text-align:right;"> 0.1707790 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 980_1_2.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.251 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.3797412 </td>
   <td style="text-align:right;"> 0.4075463 </td>
   <td style="text-align:right;"> 0.2709355 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 980_1_3.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.119 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.3119756 </td>
   <td style="text-align:right;"> 0.1546872 </td>
   <td style="text-align:right;"> 0.1577793 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 982_2_3.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.147 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.6152616 </td>
   <td style="text-align:right;"> 0.2906246 </td>
   <td style="text-align:right;"> 0.2745852 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 982_2_4.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.085 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.4999167 </td>
   <td style="text-align:right;"> 0.3045401 </td>
   <td style="text-align:right;"> 0.2337278 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 994_1_1.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.177 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.5773281 </td>
   <td style="text-align:right;"> 0.3858711 </td>
   <td style="text-align:right;"> 0.2964134 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 994_1_2.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.234 </td>
   <td style="text-align:right;"> 0.8000000 </td>
   <td style="text-align:right;"> 0.2296237 </td>
   <td style="text-align:right;"> 0.3995842 </td>
   <td style="text-align:right;"> 0.4158020 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 982_1_3.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.489 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.4489525 </td>
   <td style="text-align:right;"> 0.4375283 </td>
   <td style="text-align:right;"> 0.3552338 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 982_1_4.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.121 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.4093854 </td>
   <td style="text-align:right;"> 0.2147404 </td>
   <td style="text-align:right;"> 0.1976451 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 757_2_2.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.073 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.1239309 </td>
   <td style="text-align:right;"> 0.1248079 </td>
   <td style="text-align:right;"> 0.0917983 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 758_1_1.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.124 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.3395365 </td>
   <td style="text-align:right;"> 0.2105049 </td>
   <td style="text-align:right;"> 0.1798740 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 982_2_1.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.022 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.1881976 </td>
   <td style="text-align:right;"> 0.1196183 </td>
   <td style="text-align:right;"> 0.0938176 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 758_1_2.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.298 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.3030106 </td>
   <td style="text-align:right;"> 0.2726333 </td>
   <td style="text-align:right;"> 0.2297746 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 758_1_3.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.240 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.2883087 </td>
   <td style="text-align:right;"> 0.4288145 </td>
   <td style="text-align:right;"> 0.2506444 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 982_1_1.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.420 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.5840516 </td>
   <td style="text-align:right;"> 0.4051780 </td>
   <td style="text-align:right;"> 0.3636710 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 982_1_2.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.070 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.0989643 </td>
   <td style="text-align:right;"> 0.1395528 </td>
   <td style="text-align:right;"> 0.0884929 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 982_2_2.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.012 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.0852925 </td>
   <td style="text-align:right;"> 0.1085789 </td>
   <td style="text-align:right;"> 0.0628315 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 994_2_3.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.112 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.2592133 </td>
   <td style="text-align:right;"> 0.2491712 </td>
   <td style="text-align:right;"> 0.1664598 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 756_1_3.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.988 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.9497630 </td>
   <td style="text-align:right;"> 0.9242942 </td>
   <td style="text-align:right;"> 0.9462835 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 989_1_1.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.879 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.7748558 </td>
   <td style="text-align:right;"> 0.8681318 </td>
   <td style="text-align:right;"> 0.8612661 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 994_2_1.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.850 </td>
   <td style="text-align:right;"> 0.8000000 </td>
   <td style="text-align:right;"> 0.9166502 </td>
   <td style="text-align:right;"> 0.7574031 </td>
   <td style="text-align:right;"> 0.8310133 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 769_1_2.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.032 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.0411933 </td>
   <td style="text-align:right;"> 0.1430631 </td>
   <td style="text-align:right;"> 0.0654277 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 766_2_2.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.204 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.4386853 </td>
   <td style="text-align:right;"> 0.4356677 </td>
   <td style="text-align:right;"> 0.2695882 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 780_2_2.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.151 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.3142257 </td>
   <td style="text-align:right;"> 0.2625172 </td>
   <td style="text-align:right;"> 0.4127049 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 780_1_2.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.079 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.2048144 </td>
   <td style="text-align:right;"> 0.1127005 </td>
   <td style="text-align:right;"> 0.1104923 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 780_1_1.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.396 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.4244571 </td>
   <td style="text-align:right;"> 0.4439179 </td>
   <td style="text-align:right;"> 0.5468630 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 768_1_1.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.108 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.1547524 </td>
   <td style="text-align:right;"> 0.2058445 </td>
   <td style="text-align:right;"> 0.1285129 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 773_2_1.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.942 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.9280066 </td>
   <td style="text-align:right;"> 0.8514320 </td>
   <td style="text-align:right;"> 0.9111289 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 780_2_4.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.302 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.6583041 </td>
   <td style="text-align:right;"> 0.6602639 </td>
   <td style="text-align:right;"> 0.6359112 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 767_2_3.txt </td>
   <td style="text-align:left;"> HV </td>
   <td style="text-align:right;"> 0.492 </td>
   <td style="text-align:right;"> 0.0454545 </td>
   <td style="text-align:right;"> 0.2876439 </td>
   <td style="text-align:right;"> 0.5135296 </td>
   <td style="text-align:right;"> 0.3346570 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 770_2_2.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.893 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.9282741 </td>
   <td style="text-align:right;"> 0.8327190 </td>
   <td style="text-align:right;"> 0.8942675 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 766_1_4.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.858 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.4341256 </td>
   <td style="text-align:right;"> 0.6155569 </td>
   <td style="text-align:right;"> 0.7076899 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 774_1_1.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.982 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.9660126 </td>
   <td style="text-align:right;"> 0.9576942 </td>
   <td style="text-align:right;"> 0.9571959 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 771_1_1.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.975 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.8549635 </td>
   <td style="text-align:right;"> 0.9221686 </td>
   <td style="text-align:right;"> 0.9188022 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 768_1_3.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.875 </td>
   <td style="text-align:right;"> 0.8666667 </td>
   <td style="text-align:right;"> 0.8144125 </td>
   <td style="text-align:right;"> 0.8513563 </td>
   <td style="text-align:right;"> 0.8518589 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 766_2_3.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.958 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.7984063 </td>
   <td style="text-align:right;"> 0.8189564 </td>
   <td style="text-align:right;"> 0.8746099 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 779_2_4.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.953 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.9452917 </td>
   <td style="text-align:right;"> 0.9173815 </td>
   <td style="text-align:right;"> 0.9346875 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 766_1_2.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.861 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.7236454 </td>
   <td style="text-align:right;"> 0.6504859 </td>
   <td style="text-align:right;"> 0.7895521 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 767_2_4.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.255 </td>
   <td style="text-align:right;"> 0.8000000 </td>
   <td style="text-align:right;"> 0.2361473 </td>
   <td style="text-align:right;"> 0.2724249 </td>
   <td style="text-align:right;"> 0.3908930 </td>
   <td style="text-align:left;"> HV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 774_1_3.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.401 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.3690218 </td>
   <td style="text-align:right;"> 0.5977038 </td>
   <td style="text-align:right;"> 0.5727006 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 780_2_3.txt </td>
   <td style="text-align:left;"> PH </td>
   <td style="text-align:right;"> 0.853 </td>
   <td style="text-align:right;"> 0.9230769 </td>
   <td style="text-align:right;"> 0.7862459 </td>
   <td style="text-align:right;"> 0.6164631 </td>
   <td style="text-align:right;"> 0.7946965 </td>
   <td style="text-align:left;"> PH </td>
  </tr>
</tbody>
</table>






