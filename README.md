# Compare-the-GDP-rate-of-Aus-India-and-Usa-
GDP per capita is gross domestic product divided by midyear population. GDP is the sum of
gross value added by all resident producers in the economy plus any product taxes and minus any
subsidies not included in the value of the products. It is calculated without making deductions for
depreciation of fabricated assets or for depletion and degradation of natural resources. Gross
domestic product (GDP) is the ordinary measurement of the value added created through the
production of goods &amp; services in a country during a certain time period. It also measures the
income earned from that production (manufacturing, development), or the total amount spent on
final goods and services. The proposed dataset for the visualization by R and Tableau by the
metadata of Australia, USA and India. We have concluded that log of GDP&#39;s of India, Usa and
Australia are integrated of order one that make linear structure for forecasting. We then
performed Engle-Granger two step test for integration. After fitting the integrating regression
models, the residuals were tested for stationary by means of ADF test. The results imply that
there is an equilibrating mechanism that keeps the countries&#39; GDP from drifting too far apart
from each other.
Data Source World Bank
Last Updated Date 4/24/2019
#### Scope of the project: By using this project we have compare the GDP rate of Aus, India andUsa as well as we forecast the GDP rate of Aus, India and Usa.
By Visualization we have
compare the trend and found the in which year GDP rate is lowest and Hightest. The project has
been integrate the uncomplicated, intermediary, questions that are solved and figerout the trend of
GDP in USA, India and Australia and the forecast GDP of USA, India and Australia?

Data Preparation
The data is sourced by the World Bank database, which in excel format. I have changed excel to
csv. It has 3 variables and 57 rows, that contain GDP data yearly wise 1960 to 2016.In this data
set has no missing values. So now I have exported the data set to R tool that help to analysis and
visualization the data. A metadata is a dataset that contain all information about other data, thus
this dataset has taken by Worldbank. In this data set contains attributes associated with
Usa,India, Australia and all other countries.The visualization by R project , I have use only
Usa,India, Australia GDP to compare the Gdp and Forcast the GDP . that data include year 1960
to 2017, GDP per capita growth.
The dataset will contain only 3 variables, that all are quantitative variables and The dataset
contains information from 1960 to 2017. In data set no missing values. The data Australia ,India,
Usa was sourced from the Data World Bank website.
I have investigated the dataset of GDP by &quot; WorldBank&quot; . I have focused on the GDP per
capita and tryed to made a model using the data from the dataset. I l also compare the the total

GDPs of Aus,Ind and Usa. I have use of EDA and visualization, build model, forecasting,used
Arima model.

#### Visualization:
The all figures shows the relation between GDP growth with time and year ( y axis ).By R code
we can Assess the Relationship between the Economic Growths of Aus, Usa, and India. The
statistics indicate disparity in GDP, with Usa having the highest GDP. The standard deviation of
Usa GPD is the highest, indicating substantial fluctuation over the time period under study. First
step of our analysis consisted of determining the order of integration of the three variables by
means of the ADF test. The test was first performed with only a drift term and then with constant
and trend for the series in levels and first difference. We have concluded that log of GDP&#39;s of
India, Usa and Australiya are integrated of order one that make linear structure for forcasting.
We then performed Engle-Granger two step test for cointegration. After fitting the cointegrating
regression models, the residuals were tested for stationarity by means of ADF test. The results
imply that there is an equilibrating mechanism that keeps the countries&#39; GDP from drifting too
far apart from each other.



#title: "To Assess the Relationship between the Economic Growths of Aus,India and USA"




##Executive Summary 

#The statistics indicate disparity in GDP, with USA having the highest GDP. The standard deviation of USA GPD is the highest, indicating substantial fluctuation over the time period under study.First step of our analysis consisted of determining the order of integration of the three variables by means of the ADF test. The test was first performed with only a drift term and then with constant and trend for the series in levels and first difference. *We concluded that log of GDP's of India, USA and Aus are integrated of order one.* 


#We then performed Engle-Granger two step test for cointegration. After fitting the cointegrating regression models, the residuals were tested for stationarity by means of ADF test. The results imply that there is an equilibrating mechanism that keeps the countries' GDP from drifting too far apart from each other. 


#When two variables are cointegrated, there should be granger causality in at least one direction.The results of Toda-Yamamoto Granger causality test show strong evidence of causality running from India to Aus. But Aus's GDP does not granger-cause India's GDP. We also found that *USA and India together granger-cause Aus's GDP. One possible reason for this could be the enormous volumes of exports from USA and India to Aus. Exports from Aus to these countries are very small as compared to the imports from these countries.* Based on the results of the causality analysis, we can conclude that India and USA affect Aus's economic growth more than Aus's economic growth affect those countries'.


#In the end, an VECM (Vector Error Correction Model) was fitted. The fit looks good with a MAPE of 0.64 except that the model is slightly over-predicting. 

#Import requuired libraries

library(tseries)
library(vars)
library(urca)
library(forecast)
library(zoo)
library(ggplot2)



gdp<- read.csv("mainDataOFGDP.csv", header = T)
dim(gdp)
gdp<- ts(gdp, start = 1960, frequency = 1)
class(gdp)
head(gdp)

####Summary Statistics (in millions of dollars)
summary(gdp)

####Plots



library(ggplot2); library(ggfortify)
autoplot(gdp)
autoplot(gdp, facets=F)


#We can see a clear increasing trend in GDP of all countries and hence the series cannot be stationary. We will difference the series apply ADF test to check if the series is stationary or not.

#We will try transforming by taking natural logarithm.

#India vs. Aus
plot(gdp[, 2],gdp[,1],pch=20 )     #this is linear
plot(log(gdp[, 2]),log(gdp[,1]),pch=20 )

#India vs. USA
plot(gdp[, 3],gdp[,1],pch=20 )

#various transformations
plot((log(gdp[, 3])),(log(gdp[,1])),pch=20 )
plot((gdp[, 3])^(1/2),(gdp[,1]),pch=20 )   #this is  linear


###

###After making required transformations, scatter plots between the dependent and independent variable seems linear. Scatter plot of India VS. Aus was linear therefore, no transformation was required. Scatter plot of India's and USA's GDP had a slight curve. So, various transformations like log, square root, cube root were tried. Transformation in which USA's GDP was transformed using square root was found to be linear. 


##**2nd attempt** -- Taking log transformations of all the variables as we were getting series which were integrated of different orders. After log transformations all are I(1).

###Now, we will check the stationarity of India's GDP, Aus's GDP and square root of USA's GDP. 

####Checking stationarity

#S_USA<- sqrt(gdp[ ,3])
#S_USA<- ts(S_USA, start = 1960)

autoplot(log(gdp[,3]), ts.colour = "blue",color="red", main= "Time Series Plot of log of USA's GDP",lty=2)
autoplot(log(gdp[ ,c(1,2,3)]), facets = F)
#code for making b/w plot
autoplot(log(gdp[ ,c(1,2,3)]), facets = F, cex = 2)+xlab("Year")+ylab("Log of GDP")+theme(legend.text = element_text(size=10, face="bold"), axis.text.x = element_text(face = "bold", size = 10), axis.text.y = element_text(face = "bold", size = 10),axis.title.y = element_text(face = "bold", size = 11), axis.title.x= element_text(face = "bold", size = 12))+aes(linetype = series)+ scale_linetype_manual(values = c(1,2,3))+scale_color_manual(values = c("black", "black","black"))+scale_size_manual(values = c(10,10,10))+geom_line(size=0.72)


###Log Transform
lgdp<- log(gdp)
summary(ur.df((lgdp[,1]), type="drift", selectlags = "AIC"))
summary(ur.df(lgdp[,1], type="trend", selectlags = "AIC"))
(ur.df(diff(lgdp[,1]), type = "trend", selectlags = "AIC"))
(ur.df(diff(lgdp[,1]), type="drift", selectlags = "AIC"))

(ur.df((lgdp[,2]), type="drift", selectlags = "AIC"))
(ur.df((lgdp[,2]), type="trend", selectlags = "AIC"))
(ur.df(diff(lgdp[,2]), type = "drift", selectlags = "AIC"))
(ur.df(diff(lgdp[,2]), type = "trend", selectlags = "AIC"))

(ur.df((lgdp[,3]), type="drift", selectlags = "AIC"))
(ur.df((lgdp[,3]), type="trend", selectlags = "AIC"))
(ur.df(diff(lgdp[,3]), type = "drift", selectlags = "AIC"))
(ur.df(diff(lgdp[,3]), type = "trend", selectlags = "AIC"))




####How to interpret the results from the 
#Clearly, these series are non-stationary. 



summary(ur.df(diff(tgdp[,1], differences = 1), type="none"))
summary(ur.df(diff(tgdp[,2], type = "none")))
summary(ur.df(diff(tgdp[,3], differences = 1), type="none"))

###The p-value is less than 0.05. We can reject the null hypothesis that the series is non-stationary or that this series contains a unit root. First differences of log og GDP of India, USA and Aus are stationary, hence they are Integrated of order 1. 

####Plots of stationary series after differencing.


par(mfrow=c(3,1), oma=c(0,5,0,5))
plot(diff(log(gdp[,1]), differences = 1), main="Time series plot of Dlog(GDP_India)", type="l")
plot(diff(log(gdp[,2]), differences = 1), main="Time series plot of Dlog(GDP_Aus)", type="l")
plot(diff(log(gdp[,3]), differences = 2), main="Time series plot of Dlog(GDP_USA)", type="l")



####Testing Cointegration

##We will use the Engle Granger 2 step Methodology to test for cointegration as described in *Analysis of Integrated and Co-Integrated Time Series with R*. The long-run relationship for each of the series are simply estimated by Ordinary Least Squares (OLS). The residuals from these three long-run relationships are stored as objects. An Augmented Dickey- Fuller (ADF) Test is applied to the residuals of each equation for testing whether the variables are cointegrated or not. We will use the critical values found in MacKinnon (1991) or Engle and Yoo (1987).


lgdp<- log(gdp)
BoxCox.lambda(gdp)
tgdp<- BoxCox(gdp, lambda= 0.04184914 )


ilm<- summary(lm(India~ Aus + USA, data= lgdp))
clm<- summary(lm(USA~ Aus + India, data= lgdp))
plm<- summary(lm(Aus~ India + USA, data= lgdp))

error.i<- ts(resid(ilm), start = 1960, frequency = 1)
error.c<- ts(resid(clm), start = 1960, end= 2016, frequency = 1)
error.p<- ts(resid(plm), start = 1960, end= 2016, frequency = 1)

{plot(resid(ilm), type="l") 
abline(h=0)}
{plot(resid(clm), type="l")
abline(h=0)}  
{plot(resid(plm), type="l")
abline(h=0)}

summary(ur.df(error.i, selectlags = "AIC", type = "drift"))
summary(ur.df(error.p, selectlags = "AIC", type = "drift"))
summary(ur.df(error.c, selectlags = "AIC", type = "drift"))

punitroot(c(-3.4561, -3.5721,-2.0123), N=57, trend = "nc", statistic = "t")
#0.0008411547 0.0005822858 0.0431940180
punitroot(c(-3.4234,-3.5411, -1.9952), N=57, trend = "c", statistic = "t")


###Critical values from MacKinnon 1991


#beta+ beta1/T + beta2/T^2

#no constant
-1.9393- (0.398/56) - (10.04/56^2)  #-1.949609
#no trend
-3.7429  -8.352/56-13.41/56^2  #-3.896319
#with trend
-4.1193 - (12.024/56)-(13.13/56)   #-4.568479

#beta+ beta1/T + beta2/T^2 + beta3/T^3
#constant term included
-3.74066-(8.5631/56)-(10.852/56^2)+(27.982/56^3) #-3.896874
#trend term included
-4.11890-(11.8922/56)-(19.031/56^2)+(77.332/56^3) #-4.336889


test.i@teststat #-3.456054
test.p@teststat
test.c@teststat

##Test statistic is greater than the critical value of -1.949609, calculated from *Critical Values for Cointegrating Tests* by James G MacKinnon (1991). Hence, we reject the null hypothesis of a unit root and conclude that the residuals are stationary. 

####Model Building
##we will firstly split the data into training and testing dataset.Methodology employed is as described in *Forecasting: Principles and Practice by Rob J Hyndman & George*

train<- window(lgdp, end = 2013)
test<- window(lgdp, start = 2014)

fit<- auto.arima(train[,1], xreg = train[,c(2,3)], approximation = F, stepwise = F)
Box.test(residuals(fit), type="Ljung" )
ggAcf(residuals(fit)); tsdisplay(arima.errors(fit), main="Arima Errors")
tsdisplay(residuals(fit))
fcast<- forecast(fit, xreg= test[,c(2,3)], h=3)
fcast; 
autoplot(fcast, col="black", lty = 3)+autolayer(test[,1], series = "test data", type="p")
autoplot(fcast)


##ljung box test shows that the residuals are uncorrelated

####Comaparing actual to fitted values


{plot(train[,1], ylab = expression('Log of GDP'[India]), xlab="Year")
lines(smooth.spline(fitted(fit)), col="red")
#adding a legend
legend(1960,28, legend = c("Actual value of LGDP", "Fitted Value of LGDP"), col=c("black", "red"), lty=1, cex= 0.8, text.font = 4, bg="grey" )}

#comparing forecasts
{plot(lgdp[,1],  lty= 1, type="l", lwd=1, ylim=c(24,29), xlab="Year", ylab=expression('Log of GDP'[India]))
lines(2014:2016, fcast$mean, pch= 20, col="blue", type="b", lty=3)
lines(2014:2016, fcast$lower[,1], type = "l", lty= 3)
lines(2014:2016, fcast$lower[,2], type = "l", lty= 5)
lines(2014:2016, fcast$upper[,1], type = "l", lty= 3)
lines(2014:2016, fcast$upper[,2], type = "l", lty= 5)
legend("topleft", legend = c("Actual","point forecasts","80% conf. int.","90% conf. int."), lty=c(1,3,3,5), cex= 0.8, text.font = 4, pch= c(NA, 20,NA, NA ), bg="grey", col=c("black", "blue","black","black"))}
legend(1960,28, legend = "point forecasts", pch=20, text.font = 3, col="blue")


####Accuracy and Model Diagnostics

accuracy(fcast, lgdp[,1])
accuracy(fcast, lgdp[,1])["Test Set","MAPE"]
checkresiduals(residuals(fit))


####Error Correction Model
#The necessary first differences of the series and its lagged values are
#created, as well as the series for the error term lagged by one period.


train1<- ts(embed(diff(train), dimension = 2), start = 1962)
colnames(train1)<- c("DInd","DPak","DChn","D1Ind","D1Pak","D1Chn") 
error_ecm<- window(lag(error.i, k=-1), start = 1962)
  
ecm <??? lm( DInd ~ error_ecm[-53] + D1Ind + D1Pak + D1Chn , data = train1 )


####Modeling GDP with deterministic and stochastic trends.

#Aus
train_p<- window(gdp[,2], end=2014)
test_p<- window(gdp[,2], start= 2015)
#deterministic trend
trend<- seq_along(train_p)
model1<- auto.arima(train_p, d=0, xreg=1960:2014)
summary(model1)

#stochastic trend
model2<- auto.arima(train_p, d=1)
summary(model2)
autoplot(train_p)+autolayer(fitted(model1))
autoplot(train_c)+autolayer(fitted(model3))

train_c<- window(gdp[,3], end=2014)
test_c<- window(gdp[,3], start= 2015)
model3<- auto.arima(train_c, d=1)
summary(model3)

model4<- auto.arima(train_p, stepwise = F, approximation = F)
autoplot(train_p)+autolayer(fitted(model4))
f_p<- forecast(model4, h=3)
accuracy(f_p, test_p)
f_p; test_p

model5<- auto.arima(train_c, stepwise = F, approximation = F)
autoplot(train_c)+autolayer(fitted(model5))
f_c<- forecast(model5, h=3)
accuracy(f_c, test_c)
f_c; test_c

f_c$mean[3]; f_p$mean[3]
xreg1<- data.frame(Aus=log(f_p$mean[3]), USA=log(f_c$mean[3]))

forecast(fit, xreg = xreg1  );test

tail(gdp[,1]); exp(28.47464)



####Causality tests
#Null hypothesis would be that Dln(Chn) and Dln(Pak) does not granger cause Dln(Ind)


#library(vars)
#series must be stationary, make stationary by taking first difference
DlnChn<- diff(lgdp[,3])
DlnPak<- diff(lgdp[,2])
DlnInd<- diff(lgdp[,1])

DlnGDP<- cbind(DlnInd, DlnPak, DlnChn)
D2<- cbind(DlnInd, DlnPak)
D3<- cbind(DlnInd, DlnChn)

View(DlnGDP)

#Granger causality test is very sensitive to lag length therefore we choose
#optimal lag length by AIC 
VARselect(DlnGDP, lag.max= 5, type="const")
vargdp<- VAR(DlnGDP, type="const", p=1)
vard2<- VAR(D2, type="const", ic="AIC")
vard3<- VAR(D3, type="const", ic="AIC")
vard2$p
vard3

#causality test
#granger and instantaneous
#DlnChn and DlnPak does not cause DlnInd
causality(vargdp, cause = c("DlnInd"))
causality(vard2, cause="DlnInd")
causality(vargdp, cause = "DlnPak")

