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
