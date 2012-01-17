## Using the sleep dataset

fn = system.file('gender.csv', package="bootES")
myData = read.csv(fn, strip.white=TRUE, header=TRUE)
myData.boot <- bootES(myData, data.col="Meas3", grp.col="Condition",
                      effect.type="r",
                      contrasts = c(A = -40, B = -10, C = 50))

## 1.1.r One Group, One Measure on One Dependent Variable
duration.boot <- bootES(faithful, data.col='eruptions')

## 1.2.r One Group, Two Measures on Two Different Dependent Variables
slope.boot <- bootES(faithful, effect.type="slope")
cor.boot <- bootES(faithful, effect.type="r")

## 1.2.d One Group, Two Repeated Measures [counterbalanced ordering]
data(sleep)

## 1.3+ One Group, Multiple Repeated Measures

## 2. Two Groups

## 2.1 Two Groups, One Measure

## 2.2 Two Groups, Two Repeated Measures

## 2.3+ Two Groups, Multiple Repeated Measures

## 3.1 Three or More Groups, One Measure

## 3.2 Three or More Groups, Two Repeated Measures

## 3.3+ Three or More Groups, Multiple Repeated Measures




