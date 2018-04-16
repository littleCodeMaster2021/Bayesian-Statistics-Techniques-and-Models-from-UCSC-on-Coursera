
rm(list=ls())


#With: Evaluate an R expression in an environment constructed from data, possibly modifying the original data.
#with(data, expr, …)
with(warpbreaks, table(wool, tension))
#glm is used to fit generalized linear models, specified by giving a symbolic description of the 
#linear model and a description of the error distribution family and link function to be used in the model.
with(data.frame(u = c(5,10,15,20,30,40,60,80,100),
                lot1 = c(118,58,42,35,27,25,21,19,18),
                lot2 = c(69,35,26,21,18,16,13,12,12)),
     list(summary(glm(lot1 ~ log(u), family = Gamma)),
          summary(glm(lot2 ~ log(u), family = Gamma))))

#A data frame is split by row index's factors into data frames and function FUN is applied to each subset in turn.
# by(data, INDICES, FUN, ..., simplify = TRUE)
# INDICES warpbreaks[,"tension"]has three levels: "L" "M" "H"
by(warpbreaks[, 1:2], warpbreaks[,"tension"], summary)

by(warpbreaks, warpbreaks[,"tension"],
function(x) lm(breaks ~ wool, data = x))

###############################apply(X, Margin, FUN, ...)#####Margin: 1: row, 2: column#################
x<-matrix(1:12,ncol=3)
apply(x,1,sum) # sum of rows
x <- cbind(x1 = 3, x2 = c(4:1, 2:5))
myFUN<- function(x, c1, c2) {
     c(sum(x[c1],1), mean(x[c2])) 
}

# 让数据框的x1列加1，并计算出x1,x2列的均值。
# Method 1
apply(x,1,myFUN,c1='x1',c2=c('x1','x2'))
# Method 2
data.frame(x1=x[,1]+1,x2=rowMeans(x))
# Method 3
df<-data.frame()
  for(i in 1:nrow(x)){
     row<-x[i,]
    df<-rbind(df,rbind(c(sum(row[1],1), mean(row))))}

# 清空环境变量
 rm(list=ls())

# 封装fun1
 fun1<-function(x){
     myFUN<- function(x, c1, c2) {
         c(sum(x[c1],1), mean(x[c2])) 
       }
     apply(x,1,myFUN,c1='x1',c2=c('x1','x2'))
   }

# 封装fun2
 fun2<-function(x){
     df<-data.frame()
     for(i in 1:nrow(x)){
         row<-x[i,]
         df<-rbind(df,rbind(c(sum(row[1],1), mean(row))))
       }
   }

# 封装fun3
 fun3<-function(x){
     data.frame(x1=x[,1]+1,x2=rowMeans(x))
   }

# 生成数据集
 x <- cbind(x1=3, x2 = c(400:1, 2:500))

# 分别统计3种方法的CPU耗时。
 system.time(fun1(x))

 system.time(fun2(x))


 system.time(fun3(x))
