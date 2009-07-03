#A "MMIXclass" class has been created, in order to print, summary, and plot the
#outputs of functions of the package "MMIX" in a suitable way


print.MMIXclass<-function(x, ...){
  print(x[[1]])
}

summary.MMIXclass<-function(object, ...){
  print.default(object[[2]])
}

plot.MMIXclass<-function(x, ...){
  barplot(height=sort(x[[3]]),names.arg=names(coef)[rank(x[[3]])],
  space=1,width=0.1,cex.names=0.8,horiz=TRUE,las=2)
  title(main=rownames(x[[3]]))
}
