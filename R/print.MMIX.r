#A "MMIXclass" class has been created, in order to print, summary, and plot the
#outputs of functions of the package "MMIX" in a suitable way


print.MMIXclass<-function(x, ...){
  print(x[[1]],...)
}

summary.MMIXclass<-function(object, ...){
  print.default(object[[2]],...)
}

plot.MMIXclass<-function(x, ...){
  barplot(height=sort(x[[3]]),names.arg=names(coef)[rank(x[[3]])],horiz=TRUE,...)
  title(main=rownames(x[[3]]))
}
