x <- data_frame(x = 1:10, y = x)

ggplot(x,aes(x,y)) + geom_point()
