# http://stats.stackexchange.com/questions/115020/best-approach-in-r-for-interpolating-and-curve-fitting-a-tiny-dataset
# http://cran.r-project.org/web/packages/ggdendro/vignettes/ggdendro.html
# https://rpubs.com/mcgill_linguistics/63173
# http://rforpublichealth.blogspot.it/2014/12/animations-and-gifs-using-ggplot2.html
require(ggplot2)

p.dp.1 <- ggplot(data = DP,
                 aes(x = Conc,
                     y = Activity)) + 
  ylab("Activity") + 
  xlab("Concentration") +
  geom_point() + 
  theme_bw(base_size = 8) +
  geom_smooth(method = "lm", 
              aes(colour = "Linear")) +
  scale_color_manual(name = "Fits",
                     breaks = c("Linear"),
                     values = c("blue"))


require(ggplot2)

p.dp.1 <- ggplot(data = DP,
                 aes(x = Conc,
                     y = Activity)) + 
  ylab("Activity") + 
  xlab("Concentration") +
  geom_point() + 
  theme_bw(base_size = 8) +
  geom_smooth(method = "lm", 
              aes(colour = "Linear")) +
  scale_color_manual(name = "Fits",
                     breaks = c("Linear"),
                     values = c("blue"))

libary(ggplot2)
d = data.frame(DPConc= c(0, 83, 166, 416),
               DPActivity=c(100, 67.71, 6.3, 16.55))
ggplot(d, aes(x=DPConc, y=DPActivity)) + geom_point(size=3) + geom_smooth() + theme_minimal()

n <- 100
x <- seq(n)
y <- rnorm(n, 50 + 30 * x^(-0.2), 1)
Data <- data.frame(x, y)

plot(y ~ x, Data)

# fit a loess line
loess_fit <- loess(y ~ x, Data)
lines(Data$x, predict(loess_fit), col = "blue")

# fit a non-linear regression
nls_fit <- nls(y ~ a + b * x^(-c), Data, start = list(a = 80, b = 20, 
                                                      c = 0.2))
lines(Data$x, predict(nls_fit), col = "red")

library(ggplot2)
ggplot(Data, aes(x,y)) + geom_point() + geom_smooth()