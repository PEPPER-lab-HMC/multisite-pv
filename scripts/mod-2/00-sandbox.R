## Testing sandbox
a = 1.5
b = 8
c = 0.3
x = seq(0, 1, 0.05)
y1 = a*exp(-b*x) + c

tpf <- function(x, a, b, c, d) {
  a*exp(-b*x) + c + d*(x-cp)
}

d = -0.2
cp = 0.6
y2 = c + d*(x-cp)
plot(x, y1, ylim = c(0, 1.3))
points(x, y2, col = "red")
abline(h = 0.3)
abline(v = 0.6)

tpf(0.6, a, b, c, d)
