using integrali

int = Interval(0, 5) 
f = x -> sin(x)/x

compute_integral_gl(f, int)

find_best_split(f, int)