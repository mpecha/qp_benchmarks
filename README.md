# QP benchmarks
## testfiles
### ex1
Displacement of a string with fixed ends subject to the lower bound. Solved as finite differences discretization of
```
-u''(x) = -15,  x in [0,1]
u(0) = u(1) = 0
s.t. u(x) > sin(4*pi*x -pi/6)/2 -2
```

### ex2
Displacement of a string with fixed ends subject to the lower bound on the first half of the domain. Solved as finite differences discretization of
```
-u''(x) = -15,  x in [0,1]
u(0) = u(1) = 0
s.t. u(x) > sin(4*pi*x -pi/6)/2 -2, x in [0,1/2]
```

### jbearing2
Pressure distribution in a journal bearing from MINPACK-2 (epsilon = 0.1, b = 10). See http://ftp.mcs.anl.gov/pub/tech_reports/reports/P153.pdf, section 4.2, page 33.

The problem has only a lower bound. The provided upper bound is set to constant vector = 1000 and can be ignored.

