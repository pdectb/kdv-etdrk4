# Soliton.
sech = (x) -> 2 / (exp(x) + exp(-x))
soliton = (A, x1, x) -> A * (sech(sqrt(A/12)*(x+x1))).pow(2)

# Grid and initial condition.
N = 256
x = 2*pi/N * linspace(-N/2, N/2-1, N) #;
u0 = soliton(800, 1, x)+soliton(200, 0, x) #;

plot x, u0,
    xlabel: "x"
    ylabel: "u0(x)"
    height: 160
    series: 
        shadowSize: 0
        color: "black"
        lines: lineWidth: 1

# KdV step (imported from step.coffee)
etdrk4 = new $blab.Etdrk4
    N: 256
    h: 4e-5
    M: 64 # no. pts for complex means
    dispersion: (z) -> j*z.pow(3) # KdV dispersion.

# Evolve initial condition
v = etdrk4.fft u0 #;
for n in [1..200] #;
    {u, v} = etdrk4.computeUV(v)        

plot x, u,
    xlabel: "x"
    ylabel: "u(x)"
    height: 160
    series: 
        shadowSize: 0
        color: "black"
        lines: lineWidth: 1

