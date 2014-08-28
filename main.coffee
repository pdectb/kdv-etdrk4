# Original code:
# kdv.m - Solve KdV equation by Fourier spectral/ETDRK4 scheme
#         A.-K. Kassam and L. N. Trefethen 4/03
# http://people.maths.ox.ac.uk/trefethen/publication/PDF/2005_111.pdf

# helper functions
zeros = (x) -> nm.zeros(1, x.length)[0]  # Zero vector same length as x
r2c = (x) -> complex x, zeros(x)  # Real to complex (for fft input)
fft = (x) -> r2c(x).fft()
ifftRe = (x) -> x.ifft().x  # Real part of ifft 
sech = (x) -> 2 / (exp(x) + exp(-x))
soliton = (A, x1) -> 3*A.pow(2) * (sech(.5*A*(x+x1))).pow(2)

# Set up grid and two-soliton initial data:
N = 256
x = 2*pi/N * linspace(-N/2, N/2-1, N)
u = soliton(4.5,1.5)+soliton(2.9,-1.5)

eplot x, u,
    xlabel: "x"
    ylabel: "y"




# Precompute ETDRK4 scalar quantities (Kassam-Trefethen):
h = 1e-3 # time step (>> than used to find steady state in Octave)
k = linspace(0, N/2-1, N/2)
    .concat([0])
    .concat(linspace(-N/2+1, -1, N/2-1)) # wavenumbers
L = j*k.pow(3) + r2c(1/16*(1 - 0.2*k.pow(2))) # r2c: quirk. Add j0.
E = exp(h*L) 
E2 = exp(h*L/2)
M = 64 # no. pts for complex means
r = exp(2*j*pi*((linspace 1, M, M)-0.5)/M) # roots of unity

cross = (A, B) -> (a+b for b in B for a in A)

LR = complex(cross(h*L.x, r.x), cross(h*L.y, r.y))

mean = (A) -> 
    w = A.size()[1]
    R = (z.sum() for z in A.x)/w
    I = (z.sum() for z in A.y)/w
    complex R, I

Q  = h*mean((exp(LR/2)-1)/LR)
f1 = h*mean((-4-LR+exp(LR)*(4-3*LR+LR.pow(2)))/LR.pow(3))
f2 = h*mean((4+2*LR+exp(LR)*(-4+2*LR))/LR.pow(3))
f3 = h*mean((-4-3*LR-LR.pow(2)+exp(LR)*(4-LR))/LR.pow(3))
g = -0.5*j*k

# Time-stepping by ETDRK4 formula (Cox-Matthews):

t = 0
kk = 0
v=fft(u)
Nx = (x) -> g*fft(ifftRe(x).pow(2))

y = []

simStep = (step) ->
    t = t+h
    Nv = Nx(v)
    a = E2*v+Q*Nv
    Na = Nx(a)
    Nb = Nx(E2*v+Q*Na)
    Nc = Nx(E2*a+Q*(2*Nb-Nv))
    v = E*v+(Nv*f1+(Na+Nb)*f2+Nc*f3)
    [v.x[0], v.y[0]] = [0, 0] # zero DC
    [v.x[N/2+1], v.y[N/2+1]] = [0, 0] # squash instabilty
    y.push ifftRe(v) if step % 50 is 0
    t+h/2 < .2  # Simulation run condition
    
fig = figure
    xlabel: "x"
    ylabel: "y"
    
simEnd = -> eplot x, y, fig: fig  #;



new $blab.Simulation {step: simStep, end: simEnd}

#!end (coffee)

