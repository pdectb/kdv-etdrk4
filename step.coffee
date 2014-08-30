# Original code:
# kdv.m - Solve KdV equation by Fourier spectral/ETDRK4 scheme
#         A.-K. Kassam and L. N. Trefethen 4/03
# http://people.maths.ox.ac.uk/trefethen/publication/PDF/2005_111.pdf

class $blab.Etdrk4

    # helper functions
    zeros = (x) -> nm.zeros(1, x.length)[0]  # Zero vector same length as x
    r2c = (x) -> complex x, zeros(x)  # Real to complex (for fft input)
    ifftRe = (x) -> x.ifft().x  # Real part of ifft 
    
    constructor: (@spec) ->

        N = @spec.N
        h = @spec.h
        M = @spec.M
        dispersion = @spec.dispersion   
                
        # Precompute ETDRK4 scalar quantities (Kassam-Trefethen):

        k = linspace(0, N/2-1, N/2)
            .concat([0])
            .concat(linspace(-N/2+1, -1, N/2-1)) # wavenumbers

        L = dispersion(k)
        @E = exp(h*L) 
        @E2 = exp(h*L/2)

        r = exp(j*2*pi*((linspace 1, M, M)-0.5)/M) # roots of unity

        LR = complex(@cross(h*L.x, r.x), @cross(h*L.y, r.y))

        @Q  = h*@mean((exp(LR/2)-1)/LR)
        @f1 = h*@mean((-4-LR+exp(LR)*(4-3*LR+LR.pow(2)))/LR.pow(3))
        @f2 = h*@mean((4+2*LR+exp(LR)*(-4+2*LR))/LR.pow(3))
        @f3 = h*@mean((-4-3*LR-LR.pow(2)+exp(LR)*(4-LR))/LR.pow(3))
        @g = -0.5*j*k

    cross: (A, B) -> (a+b for b in B for a in A)

    mean: (A) -> 
        w = A.size()[1]
        R = (z.sum() for z in A.x)/w
        I = (z.sum() for z in A.y)/w
        complex R, I

    fft: (x) -> r2c(x).fft()

    Nx: (x) -> @g*@fft(ifftRe(x).pow(2))

    computeUV: (v) ->
        Nv = @Nx(v)
        a = @E2*v+@Q*Nv
        Na = @Nx(a)
        Nb = @Nx(@E2*v+@Q*Na)
        Nc = @Nx(@E2*a+@Q*(2*Nb-Nv))
        v = @E*v+(Nv*@f1+(Na+Nb)*@f2+Nc*@f3)
        u = ifftRe(v)
        {u, v}
