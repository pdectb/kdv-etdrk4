
# helper functions
sech = (x) -> 2 / (exp(x) + exp(-x))
soliton = (A, x1, x) -> 3*A.pow(2) * (sech(.5*A*(x+x1))).pow(2)

# Set up grid and two-soliton initial data:
N = 256
x = 2*pi/N * linspace(-N/2, N/2-1, N)
u0 = soliton(sqrt(800/3), 1, x)+soliton(sqrt(200/3), 0, x)

# KdV step
etdrk4 = new $blab.Etdrk4
    N: 256
    h: 4e-5
    M: 64 # no. pts for complex means
    dispersion: (z) -> j*z.pow(3)

lineChart = new $blab.LineChart
    id: "solitons"
    xLabel: "x"
    yLabel: "u"
    xLim: [-pi, pi]
    yLim: [0, 1000]
    xTicks: 7
    yTicks: 5
    click: (x, y) -> initSoliton(x, y) 

# Animation

animateId = null

animate = (snapshotFunction, numSnapshots, delay=10) ->
    stopAnimation()
    n = 1
    snapshot = ->
        snapshotFunction()
        n++
        stopAnimation(animateId) if n>numSnapshots
    animateId = setInterval (-> snapshot()), delay

stopAnimation = (animateId) ->
    clearTimeout animateId if animateId
    animateId = null

$("#kdv-stop-button").on "click", -> stopAnimation(animateId)

v = etdrk4.fft u0

snapshot = ->
    {u, v} = etdrk4.computeUV(v)        
    lineChart.plot(x, u)

animate snapshot, 200, 10
    
