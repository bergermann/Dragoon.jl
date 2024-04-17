
using Distributed, JLD2

@everywhere using Dragoon


freqs = genFreqs(22.025e9,50e6; length=10)
freqsplot = genFreqs(22.025e9,150e6; length=1000)


n = 20
d0 = findpeak(22.025e9,n)
p0 = dist2pos(ones(n)*d0);

b = AnalyticalBooster(d0); b.tand = 0.;


sigx = [1e-6,10e-6,50e-6,100e-6];
N = 1000; Term = zeros(Float64,N,4,length(sigx));

for s in eachindex(sigx)
    print("\nσx = $(sigx[s]), i: 0")

    for i in 1:N
        (i%(N/10) == 0) && print(",",i)

        move(b,p0+randn(n)*sigx[s]; additive=false)

        hist = initHist(b,100,freqs,ObjAnalytical)
        b.summeddistance = 0.

        trace, term = nelderMead(b,hist,freqs,
            1.,1+2/b.ndisk,0.75-1/(2*b.ndisk),1-1/(b.ndisk),1e-6,
            ObjAnalytical,
            # InitSimplexRegular(5e-5),
            InitSimplexCoord(1e-5),
            DefaultSimplexSampler,
            UnstuckNew(InitSimplexRegular(5e-5),true,-10000);
            maxiter=2000,
            showtrace=false,
            unstuckisiter=true,
            resettimer=true,
            returntimes=true)

        Term[i,1,s] = term[1]
        Term[i,2,s] = term[2][1].value
        Term[i,3,s] = term[2][2]
        Term[i,4,s] = term[2][3].value
    end
end; plotTerm(Term,sigx; maxt=10000,maxdist=10000);

p1, p2 = plotTerm(Term,sigx; maxt=1000,maxdist=500,
    annpos=[(0.95,0.2),(0.95,0.2),(0.95,0.2),(0.95,0.6)])
savefig(p1,"20discs_nm_time.pdf")
savefig(p2,"20discs_nm_dist.pdf")



# ==============================================================================
# simulated annealing



N = 1000; Term = zeros(Float64,N,4,length(sigx));

T = collect(range(1; length=2001,stop=0));
for s in eachindex(sigx)
    print("\nσx = $(sigx[s]), i: 0")

    for i in 1:N
        (i%(N/10) == 0) && print(",",i)

        move(b,p0+randn(n)*sigx[s]; additive=false)

        hist = initHist(b,100,freqs,ObjAnalytical)
        b.summeddistance = 0.

        trace, term = simulatedAnnealing(b,hist,freqs,T,25e-6,
            ObjAnalytical,
            UnstuckDont;
            maxiter=length(T),
            nreset=50,
            nresetterm=2,
            showtrace=false,
            unstuckisiter=true,
            traceevery=1,
            resettimer=true,
            returntimes=true);

        Term[i,1,s] = term[1]
        Term[i,2,s] = term[2][1].value
        Term[i,3,s] = term[2][2]
        Term[i,4,s] = term[2][3].value
    end
end; plotTerm(Term,sigx);


p1, p2 = plotTerm1(Term,sigx)
savefig(p1,"20discs_sa_time.pdf")
savefig(p2,"20discs_sa_dist.pdf")



# ==============================================================================
# steepest descent



N = 100; Term = zeros(Float64,N,4,length(sigx));

for s in eachindex(sigx)
    print("\nσx = $(sigx[s]), i: 0")

    for i in 1:N
        (i%(N/10) == 0) && print(",",i)

        move(b,p0+randn(n)*sigx[s]; additive=false)

        hist = initHist(b,100,freqs,ObjAnalytical)
        b.summeddistance = 0.

        trace, term = linesearch(b,hist,freqs,10e-6,
            ObjAnalytical,
            SolverSteep,
            Derivator1(1e-6,"double"),
            StepNorm("unit"),
            SearchExtendedSteps(10),
            UnstuckRandom(1e-5,-13000);
            ϵgrad=0.,maxiter=Int(1000),showtrace=false,
            resettimer=true,returntimes=true);        

        Term[i,1,s] = term[1]
        Term[i,2,s] = term[2][1].value
        Term[i,3,s] = term[2][2]
        Term[i,4,s] = term[2][3].value
    end
end; plotTerm(Term,sigx);

bu = copy(Term)

p1, p2 = plotTerm(Term,sigx; annpos=[(0.95,0.2),(0.95,0.2),(0.35,0.2),(0.35,0.2)])
savefig(p1,"20discs_sd_time.pdf")
savefig(p2,"20discs_sd_dist.pdf")





# ==============================================================================
# hybrid method



N = 100; Term = zeros(Float64,N,4,length(sigx));

for s in eachindex(sigx)
    print("\nσx = $(sigx[s]), i: 0")

    for i in 1:N
        (i%(N/10) == 0) && print(",",i)

        move(b,p0+randn(n)*sigx[s]; additive=false)

        hist = initHist(b,100,freqs,ObjAnalytical)
        b.summeddistance = 0.

        trace, term = linesearch(b,hist,freqs,1e-6,
            ObjAnalytical,
            SolverHybrid("inv",0,1e-6,1),
            Dragoon.Derivator2_(1e-6,1e-6,"double"),
            StepNorm("unit"),
            SearchExtendedSteps(200),
            # UnstuckDont;
            UnstuckRandom(10e-6,-13000);
            ϵgrad=0.,maxiter=Int(200),showtrace=false,
            resettimer=true,returntimes=true);        

        Term[i,1,s] = term[1]
        Term[i,2,s] = term[2][1].value
        Term[i,3,s] = term[2][2]
        Term[i,4,s] = term[2][3].value
    end
end; plotTerm(Term,sigx);


p1, p2 = plotTerm(Term,sigx; annpos=[(0.95,0.2),(0.95,0.2),(0.35,0.2),(0.35,0.2)])
savefig(p1,"20discs_hy_time.pdf")
savefig(p2,"20discs_hy_dist.pdf")




















































N = 100; B = zeros(Float64,N,length(freqsplot),length(sigx));

for s in eachindex(sigx)
    print("\nσx = $(sigx[s]), i: 0")

    for i in 1:N
        (i%(N/10) == 0) && print(",",i)

        move(b,dist2pos(pos2dist(p0)+randn(n)*sigx[s]); additive=false)

        B[i,:,s] = getBoost1d(b,freqsplot)
    end
end;

move(b,p0; additive=false)
B0 = getBoost1d(b,freqsplot)

p1 = plot(freqsplot/1e9,[B0,B0,B0,B0]/1e3; legend=false,linesize=5,color=:red,
    layout=grid(2,2))

ylims!(p1,(-10,510)); xticks!(p1,[22.0,22.05])

for i in 1:4
    plot!(p1[i],freqsplot/1e9,B[:,:,i]'/1e3; c=:black,linesize=0.1,alpha=0.3)
end;

plot!(p1,freqsplot/1e9,[B0,B0,B0,B0]/1e3; linesize=5,color=:red)

plot!(p1[1]; xformatter=:none,bottom_margin=(-3.5,:mm),right_margin=(0.,:mm))
plot!(p1[2]; xformatter=:none, yformatter=:none,bottom_margin=(-3.5,:mm),
    left_margin=(-3.5,:mm))
plot!(p1[3]; top_margin=(0.,:mm),xlabel=" ",ylabel=" ")
plot!(p1[4]; yformatter=:none,left_margin=(-3.5,:mm))

for s in eachindex(sigx)
    annotate!(p1[s],21.95,400,
        (
            L"\sigma_x="*"$(round(Int,sigx[s]/1e-6)) µm\n",
            8,:black,:left
        ))
end

annotate!(p1[3],22.11,-100,(L"f"*" [GHz]",12,:black,:center))
annotate!(p1[3],21.9176,500,(L"f_{\beta,F}\times 10^3",10,90.,:center,:black))

display(p1)