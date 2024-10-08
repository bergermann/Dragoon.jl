
using Pkg; Pkg.update()


using Dragoon
using Plots, Plots.Measures
using Dates
using JLD2

include(joinpath(pwd(),"examples\\analysis\\tools\\tools.jl"));

@load "examples\\full_20_24.0_6.0e-5.jld2"
# @load "examples\\full_20_24.0.jld2"

f = 25

freqs = genFreqs(f*1e9+25e6,100e6; n=100);
freqsplot = genFreqs(f*1e9+25e6,150e6; n=1000);

booster = AnalyticalBooster(P0[f]; tand=6e-5)
# booster = AnalyticalBooster(P0[f]; tand=0)
obj = ObjAnalytical
hist = initHist(booster,2*(booster.ndisk^2),freqs,obj);

plot(freqsplot/1e9,getBoost1d(booster,freqsplot))


ref0 = getRef1d(booster,freqs);
# ref0 = [
#   0.9240607675676515 - 0.38224559885263165im,
#   0.7380980822551304 - 0.674693427395952im,
#   0.462948064022487 - 0.8863854071552547im,
#   0.14786245633322903 - 0.9890079342488161im,
#   -0.15586000240126052 - 0.9877791553031501im,
#   -0.41490614421034666 - 0.9098642159665684im,
#   -0.617404400711767 - 0.7866459216074747im,
#   -0.7658156811116661 - 0.6430601391498648im,
#   -0.8687939779605591 - 0.49517373098685186im,
#   -0.9360561794775265 - 0.3518505774642497im,
#   -0.9761192037376331 - 0.21723558662114203im,
#   -0.9956895298891584 - 0.09274890872183986im,
#   -0.9997664170139217 + 0.02161276037728511im,
#   -0.9919648735330898 + 0.12651359482877816im,
#   -0.9748503346892108 + 0.2228605504711954im,
#   -0.9502174864493367 + 0.3115874330679862im,
#   -0.9193032370696224 + 0.393549943861387im,
#   -0.8829429833546631 + 0.46948023189991367im,
#   -0.8416837189379979 + 0.5399708485415708im,
#   -0.7958656259972035 + 0.6054732903737947im,
#   -0.7456820601606868 + 0.6663019324259252im,
#   -0.6912243153190141 + 0.7226402603728616im,
#   -0.6325172274726643 + 0.7745462910312649im,
#   -0.5695479277801538 + 0.8219581242139662im,
#   -0.5022919227600631 + 0.8646981116725166im,
#   -0.4307376970992456 + 0.9024771666350553im,
#   -0.35491058222663224 + 0.9349002506275896im,
#   -0.2748980884521972 + 0.9614733698679944im,
#   -0.19087567409544032 + 0.9816142200674294im,
#   -0.10313357313905225 + 0.9946675153495173im,
#   -0.012103106793501572 + 0.9999267547205412im,
#   0.08161873838943179 + 0.9966636250730342im,
#   0.1772493286384882 + 0.9841659796483582im,
#   0.27380786206107666 + 0.9617844117439002im,
#   0.3701142439384575 + 0.9289862466333114im,
#   0.46480279904208954 + 0.8854142296138187im,
#   0.5563540269667009 + 0.8309453632327204im,
#   0.6431457015941184 + 0.7657438256499466im,
#   0.7235225589439993 + 0.6903007364179202im,
#   0.7958805407946918 + 0.6054536850861058im,
#   0.8587587180131229 + 0.512380195008022im,
#   0.9109293743169278 + 0.4125623286323685im,
#   0.9514758477257652 + 0.307723432962967im,
#   0.9798484268083656 + 0.1997424854185147im,
#   0.9958914606454912 + 0.09055494803379414im,
#   0.9998389086057622 - 0.017948728033846926im,
#   0.992280221622198 - 0.12401597388001956im,
#   0.9741024265108893 - 0.2261071928657065im,
#   0.946416901645635 - 0.32294743888049116im,
#   0.9104801519376988 - 0.41355276921755146im,
# ];
plot(freqs/1e9,real.(ref0); c=:blue,label="real"); plot!(freqs/1e9,imag.(ref0); c=:red,label="imag")
plot(freqs/1e9,abs.(ref0))
obj1 = ObjRef1dSquare(ref0);
obj2 = ObjRef1dS(ref0,x->x);
obj3 = Dragoon.ObjRef1dSP(ref0,x->x);


d0 = findpeak1d((freqs[1]+freqs[end])/2,booster.ndisk;
    eps=booster.epsilon,tand=booster.tand,granularity=10_000,deviation=0.3)
move(booster,dist2pos(ones(booster.ndisk)*d0); additive=false);

plot(freqsplot/1e9,getBoost1d(booster,freqsplot)/1e3; xlabel="Frequency [GHz]",ylabel="Boostfactor β² × 10³",label="")
plot!(twinx(),freqsplot/1e9,abs.(getRef1d(booster,freqsplot)); ylabel="|S_11|",c=:red,label="")

move(booster,dist2pos(ones(booster.ndisk)*d0); additive=false);
trace = nelderMead(booster,hist,freqs,
            1.01,1+2/booster.ndisk,0.75-1/2booster.ndisk,1-1/booster.ndisk,1e-12,
            obj2,
            # ObjAnalytical,
            InitSimplexRegular(1e-4),
            DefaultSimplexSampler,
            # UnstuckNew(InitSimplexRegular(1e-4),true,5);
            UnstuckDont;
            maxiter=Int(4e2),
            showtrace=true,
            showevery=Int(1e0),
            unstuckisiter=true,);

            
plot(freqsplot/1e9,getBoost1d(booster,freqsplot))
plot(freqs/1e9,abs.(ref0)); plot!(freqsplot/1e9,abs.(getRef1d(booster,freqsplot)))



# move(booster,dist2pos(ones(booster.ndisk)*d0); additive=false);
# trace = simulatedAnnealing(booster,hist,freqs,
#             10e-6,
#             TempLinear(0.1,100_001),
#             obj2,
#             UnstuckDont;
#             maxiter=Int(100_001),
#             nreset=500,
#             showtrace=true,
#             showevery=1000,
#             unstuckisiter=true,
#             traceevery=1000);
            
plot(freqsplot/1e9,getBoost1d(booster,freqsplot))
plot(freqs/1e9,abs.(ref0)); plot!(freqsplot/1e9,abs.(getRef1d(booster,freqsplot)))
plot(1:20,pos2dist(P0[f]); seriestype=:scatter); plot!(1:20,pos2dist(booster.pos); seriestype=:scatter)




for i in eachindex(trace)
    println(i,": ",getSimplexInnerSize(trace[i].x))
end