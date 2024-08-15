
using Pkg; Pkg.update()
using Dragoon, BoostFractor
using HDF5, DataFrames

include("tools/tools.jl");

data = prepareDataAll1d(getPath(),0.4; filterin="REF");
sortData!(data)

ref1 = ref1d([
    0.0072448589080291290,
    0.0070711997431429300,
    0.0071246144242139010,
    0.0073342065660182530,
    0.0072327454408213510,
    0.0072511235941711350,
    0.0071632103241744970,
    0.0072023665514860190,
    0.0073048489540810120,
    0.0071862085372421955,
    0.0070924346714088270,
    0.0073231298843525620,
    0.0072555088582221690,
    0.0072399407925294880,
    0.0072071841247186580,
    0.0072590448084468580,
    0.0072497272553200550,
    0.0070498810897688120,
    0.0069997371085350490,
    0.0071354746829380430,
],data.freqs; eps=data.s.eps,tand=data.s.tand)

ref2 = ref1d([
    0.0069864996648397000,
    0.0072560130634651010,
    0.0072587766881647050,
    0.0072180129444167985,
    0.0071810463115898690,
    0.0072643048247736110,
    0.0072005219024528160,
    0.0072281469176733665,
    0.0072039644723843060,
    0.0072056913829301470,
    0.0072293906115958880,
    0.0071775904582249700,
    0.0072014504364220360,
    0.0073195415726201870,
    0.0071448147947109080,
    0.0073325781037591470,
    0.0071798833031334620,
    0.0073681846490424320,
    0.0066949319976553630,
    0.0074786250621039780,
]
,data.freqs; eps=data.s.eps,tand=data.s.tand)

ref3 = ref1d([
    0.0069424515103718450,
    0.0070364724794998870,
    0.0078732484028712010,
    0.0070004748835278420,
    0.0073614019954116670,
    0.0071549044978148910,
    0.0070875765947395550,
    0.0074426096458853840,
    0.0073301995131451340,
    0.0074220155045692520,
    0.0071730390572375410,
    0.0070696702277422010,
    0.0073113174166905695,
    0.0065711615775509440,
    0.0095807997342139400,
    0.0048028172495510070,
    0.0079482239556944450,
    0.0064862237924530890,
    0.0076375922093132050,
    0.0068707736863909690,
],data.freqs; eps=data.s.eps,tand=data.s.tand)

match = zeros(Int,length(data))
for i in eachindex(match)
    r1 = sum(abs.(data.ref[:,i]-ref1))
    r2 = sum(abs.(data.ref[:,i]-ref2))
    r3 = sum(abs.(data.ref[:,i]-ref3))

    match[i] = argmin([r1,r2,r3])
end
obj_ = zeros(length(data));
for i in eachindex(obj_); obj_[i] = -minimum(data.boost[:,i]); end
data = data[findall(x->x<-14_000,obj_)]

initdist = findpeak1d(data.s.f0,20)
d0 = initdist*ones(data.s.ndisk); p0 = dist2pos(d0);
freqsplot = genFreqs(22.025e9,150e6; n=1000);


histogram(data.obj; xlabel="Reference Match ∑|R-R_ref|",title="Distribution of Converged States (20 Discs, ε=24)",
    label=false,)

showDist(data,1000; xlabel="Disc Index", ylabel="d_i [mm]")
# hline!([initdist]*1e3)


showDistribution(data,d0)




out = findOutliers(data,6e-3; showdistribution=true); length(out)
showDist(out)
showFields(out,data.freqs,freqsplot)


obj1 = obj_[findall(isequal(1),match)]
data1 = data[findall(isequal(1),match)]

obj2 = obj_[findall(isequal(2),match)]
data2 = data[findall(isequal(2),match)]

obj3 = obj_[findall(isequal(3),match)]
data3 = data[findall(isequal(3),match)]