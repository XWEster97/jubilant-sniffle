#coding=utf-8
ENV["GKS_ENCODING"] = "utf-8"
using DataFrames,Plots,JLD,IterativeSolvers,LaTeXStrings
using Plots.PlotMeasures
pyplot()
pyplot(tickfont=font("serif"), titlefont=font("serif"))

function ht(a,b)
    return hcat(a,b)
end

result=load("data/initgitrun.jld")
result1=result["MRLSBCG"]
result2=result["Diagonal"]
malloc1=[getproperty(x,:malloc) for x in result1[!,:memallocs]]
malloc2=[getproperty(x,:malloc) for x in result2[!,:memallocs]]
history_data1=[getproperty(x,:data) for x in result1[!,:history]]
history_data2=[getproperty(x,:data) for x in result2[!,:history]]
resnorm_data1=[x[:resnorm] for x in history_data1]
resnorm_data2=[x[:resnorm] for x in history_data2]
history_iters1=[getproperty(x,:iters) for x in result1[!,:history]]
history_iters2=[getproperty(x,:iters) for x in result2[!,:history]]
realloc1=[getproperty(x,:pause) for x in result1[!,:memallocs]]
realloc2=[getproperty(x,:pause) for x in result2[!,:memallocs]]


realloc=ht(realloc1,realloc2)
x1=[48*i^4 for i=5:5:size(result1,1)*5]
x2=[48*36*i^2 for i=5:5:size(result1,1)*5]
x=ht(x1,x2)
malloc=ht(malloc1,malloc2)
history_iters=ht(history_iters1,history_iters2)
P1=plot(x,ht(result1[!,:t],result2[!,:t]),ylabel=" Execution time(s)")
P2=plot(x,ht(result1[!,:bytes],result2[!,:bytes]),ylabel="Memoray allocation(bytes)")
P3=plot(x,ht(result1[!,:gctime],result2[!,:gctime]),ylabel="Garbage Collection(s)")
P4=plot(x,history_iters,ylims = (40, 200),ylabel="# of iterations")
P5=plot(x,malloc,xlabel = "# of regression coefficients",ylabel="Requested memory")
# P5=plot(x,realloc,xlabel = "# of regression coefficients",ylabel="Requested memory")
P_pro=plot(P1, P2, P3, P4, P5,layout = grid(5, 1, heights=[0.2 ,0.2, 0.2, 0.2, 0.2]),
    size=(630,900),legend = false,left_margin = 2mm,
    xtickfontsize=10,ytickfontsize=10,right_margin = 2mm)

# plot(Plots.fakedata(50, 5), w = 3)
savefig(P_pro,"data/initgitrun.eps")

pyplot(tickfont=font("serif"), titlefont=font("serif"))
result1[!,:history][1].data[:resnorm]

history_data1=[getproperty(x,:data) for x in result1[!,:history]]
history_data2=[getproperty(x,:data) for x in result2[!,:history]]
