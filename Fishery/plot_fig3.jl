# title: Stackelberg Evolutionary Game Theory - Fish application 
# author: alexander stein

#################
### libraries ###
#################

using Plots


#########################
### Define the system ###
#########################

# Parameterisation
Kmax = 10000    # for the g-function approach
r = 1
sigmaK = 0.5
sigmaH = 0.5

c = 500.0         # for the manager

# Define helping function
K(v) = Kmax*exp(-v^2/sigmaK^2)
H(v,m) = m*exp(-v^2/sigmaH^2)

# Define the G-function
G(v,u,x,m) = r*( 1 - x /(K(v)) ) - H(v,m)

# Define the quality function 
Q(m) = H(u,m)*x - c*m

# Define the ecological equilibria
xstar(u,m) = exp(-u^2/sigmaH^2 - u^2/sigmaK^2)*Kmax*(r*exp(u^2/sigmaH^2)-m)/r

# Define the evolutionary equilibria
ustar0(m) = 0
ustar1(m) = sigmaH*sqrt( log( (m*sigmaH^2+m*sigmaK^2)/(r*sigmaH^2) ) )
ustar2(m) = -sigmaH*sqrt( log( (m*sigmaH^2+m*sigmaK^2)/(r*sigmaH^2) ) )

# Define the best response
mstar(u) =  exp(u^2/sigmaH^2) * (Kmax-c*exp(u^2/sigmaH^2+u^2/sigmaK^2))*r/(2*Kmax)

############
### Main ###
############

function plot_bifurcation(mc)
    mrange_before = LinRange(0,mc,100)
    mrange_after = LinRange(mc,r,100)
    
    ustar0_vec_before = [ustar0(m) for m in mrange_before]
    ustar0_vec_after = [ustar0(m) for m in mrange_after]
    ustar1_vec = [ustar1(m) for m in mrange_after]
    ustar2_vec = [ustar2(m) for m in mrange_after]

    urange = LinRange(-0.6, +0.6, 100)
    mstar_vec = [mstar(u) for u in urange]

    plt = plot(mrange_before, ustar0_vec_before, label="equilibria u*(m)", linewidth=3.0, color=:black, linestyle=:solid)
    plot!(mrange_after, ustar0_vec_after, label=nothing, linewidth=3.0, color=:black, linestyle=:dash)
    plot!(mrange_after, ustar1_vec, label=nothing, linewidth=3.0, color=:black, linestyle=:solid)
    plot!(mrange_after, ustar2_vec, label=nothing, linewidth=3.0, color=:black, linestyle=:solid)

    plot!(mstar_vec, urange, label="best response m*(u)", linewidth=3.0, color=:orange, linestyle=:solid)

    #scatter!([0.475],[0.0],label="",marker = (:circle, 8, 0.6, :blue))
    scatter!([0.475],[0.0],label="",marker = (:circle, 8, 0.6, :blue), text = "Stackelberg                 \n Nash          \n \n \n", mode="markers+text", textposition="top center")

    plot!(xlabel="Harvesting effort (m)")
    plot!(ylabel="Catchability (u)")
    plot!(xtickfontsize=14,ytickfontsize=14,xguidefontsize=16,yguidefontsize=16)
    plot!(legend_position=(0.15,0.85))

    #annotate!(xpos, ypos, "my text", :color)
    #plot!(title="B", titleloc= :left, titlefont=26)
    #plot!(legend=nothing)

    savefig(plt, "bifurcation_fish")

    return plt
end

mc=0.5
plt_bifurcation = plot_bifurcation(mc)
plt_bifurcation

function plot_Gfunction()
    m0 = 0.25
    xstar0 = 7500
    ustar0 = 0

    m1 = 0.75
    xstar1 = 3333
    ustar1 = 0.318381

    vrange = LinRange(-0.6, +0.6, 100)

    Gv0 = [G(v,ustar0,xstar0,m0) for v in vrange]
    Gv1 = [G(v,ustar1,xstar1,m1) for v in vrange]

    plt = plot()
    plot!(vrange, Gv0, label="m=0.25", linewidth=3.0, color=:blue, linestyle=:solid)
    plot!(vrange, Gv1, label="m=0.75", linewidth=3.0, color=:red, linestyle=:solid)

    scatter!([0.0],[0.0],label="",marker = (:circle, 8, 0.6, :blue), text = "\n \n u₀*", mode="markers+text", textposition="top center")
    scatter!([0.33],[0.0],label="",marker = (:circle, 8, 0.6, :red), text = "\n u₊*", mode="markers+text", textposition="top center")
    scatter!([-0.33],[0.0],label="",marker = (:circle, 8, 0.6, :red), text = "\n u₋*", mode="markers+text", textposition="top center")


    plot!(xlabel="Catchbility (v)")
    plot!(ylabel="Growth rate (G)")
    plot!(xtickfontsize=14,ytickfontsize=14,xguidefontsize=16,yguidefontsize=16)
    #plot!(legend_position=(0.7,0.7))
    #plot!(legend=:inside)
    plot!(legend_position=(0.45,0.5))
    #plot!(title="A", titleloc= :left, titlefont=26)
    #plot_titlelocation

    savefig(plt, "GoverV_fish")
    return plt
end

plt_GoverV = plot_Gfunction()

plt = plot(plt_GoverV , plt_bifurcation, layout=grid(1,2),size=(1200,500), title=["" ""], titleloc = :left, margin=12Plots.mm)

#annotate!( 2,  2, text("A", :left, 48), subplot=2)
#annotate!( 2,  2, text("B", :left, 48), subplot=2)

annotate!((-1.85, 0.75,text("A", :left, 22)))
annotate!((-0.30, 0.75,text("B", :left, 22)))
#annotate!((0.25, 0.10 ,text("Stack/Nash", :left, 10)))

plot!(legendfont=font(12))
savefig(plt, "fish2.pdf")