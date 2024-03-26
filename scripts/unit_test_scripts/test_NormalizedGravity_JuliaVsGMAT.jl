using DrWatson
@quickactivate "SolarSailingQLawProject"
include(srcdir("Includes.jl"))
import GLMakie as GM

gmatData = readdlm(datadir("ReportFile1.txt"), ' ', '\n'; header=false)
juliaData = readdlm(datadir("ReportFile2.txt"), '\t', '\n'; header=false)

gmatCart = view(gmatData, :,2:4)
gmatTime = (gmatData[:,1].-gmatData[1,1])*86400 # convert from mjd to seconds from epoch
juliaCart = view(juliaData, :,2:4)
juliaTime = juliaData[:,1] .- juliaData[1,1] # convert to seconds from epoch

# Plot results
fig = GM.Figure()
axX = GM.Axis(
    fig[1,1],
    xlabel="Time from Epoch [s]",
    ylabel="X Difference [km]",
    title="Difference (Julia-GMAT)"
)
axY = GM.Axis(
    fig[2,1],
    xlabel="Time from Epoch [s]",
    ylabel="Y Difference [km]"
)
axZ = GM.Axis(
    fig[3,1],
    xlabel="Time from Epoch [s]",
    ylabel="Z Difference [km]"
)

GM.lines!(axX, juliaTime, juliaCart[:,1].-gmatCart[:,1])
GM.lines!(axY, juliaTime, juliaCart[:,2].-gmatCart[:,2])
GM.lines!(axZ, juliaTime, juliaCart[:,3].-gmatCart[:,3])

fig