using DrWatson
@quickactivate "SolarSailingQLawProject"
include(srcdir("Includes.jl"))
import GLMakie as GM
mu = 398600.4418

# read matlab data
readMatlab = true
if readMatlab
    mlData = readdlm(joinpath(@__DIR__, "data_test\\MatlabSimData.txt"), '\t', '\n'; header=false)
end


# Read and assign gmat
gmatData = readdlm(joinpath(@__DIR__, "data_test\\ReportFile1.txt"), ' ', '\n'; header=false)
gmatCart = view(gmatData, :,2:4)
gmatTime = (gmatData[:,1].-gmatData[1,1])*86400 # convert from mjd to seconds from epoch
gmatKep = view(gmatData, :,5:10)

# Read and assign julia
juliaData = readdlm(joinpath(@__DIR__, "data_test\\ReportFile2.txt"), '\t', '\n'; header=false)
juliaTime = juliaData[:,1] .- juliaData[1,1] # convert to seconds from epoch

# Interpolate julia simulated data to each time in the gmat data
sf = juliaTime[end] - juliaTime[1] # scaling factor
ts = juliaTime./sf # scaled times to make spline
states = juliaData[:,2:4] # positions to make spline with
dx0dt = juliaData[1, 5:7]
dxfdt = juliaData[end, 5:7]
dx0dts = sf.*dx0dt
dxfdts = sf.*dxfdt
spl = CubicSpline(ts, states, dx0dts, dxfdts) # creating the spline

# For each time present in the GMAT trjectory [s from epoch], interpolate the julia propagation to that point
interpJuliaData = Matrix{Float64}(undef, size(gmatData, 1), 6)
juliaKep = Matrix{Float64}(undef, size(gmatData, 1), 6)
juliaCart = Matrix{Float64}(undef, size(gmatData, 1), 3)
for idx in eachindex(gmatTime)
    t = gmatTime[idx]
    gts = t/sf 
    rtmp = interpolate(spl, gts)

    drdts = getPositionPartials(spl, gts)
    dtdts = 1/sf
    vtmp = drdts*dtdts

    interpJuliaData[idx, :] = [rtmp[1], rtmp[2], rtmp[3], vtmp[1], vtmp[2], vtmp[3]];
    juliaKep[idx, :] .= rv2coe(rtmp, vtmp, mu)
    juliaCart[idx, :] .= rtmp

end

# Plot Difference
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

xdiff = juliaCart[:,1].-gmatCart[:,1]
ydiff = juliaCart[:,2].-gmatCart[:,2]
zdiff = juliaCart[:,3].-gmatCart[:,3]
GM.lines!(axX, gmatTime, xdiff)
GM.lines!(axY, gmatTime, ydiff)
GM.lines!(axZ, gmatTime, zdiff)

# Plot orbital elements
fig2 = GM.Figure()
axa = GM.Axis(
    fig2[1,1],
    xlabel="Time[s] from epoch",
    ylabel="Semi-major axis [km]"
)
axe = GM.Axis(
    fig2[2,1],
    xlabel="Time[s] from epoch",
    ylabel="Eccentricity"
)
axi = GM.Axis(
    fig2[3,1],
    xlabel="Time[s] from epoch",
    ylabel="Inclination [deg]"
)
axape = GM.Axis(
    fig2[1,2],
    xlabel="Time[s] from epoch",
    ylabel="Arg. Perigee [deg]"
)
axran = GM.Axis(
    fig2[2,2],
    xlabel="Time[s] from epoch",
    ylabel="RAAN [deg]"
)
axtru = GM.Axis(
    fig2[3,2],
    xlabel="Time[s] from epoch",
    ylabel="True anomaly [deg]"
)

GM.lines!(axa, gmatTime, juliaKep[:,1]-gmatKep[:,1])
GM.lines!(axe, gmatTime, juliaKep[:,2]-gmatKep[:,2])
GM.lines!(axi, gmatTime, juliaKep[:,3]*180/pi-gmatKep[:,3])
GM.lines!(axape, gmatTime, juliaKep[:,4]*180/pi-gmatKep[:,4])
GM.lines!(axran, gmatTime, juliaKep[:,5]*180/pi-gmatKep[:,5])
GM.lines!(axtru, gmatTime, juliaKep[:,6]*180/pi-gmatKep[:,6])

# If we are running the matlab data, plot that 
if readMatlab
    fig3 = GM.Figure()
    axX = GM.Axis(
        fig3[1,1],
        xlabel="Time from Epoch [s]",
        ylabel="X Difference [km]",
        title="Difference (Julia-MATLAB)"
    )
    axY = GM.Axis(
        fig3[2,1],
        xlabel="Time from Epoch [s]",
        ylabel="Y Difference [km]"
    )
    axZ = GM.Axis(
        fig3[3,1],
        xlabel="Time from Epoch [s]",
        ylabel="Z Difference [km]"
    )
    gmDataLong = readdlm(joinpath(@__DIR__, "data_test\\ReportFile1_long.txt"), ' ', '\n'; header=false)

    xd = juliaData[:,2].-mlData[:,1]
    yd = juliaData[:,3].-mlData[:,2]
    zd = juliaData[:,4].-mlData[:,3]
    GM.lines!(axX, juliaTime, xd)
    GM.lines!(axY, juliaTime, yd)
    GM.lines!(axZ, juliaTime, zd)

end
fig3