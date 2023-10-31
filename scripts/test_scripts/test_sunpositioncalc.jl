using DrWatson
@quickactivate "SolarSailingQLawProject"

include(srcdir("Includes.jl"))
epoch = "2018-06-06T20:45:00"
et = utc2et(epoch)

tspan = 1:24*3600:3*360*24*3600
k=1
dist = zeros(length(collect(tspan)))
for i in tspan
    (~, dist[k]) = sunToEarth(et+i)
    global k+=1
end

plot(collect(tspan)./(24*3600), dist)
xlabel!("Time(days)"); ylabel!("Distance(km)")
title!("Distance to sun (Epoch: $epoch)")
