using DrWatson
@quickactivate "SolarSailingQLawProject"
include(srcdir("Includes.jl"))

function main()
    x = -10:.01:10

    println("test: abs vs smoothed abs\n")
    println("btime for abs")
    @btime absx = [abs(x[i]) for i in eachindex(x)]

    println("btime for abs_smooth")
    @btime absx_smooth = [abs_smooth(x[i]) for i in eachindex(x)]

    println("btime for FD on abs")
    @btime [ForwardDiff.derivative(abs, x[i]) for i in eachindex(x)];

    println("btime for FD on abs_smooth")
    @btime [ForwardDiff.derivative(abs_smooth, x[i]) for i in eachindex(x)]

    println("end test\n")

end
main()