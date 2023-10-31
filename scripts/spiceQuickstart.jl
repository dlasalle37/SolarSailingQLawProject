using DrWatson
@quickactivate "SolarSailingQLawProject"

using SPICE
using Downloads: download

alreadyDownloaded = true
if !alreadyDownloaded
    const LSK = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls"
    const SPK = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp"

    # Download kernels
    download(LSK, "naif0012.tls")
    download(SPK, "de440.bsp")  
end

# Load leap seconds kernel
furnsh("naif0012.tls")

# Convert the calendar date to ephemeris seconds past J2000
et = utc2et("2018-02-06T20:45:00")

# Load a planetary ephemeris kernel
furnsh("de440.bsp")

# Get the position of Mars at `et` w.r.t. Earth
spkpos("mars_barycenter", et, "J2000", "none", "earth")

# get position of earth wrt sun (point from sun to earth for sunlight)
(a, c) = spkpos("earth", et, "J2000", "none", "SSB")