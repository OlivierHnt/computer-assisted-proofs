#!/usr/bin/env julia
# Download the frozen files listed with a "url" in data/hashes.json that are not
# present locally, and verify the SHA-256 of every listed file.
#
# Run from the package root, in any environment that has JSON, e.g.
#   julia --project=checker tools/fetch_data.jl
using Downloads, SHA, JSON

root = length(ARGS) ≥ 1 ? ARGS[1] : "."
dir  = joinpath(root, "data")
h    = JSON.parsefile(joinpath(dir, "hashes.json"))

for f in sort(collect(keys(h)))
    p = joinpath(dir, f)
    if !isfile(p)
        haskey(h[f], "url") || error("$f is missing and has no url in hashes.json")
        println("downloading ", f)
        Downloads.download(h[f]["url"], p)
    end
    bytes2hex(open(sha256, p)) == h[f]["sha256"] || error("sha256 mismatch for $f")
    println("pass sha256 ", f)
end
