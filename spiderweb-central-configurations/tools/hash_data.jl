#!/usr/bin/env julia
# Write data/hashes.json: the SHA-256 of every frozen file in data/.
#
# Run from the package root, in any environment that has JSON, e.g.
#   julia --project=checker tools/hash_data.jl
#
# Hashes every file in data/ except hashes.json and dotfiles, and keeps any
# "url" already recorded. A file listed with a url but absent locally (kept out
# of git, say) is left as it is. Any other tool writing the same JSON will do.
using SHA, JSON

root = length(ARGS) ≥ 1 ? ARGS[1] : "."
dir  = joinpath(root, "data")
out  = joinpath(dir, "hashes.json")

old = isfile(out) ? JSON.parsefile(out) : Dict{String,Any}()
new = Dict{String,Any}()

for f in sort(readdir(dir))
    p = joinpath(dir, f)
    (isfile(p) && f != "hashes.json" && !startswith(f, ".")) || continue
    e = Dict{String,Any}("sha256" => bytes2hex(open(sha256, p)))
    haskey(old, f) && haskey(old[f], "url") && (e["url"] = old[f]["url"])
    new[f] = e
end
for (f, e) in old                      # listed with a url, not present locally: keep
    haskey(new, f) || !haskey(e, "url") || (new[f] = e)
end

isempty(new) && error("no frozen files in $dir")
open(out, "w") do io
    JSON.print(io, new, 2)
    println(io)
end
for f in sort(collect(keys(new)))
    println(rpad(f, 24), new[f]["sha256"], haskey(new[f], "url") ? "  (url)" : "")
end
