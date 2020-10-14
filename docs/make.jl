using Documenter, CubeSatSim

makedocs(;
    modules=[CubeSatSim],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/serenity4/CubeSatSim.jl/blob/{commit}{path}#L{line}",
    sitename="CubeSatSim.jl",
    authors="Cédric Belmant",
    assets=String[],
)
