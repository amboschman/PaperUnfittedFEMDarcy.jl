using PaperUnfittedFEMDarcy
using Documenter

DocMeta.setdocmeta!(PaperUnfittedFEMDarcy, :DocTestSetup, :(using PaperUnfittedFEMDarcy); recursive=true)

makedocs(;
    modules=[PaperUnfittedFEMDarcy],
    authors="Santiago Badia <santiago.badia@monash.edu>, Anne Boschman <anne.boschman@monash.edu>, Alberto F. Martin <alberto.f.martin@anu.edu.au>, Erik Nilsson <erik.nilsson@umontpellier.fr>, Ricardo Ruiz-Baier <ricardo.ruizbaier@monash.edu>, Sara Zahedi <sara.zahedi@math.kth.se>",
    sitename="PaperUnfittedFEMDarcy.jl",
    format=Documenter.HTML(;
        canonical="https://amboschman.github.io/PaperUnfittedFEMDarcy.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/amboschman/PaperUnfittedFEMDarcy.jl",
    devbranch="main",
)
