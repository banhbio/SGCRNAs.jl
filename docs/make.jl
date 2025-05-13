using SGCRNAs
using Documenter

DocMeta.setdocmeta!(SGCRNAs, :DocTestSetup, :(using SGCRNAs); recursive=true)

makedocs(;
    modules=[SGCRNAs],
    authors="Tatsunori Osone",
    sitename="SGCRNAs.jl",
    format=Documenter.HTML(;
        canonical="https://C37H41N2O6.github.io/SGCRNAs.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
	"Quick Start" => "quickstart.md",
	"Tutorial" => "tutorial.md",
    ],
)

deploydocs(;
    repo="github.com/C37H41N2O6/SGCRNAs.jl",
    devbranch="main",
)
