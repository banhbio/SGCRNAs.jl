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
	"Reference" => "reference.md",
    ],
)

deploydocs(;
    repo="github.com/C37H41N2O6/SGCRNAs.jl.git",
    devbranch="main",
    versions = [
	"latest" => "main",
	"v1.0.0" => "v1.0.0",
    ]
)
