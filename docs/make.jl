using UnfoldCDL
using Documenter

DocMeta.setdocmeta!(UnfoldCDL, :DocTestSetup, :(using UnfoldCDL); recursive=true)

makedocs(;
    modules=[UnfoldCDL],
    authors="Shane Kuei-Hsien Chu (skchu@wustl.edu)",
    repo="https://github.com/kchu25/UnfoldCDL.jl/blob/{commit}{path}#{line}",
    sitename="UnfoldCDL.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kchu25.github.io/UnfoldCDL.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kchu25/UnfoldCDL.jl",
    devbranch="main",
)
