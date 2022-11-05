using uCDL
using Documenter

DocMeta.setdocmeta!(uCDL, :DocTestSetup, :(using uCDL); recursive=true)

makedocs(;
    modules=[uCDL],
    authors="Shane Kuei-Hsien Chu (skchu@wustl.edu)",
    repo="https://github.com/kchu25/uCDL.jl/blob/{commit}{path}#{line}",
    sitename="uCDL.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kchu25.github.io/uCDL.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kchu25/uCDL.jl",
    devbranch="main",
)
