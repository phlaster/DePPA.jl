using DePPA
using DePPA.Oligos
using DePPA.Alignments
using DePPA.Primers

using Documenter

DocMeta.setdocmeta!(DePPA, :DocTestSetup, :(using DePPA); recursive=true)

makedocs(;
    modules=[DePPA, DePPA.Oligos, DePPA.Alignments, DePPA.Primers],
    authors="phlaster <phlaster@users.noreply.github.com>",
    sitename="DePPA.jl",
    format=Documenter.HTML(;
        edit_link="master",
        assets=String[],
    ),
    warnonly=[:missing_docs],
    pages=[
        "Home" => "index.md",
        "API Reference" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/phlaster/DePPA.jl",
    devbranch="master",
    push_preview=true
)