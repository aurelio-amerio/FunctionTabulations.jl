using Tabulations
using Documenter

DocMeta.setdocmeta!(Tabulations, :DocTestSetup, :(using Tabulations); recursive=true)

makedocs(;
    modules=[Tabulations],
    authors="Aurelio Amerio <aurelio.amerio@edu.unito.it> and contributors",
    repo="https://github.com/aurelio-amerio/Tabulations.jl/blob/{commit}{path}#{line}",
    sitename="Tabulations.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "List of functions" => "functions.md",
    ],
)
