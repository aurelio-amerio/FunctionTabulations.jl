using FunctionTabulations
using Documenter

DocMeta.setdocmeta!(FunctionTabulations, :DocTestSetup, :(using FunctionTabulations); recursive=true)

makedocs(;
    modules=[FunctionTabulations],
    authors="Aurelio Amerio <aurelio.amerio@edu.unito.it> and contributors",
    repo="https://github.com/aurelio-amerio/FunctionTabulations.jl/blob/{commit}{path}#{line}",
    sitename="FunctionTabulations.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "List of functions" => "functions.md",
    ],
)
