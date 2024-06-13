using Documenter, openSOAR

makedocs(
    format = Documenter.LaTeX(),
    sitename="openSOAR Documentation", 
    modules=[openSOAR],
    authors="Athanasios Pantazides",
    pages=Any[
        "Home"=>"index.md",
        "Installation"=>"install.md",
        "Background"=>Any[
            "theory/attitude.md",
            "theory/numerical.md",
            "theory/viewing.md"
        ],
        "API"=>"api.md"
    ]
)