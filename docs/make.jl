using Documenter, OpenSOAP

makedocs(
    format = Documenter.LaTeX(),
    sitename="OpenSOAP Documentation", 
    modules=[OpenSOAP],
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