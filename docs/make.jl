using Documenter, DocumenterVitepress, OpenSOAP

makedocs(
    format = DocumenterVitepress.MarkdownVitepress(
        repo="github.com/thanasipantazides/opensoap",
        # mathengine = Documenter.KaTeX(
        #     Dict(
        #         :macros => Dict(
        #             "\\mbs" => "\\boldsymbol{#1}",
        #             "\\ket" => "\\left|#1\\right\\rangle" # macro with 1 argument
        #         )
        #     )
        # )
    ),
    sitename = "OpenSOAP Reference",
    pages = [
        "Home" => "index.md",
        "Quick start" => "quickstart.md",
        "Users Guide" => ["overview.md"],
        "Theory" => ["conops.md", "forces.md", "moments.md", "onboard.md", "targets.md"],
    ],
)
