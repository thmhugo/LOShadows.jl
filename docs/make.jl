import Pkg
Pkg.add("Documenter")
using Documenter, LOShadows, SparseArrays

push!(LOAD_PATH, "../src/")

makedocs(
	sitename = "LOShadows.jl",
	remotes = nothing,
	modules = [LOShadows],
	source = "src",
	authors = "Hugo Thomas",
	format = Documenter.HTML(prettyurls=false, sidebar_sitename=false),
	pages = [
		"About" => "index.md",
		"Usage" => [
			"usage/shadows.md",
			"usage/projectors.md"
		],
		"API" => [
			"functions" => [
				"functions/Channels.md",
				"functions/FockOperators.md",
				"functions/MatrixOperations.md",
			],
			"applications" => [
				"applications/Binning.md",
				"applications/Invariants.md"
			],
			"types" => [
				"types/DensityMatrix.md"
			],
		],
	],
)


deploydocs(
	repo = "github.com/thmhugo/LOShadows.jl.git",
	target = "build",
    branch = "gh-pages",
	versions = ["stable" => "v^", "v#.#"],
)
