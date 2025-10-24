import Pkg
Pkg.add("Documenter")
using Documenter, LOShadows, SparseArrays

push!(LOAD_PATH, "../src/")

DocMeta.setdocmeta!(BosonSampling, :DocTestSetup, :(using MyPackage); recursive=true)

makedocs(
	sitename = "LOShadows.jl",
	remotes = nothing,
	modules = [LOShadows],
	authors := "Hugo Thomas",
	format = Documenter.HTML(prettyurls=false, sidebar_sitename=false),
	pages = [
		"Usage" => [
			"usage/projectors.md",
			"usage/shadows.md"
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
