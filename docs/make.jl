using Documenter, LOShadows, SparseArrays

push!(LOAD_PATH, "../src/")

makedocs(
	sitename = "LOShadows.jl",
	remotes = nothing,
	modules = [LOShadows],
	format = Documenter.HTML(prettyurls = true),
	pages = [
		"About" => "about.md",
		"Usage" => [
			"usage/projectors.md",
			"usage/shadows.md"
		],
		"API" => [
			"types" => [
				"types/DensityMatrix.md"
			],
			"functions" => [
				"functions/FockOperators.md",
				"functions/MatrixOperations.md",
			]
		],
	],
)


deploydocs(
	repo = "github.com/thmhugo/LOShadows.jl.git",
)
