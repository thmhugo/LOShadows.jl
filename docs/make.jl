using Documenter, LOShadows, SparseArrays

push!(LOAD_PATH, "../src/")

makedocs(
	sitename = "LOShadows.jl",
	remotes = nothing,
	modules = [LOShadows],
	format = Documenter.HTML(prettyurls = true),
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
)
