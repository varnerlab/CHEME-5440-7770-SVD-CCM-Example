### A Pluto.jl notebook ###
# v0.12.4

using Markdown
using InteractiveUtils

# ╔═╡ 4d3c839c-8656-11eb-2c53-5b065c79a757
begin
	
	# setup my paths -
	_PATH_TO_ROOT = pwd()
	_PATH_TO_MODEL = joinpath(_PATH_TO_ROOT,"model")
	
	# Activate project (and download any reqd packages)
	import Pkg
	Pkg.activate(_PATH_TO_ROOT)
	Pkg.instantiate()
	
	# load packages -
	using MAT
	using Plots
	using LinearAlgebra
	
	# return -
	nothing
end

# ╔═╡ 5f523bec-8658-11eb-01d8-4daa07ab805d
html"""<style>
main {
    max-width: 1200px;
  	width: 80%;
	margin: auto;
}
"""

# ╔═╡ 871b4af6-8658-11eb-07e4-492066309dde
md"""

## Strutural analysis of __E.coli__ central carbon metabolism
Let's take a look at the reaction network that is used by __E.coli__ to convert sugar (glucose) into energy (ATP) and the carbon building blocks required for growth. The _E.coli_ model that we'll explore is taken from the publication:

[Orth JD, Fleming RM, Palsson BØ. Reconstruction and Use of Microbial Metabolic Networks: the Core Escherichia coli Metabolic Model as an Educational Guide. EcoSal Plus. 2010 Sep;4(1). doi: 10.1128/ecosalplus.10.2.1. PMID: 26443778.](https://pubmed.ncbi.nlm.nih.gov/26443778/)

"""

# ╔═╡ 256c3e14-8662-11eb-3a85-fd511cba9e3e
md"""

##### Setup my Julia environment

"""

# ╔═╡ 1443255a-8658-11eb-32f5-6b721d5863f3
# some functions that we will use later -
function make_binary_stoichiometric_matrix(stoichiometric_matrix::Array{Float64,2})::Array{Float64,2}

    # initialize -
    (number_of_rows, number_of_cols) = size(stoichiometric_matrix)
    binary_stoichiometric_matrix = zeros(number_of_rows, number_of_cols)
    for row_index = 1:number_of_rows
        for col_index = 1:number_of_cols

            # grab value -
            old_value = stoichiometric_matrix[row_index, col_index]
            if (old_value !=0.0)
                binary_stoichiometric_matrix[row_index, col_index] = 1.0
            end
        end
    end

    # return -
    return binary_stoichiometric_matrix
end

# ╔═╡ 62c08610-865f-11eb-2d62-0f0dc6c54eac
function fractional_reconstruction_of_stoichiometric_matrix(stoichiometric_matrix,U,S,V)

    # what is the rank of S?
    rank_value = rank(diagm(S))
    (number_of_rows, number_of_cols) = size(stoichiometric_matrix)
    block_array = Array{Array{Float64,2},1}()

    fraction_block = zeros(number_of_rows, number_of_cols)
    for rank_index = 1:rank_value

        fraction_block = fraction_block + S[rank_index]*(U[:,rank_index]*transpose(V[:,rank_index]))
        push!(block_array,fraction_block)
    end

    return block_array
end

# ╔═╡ 149b1c7e-8660-11eb-09d9-d1ad428024b8
md"""

#### Analysis of the Stoichiometric Array S
The stoichiometric array $S$ is the digital representation of your best understanding of __all possible__ biological reactions that __can__ occur in the system of interest. The stoichiometric matrix $S$ has row dimension equal to the number of chemical species (often denoted as $\mathcal{M}$) while the column dimension is the number of possible reactions (often denoted as $\mathcal{R}$). 

There are some interesting structural features of $S$ (such as the connectivity) that help us understand the biology, and the role of certain metabolites and reactions. To see these features, we compute the species adjacency array $A_{\mathcal{M}}$ and the reaction adjacency array $A_{\mathcal{R}}$ from the __binary__ stoichiometric matrix.
"""

# ╔═╡ f06779e6-8656-11eb-3473-d5454656a725
begin
	
	# where is the mat file?
	path_to_cobra_mat_file = joinpath(_PATH_TO_MODEL,"modelReg.mat")
    model_name = "modelReg"
	
	# load the mat file -
	file = matopen(path_to_cobra_mat_file)
    cobra_model_dictionary = read(file,model_name)
    close(file) 
end

# ╔═╡ bd4725f0-8680-11eb-2fec-a7d625bec2e6
cobra_model_dictionary

# ╔═╡ 4d37f3ee-8657-11eb-3d66-a5b3c2ddb69b
# get the stoichiometric array from the mat model file -
stoichiometric_matrix = Matrix(cobra_model_dictionary["S"])

# ╔═╡ 1d2e460c-8667-11eb-2edd-ef0755503f0c
# How many metabolites and reactions in this array?
(number_of_metabolites, number_of_reactions) = size(stoichiometric_matrix)

# ╔═╡ 8bb8295a-8665-11eb-1acb-536bf3c24c4a
md"""
	
#### Binary stoichiometric matrix
We compute the __binary__ stoichiometric array (denoted as $\bar{S}$) from the original stoichiometric array by replacing all non-zero values with 1's. Thus, we remove directionality from the array and instead focus only the connectivity of the reaction network. 

"""

# ╔═╡ 11320712-8658-11eb-06ee-dd65c9463388
# make the stm binary -
binary_stoichiometric_matrix = make_binary_stoichiometric_matrix(stoichiometric_matrix)

# ╔═╡ 115d93fc-8665-11eb-2f0b-e3d56e7e0c88
md"""

##### Species adjacency array

The species adjacency array $A_{\mathcal{M}}$ is an $\mathcal{M}\times\mathcal{M}$ that gives connection information about the chemical species in the reaction network:

$A_{\mathcal{M}} = \bar{S}\bar{S}^{T}$

where $\bar{S}$ denotes the __binary__ stoichiometric array (all non-zero values replaced by 1) and $\bar{S}^{T}$ denotes the transpose of the binary stoichiometric array. The diagonal elements of $A_{\mathcal{M}}$ give the number of reactions chemical species $i$ participates in, while the $(i,j)$ off-diagonal gives the number of shared reactions between species $i$ and species $j$.
"""

# ╔═╡ 075c9c5a-865a-11eb-3d3d-8f77297e821c
# Metabolite adj matrix -
A_m = binary_stoichiometric_matrix*transpose(binary_stoichiometric_matrix)

# ╔═╡ b2644cc0-8681-11eb-1bf0-41fe17eaed38
[diag(A_m) cobra_model_dictionary["mets"]]

# ╔═╡ dc08ae90-8681-11eb-2814-b50a1550e3d7
begin
	most_connected_m_index = sortperm(diag(A_m),rev=true)
	cobra_model_dictionary["mets"][most_connected_m_index]
end

# ╔═╡ 2caf8728-8665-11eb-2c3f-b9a83c0cd228
md"""

##### Reaction adjacency array

The reaction adjacency array $A_{\mathcal{R}}$ is an $\mathcal{R}\times\mathcal{R}$ matrix that gives connection information about the chemical reactions in the reaction network:

$A_{\mathcal{R}} = \bar{S}^{T}\bar{S}$

where $\bar{S}$ denotes the __binary__ stoichiometric array (all non-zero values replaced by 1) and $\bar{S}^{T}$ denotes the transpose of the binary stoichiometric array. The diagonal elements of $A_{\mathcal{R}}$ give the number of chemical species that particpate in reaction $i$ , while the $(i,j)$ off-diagonal gives the number of shared metabolites between reaction $i$ and reaction $j$.
"""

# ╔═╡ 2a479544-865a-11eb-3924-5f52d70461e5
# Reaction adj matrix -
A_v = transpose(binary_stoichiometric_matrix)*binary_stoichiometric_matrix

# ╔═╡ e9e9726e-866e-11eb-35d4-ebba571d4830
diag(A_v)

# ╔═╡ 08b53c6e-866f-11eb-3ac0-55d6e6c588a9
cobra_model_dictionary["rxns"]

# ╔═╡ b78fb1b6-8682-11eb-27fd-01460c2b4d59
begin
	most_connected_v_index = sortperm(diag(A_v), rev=true)
	cobra_model_dictionary["rxns"][most_connected_v_index]
end

# ╔═╡ dda1acce-865a-11eb-12f1-69c460facd2e
md"""
	
#### Singular Value Decomposition (SVD)

The Singular Value Decomposition (SVD) of a matrix is a method to decompose an array $S$ into the product of three arrays:

	
$A = U\Sigma{V}^{T}$


The decomposition is the product of a rotation $\times$ stretch $\times$ rotation operations. The SVD of the stochiometric array is used to get additional structural insight into the reaction network, see the [Palsson study](https://pubmed.ncbi.nlm.nih.gov/12900206/) or the [Palsson structural analysis lecture video](https://www.youtube.com/watch?v=8W7RBt_iCYM).

The SVD factors a matrix into a set of weight modes, where each mode represents a section of the original matrix, where the weight (importance) of that section is given by the singular values. The original matrix can be reconstructed as the weighted sum of [outer products](https://en.wikipedia.org/wiki/Outer_product):
	
$A = \sum_{i=1}^{r}\sigma_{i}\left(u_{i}\otimes{v_{i}^{T}}\right)$

where $r$ denotes the [rank](https://en.wikipedia.org/wiki/System_of_linear_equations) of the matrix $A$

"""

# ╔═╡ 47611b1c-8666-11eb-3038-9b96fe5aca80
md"""

##### Rank?
In this context, we can think of rank (denoted as $r$) as the number of blocks that we need to reconstruct the stoichiometric array. Rank is given by:

$r\leq\min\left(\mathcal{M},\mathcal{R}\right)$

"""

# ╔═╡ 2f3ef28c-8666-11eb-0880-c5caf901d230
# what is the rank?
stm_rank = rank(stoichiometric_matrix)

# ╔═╡ 7e430732-8658-11eb-2084-6d5bcd88c606
# compute the SVD of the stoichiometric array -
(U,S,V) = svd(stoichiometric_matrix);

# ╔═╡ 7d3706ec-8658-11eb-2d21-fb700858314f
# check the accuracy of the decomposition -
error_term = norm(stoichiometric_matrix - U*diagm(S)*transpose(V))

# ╔═╡ e63045c8-865d-11eb-1053-597660a1025c
# the U matrix governs metabolite combinations -
U

# ╔═╡ f0de89bc-865d-11eb-049a-37c0caf8888f
begin
	
	# what is the col index that we want to look at?
	col_index_to_look_at = 1
	
	# sort the first U col (most imporant mode in species space) -
	sort_index_col_U = sortperm(abs.(U[:,col_index_to_look_at]), rev=true)
end

# ╔═╡ 53d61990-865e-11eb-27d4-5128d87bc2e5
U[sort_index_col_U,col_index_to_look_at]

# ╔═╡ f9ae8f70-8669-11eb-03e9-5bbfbdeb078f
# What is the species ranking in the fist mode?
cobra_model_dictionary["mets"][sort_index_col_U]

# ╔═╡ 873049ea-866e-11eb-0da6-9701b44d81a4
begin
	
	# what is the col index that we want to look at?
	v_col_index_to_look_at = 1
	
	# sort the first U col (most imporant mode in species space) -
	sort_index_col_V = sortperm(abs.(V[:,v_col_index_to_look_at]), rev=true)
end

# ╔═╡ a136439e-866e-11eb-1288-296081018807
V[sort_index_col_V,v_col_index_to_look_at]

# ╔═╡ b9cc81e6-866e-11eb-3a38-ad4797f2f017
cobra_model_dictionary["rxns"][sort_index_col_V]

# ╔═╡ aae2c418-865e-11eb-1ab3-a9beeb155cc1
# The S matrix is a diagonal matrix that contains the singular values (weights in each mode)
diagm(S)

# ╔═╡ ea21f162-865e-11eb-0921-e53cab1827d1
fraction = (S./sum(S))*100.0

# ╔═╡ ccf519ae-865e-11eb-04e8-07e16f8b1949
begin
	plot(fraction, legend=false)
	xlabel!("Mode index")
	ylabel!("Percentage reconstructed")
end

# ╔═╡ c1fe9faa-865e-11eb-1414-11daed944e6e
# Compute the "blocks" to put the STM back together -
FM = fractional_reconstruction_of_stoichiometric_matrix(stoichiometric_matrix, U,S,V);

# ╔═╡ cf3fea1a-865f-11eb-2dc5-9d8b2411713f
FM[1]

# ╔═╡ e80abe46-865f-11eb-2621-7ffa528500b8
FM[end]

# ╔═╡ 0956198a-8666-11eb-1d18-91f15871c324
begin
	
	# what is the error if we sum up the blocks?
	# the max number is stm_rank
	number_of_blocks = 67
	
	# compute the reconstruction error -
	reconstruction_error_term = norm(stoichiometric_matrix - FM[number_of_blocks])
end

# ╔═╡ cf09c3d6-865f-11eb-3204-6147885196fa


# ╔═╡ cebe0a54-865f-11eb-32c1-ab3a465b2d38


# ╔═╡ ce57578c-865f-11eb-30c2-4d149e65c536


# ╔═╡ ba525e28-865d-11eb-0787-5743c4d2f3cc


# ╔═╡ ba023e32-865d-11eb-0774-11450f2ffed0


# ╔═╡ Cell order:
# ╟─5f523bec-8658-11eb-01d8-4daa07ab805d
# ╟─871b4af6-8658-11eb-07e4-492066309dde
# ╟─256c3e14-8662-11eb-3a85-fd511cba9e3e
# ╠═4d3c839c-8656-11eb-2c53-5b065c79a757
# ╠═1443255a-8658-11eb-32f5-6b721d5863f3
# ╠═62c08610-865f-11eb-2d62-0f0dc6c54eac
# ╟─149b1c7e-8660-11eb-09d9-d1ad428024b8
# ╠═f06779e6-8656-11eb-3473-d5454656a725
# ╠═bd4725f0-8680-11eb-2fec-a7d625bec2e6
# ╠═4d37f3ee-8657-11eb-3d66-a5b3c2ddb69b
# ╠═1d2e460c-8667-11eb-2edd-ef0755503f0c
# ╟─8bb8295a-8665-11eb-1acb-536bf3c24c4a
# ╠═11320712-8658-11eb-06ee-dd65c9463388
# ╟─115d93fc-8665-11eb-2f0b-e3d56e7e0c88
# ╠═075c9c5a-865a-11eb-3d3d-8f77297e821c
# ╠═b2644cc0-8681-11eb-1bf0-41fe17eaed38
# ╠═dc08ae90-8681-11eb-2814-b50a1550e3d7
# ╟─2caf8728-8665-11eb-2c3f-b9a83c0cd228
# ╠═2a479544-865a-11eb-3924-5f52d70461e5
# ╠═e9e9726e-866e-11eb-35d4-ebba571d4830
# ╠═08b53c6e-866f-11eb-3ac0-55d6e6c588a9
# ╠═b78fb1b6-8682-11eb-27fd-01460c2b4d59
# ╟─dda1acce-865a-11eb-12f1-69c460facd2e
# ╟─47611b1c-8666-11eb-3038-9b96fe5aca80
# ╠═2f3ef28c-8666-11eb-0880-c5caf901d230
# ╠═7e430732-8658-11eb-2084-6d5bcd88c606
# ╠═7d3706ec-8658-11eb-2d21-fb700858314f
# ╠═e63045c8-865d-11eb-1053-597660a1025c
# ╠═f0de89bc-865d-11eb-049a-37c0caf8888f
# ╠═53d61990-865e-11eb-27d4-5128d87bc2e5
# ╠═f9ae8f70-8669-11eb-03e9-5bbfbdeb078f
# ╠═873049ea-866e-11eb-0da6-9701b44d81a4
# ╠═a136439e-866e-11eb-1288-296081018807
# ╠═b9cc81e6-866e-11eb-3a38-ad4797f2f017
# ╠═aae2c418-865e-11eb-1ab3-a9beeb155cc1
# ╠═ea21f162-865e-11eb-0921-e53cab1827d1
# ╠═ccf519ae-865e-11eb-04e8-07e16f8b1949
# ╠═c1fe9faa-865e-11eb-1414-11daed944e6e
# ╠═cf3fea1a-865f-11eb-2dc5-9d8b2411713f
# ╠═e80abe46-865f-11eb-2621-7ffa528500b8
# ╠═0956198a-8666-11eb-1d18-91f15871c324
# ╟─cf09c3d6-865f-11eb-3204-6147885196fa
# ╟─cebe0a54-865f-11eb-32c1-ab3a465b2d38
# ╟─ce57578c-865f-11eb-30c2-4d149e65c536
# ╟─ba525e28-865d-11eb-0787-5743c4d2f3cc
# ╟─ba023e32-865d-11eb-0774-11450f2ffed0
