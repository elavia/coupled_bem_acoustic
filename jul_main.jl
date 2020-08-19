# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	
#	Penetrable BEM environment definition
#	
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# v. 13/8/20

	# User parameters
	n_cores = 2 ; # Number of extra-cores 
	

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Julia packages
	
	using Distributed
	using DelimitedFiles
	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Local packages (only in the master node)
	include("src/jul_parallel.jl" ) ;  
	include("src/jul_meshes.jl") ; 
	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Build the processor's grid
	if nprocs() != ( n_cores + 1 )
		cores = addprocs( n_cores ) ;
	end
	println("Working with : ", nprocs() ," cores.") ;
	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Distributed packages (load in all the cores)
	# Routines for all the cores
	@everywhere using LinearAlgebra
	@everywhere using SpecialFunctions
	@everywhere include("src/jul_auxiliar.jl") ;  
	Distributed.@everywhere setprecision( 128 ) ; # Bits por arbitrary precision floats
	@everywhere include("src/bem_shells.jl") ;
	@everywhere include("src/bem_fluid.jl") ;
	@everywhere include("src/bem_funciones.jl") ;
	@everywhere include("src/bem_quadrules.jl") ;
	@everywhere include("src/bem_operadores.jl") ;
	@everywhere include("src/bem_parallel.jl") ;
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
