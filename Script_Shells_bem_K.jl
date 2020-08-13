# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	
#	Script Two concentric spheres. Backscattering versus frequency
#	
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# v. 13/8/20


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Meshes corresponding to spheres with radius R1 and R2 (R1 > R2)
	meshSphereR1 = "meshes/Mesh_Sphere_Exterior_1142.stl" ; 
	meshSphereR2 = "meshes/Mesh_Sphere_Interior_2274.stl" ; 

	println("External mesh : ")
	normales1, selv1, vertex1 = GetMesh_fast( meshSphereR1 ) ;
	println("Internal mesh : ")
	normales2, selv2, vertex2 = GetMesh_fast( meshSphereR2 ) ;
    
	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 	# Frecuencias a calcular
	f_0 = 10 ; # Frecuencias en Hz
	f_1 = 38000 ;
	N = 10 ; # Cantidad de frecuencias
	FrecArray = collect( Linspace( log10(f_0), log10(f_1), N ) ) ; # Array de números de onda
	FrecArray = 10 .^(FrecArray) ;
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Parámetros de los tres medios
	c0 = 1477.4 ;
	c1 = 1.04 * c0 ; 
	c2 = 0.23 * c0 ;
	rho0 = 1026.8 ;
	rho1 = 1.04 * rho0  ; 
	rho2 = 0.00129 * rho0 ;

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	kinc = [ -1.0, 0, 0 ] ; # Incidence direction (normalized to unity)

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Distributing all the structures to all the cores
	for i in cores
		SendToProc( i, selv1 = selv1, vertex1 = vertex1, normales1 = normales1, 
			selv2 = selv2, vertex2 = vertex2, normales2 = normales2, kinc = kinc ) ;
		SendToProc( i, rho0 = rho0, rho1 = rho1, rho2 = rho2 ) ;
	end

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Calculating the far field in parallel  
	finf = bem_shell_farfield_bs_K( FrecArray, kinc, c0, c1, c2, rho0, rho1, rho2,
			selv1, vertex1, normales1, selv2, vertex2, normales2, cores ) ;

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Saving to disk the TS
	writedlm( "out/shells_bem_spheres_k.dat", [ FrecArray  TS.( finf ) ] ) ;

	@everywhere GC.gc() ;
	
