# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	
#	Julia shells scattering high-level routines
#	
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# v. 13/8/20

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Index:

# 	bem_shell_farfield_bs( pext::Array, K0::Real, K1::Real, K2::Real, rho0::Real, 
# 		rho1::Real, rho2::Real, selv1::Array, vertex1::Array, normales1::Array, 
# 		selv2::Array, vertex2::Array, normales2::Array, cores::Array )
# 	bem_shell_farfield_forward( pext::Array, K0::Real, K1::Real, K2::Real, rho0::Real, 
# 		rho1::Real, rho2::Real, selv1::Array, vertex1::Array, normales1::Array, selv2::Array, 
# 		vertex2::Array, normales2::Array, cores::Array )
# 	bem_shell_farfield_angular( pext::Array, kinc::Array, K0::Real, K1::Real, K2::Real, 
# 		rho0::Real, rho1::Real, rho2::Real, selv1::Array, vertex1::Array, normales1::Array, 
# 		selv2::Array, vertex2::Array, normales2::Array, cores::Array )
# 	bem_shell_nearfield( K0::Real, K1::Real, K2::Real, rho0::Real, rho1::Real, rho2::Real,
# 		pext::Array, pint1::Array, pint2::Array, kinc::Array, selv1::Array, vertex1::Array, 
# 		normales1::Array, selv2::Array, vertex2::Array, normales2::Array, cores::Array, TypeNumber )
# 	bem_shell_farfield_bs_K( FrecArray::Array, kinc::Array, c0::Real, c1::Real, c2::Real, 
# 		rho0::Real, rho1::Real, rho2::Real, selv1::Array, vertex1::Array, normales1::Array, 
# 		selv2::Array, vertex2::Array, normales2::Array, cores::Array )
# 	Fill_Matriz_Shell_diagonal( K0::Real, K1::Real, rho0::Real, rho1::Real, rango, 
# 		selv::Array, vertex::Array, normales::Array, TypeNumber )
# 	Fill_Matriz_Shell_offdiagonal( K::Real, rho::Real, rango, 
# 		selvColoc::Array, vertexColoc::Array, normalesColoc::Array,
# 		selvInteg::Array, vertexInteg::Array, normalesInteg::Array, TypeNumber )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	# Función que calcula el f_\infty de backscattering para un shell descripto por una mesh exterior 
	# 'selv1', 'vertex1', 'normales1' y una mesh interior 'selv2', 'vertex2', 'normales2'. Se calcula
	# para direcciones dadas por 'pext' (incidencias -pext).
	# Los cálculos se hacen en paralelo en el cluster (o en un único servidor) con los procesadores
	# dados por el arreglo 'cores' (en cores no está el proceso local "1"). 
	function bem_shell_farfield_bs( pext::Array, K0::Real, K1::Real, K2::Real, rho0::Real, rho1::Real, rho2::Real,
		selv1::Array, vertex1::Array, normales1::Array, selv2::Array, vertex2::Array, normales2::Array, cores::Array )
		# Parámetros
		NPE = size( pext )[1] ; # 
		N = size( selv1 )[1] ; # Semitamaño de la mesh exterior
		L = size( selv2 )[1] ; # Semitamaño de la mesh interior
		# Configuración de la ejecución en paralelo en el cluster
		ColPerProc = 24 ; # Cada procesador calculará 'ColPerProc' en cada tarea
		TypeNumber = ComplexF64 ; # Tipo de los números usados en los entes
		np = size( cores )[1] ; # Número de procesadores (sin el proceso local "1")
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Llenado de la matriz WKSPC con los operadores
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		println( "Filling WKSPC matrix with distributed process ... " ) ;
		# Creación de la matriz total WKSPC 
		WKSPC = Array{ TypeNumber }( undef, 2*N + 2*L, 2*N + 2*L ) ;
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Bloque "A" de 2*N x 2*N 
		Rangos = BuildRangos( N, ColPerProc ) ; # 
		RangosProc = BuildRangosProc( Rangos, np ) ; # Rangos de procesadores [1:np] salvo quizás el último
		@time for i = 1 : size( RangosProc )[1] # Un ciclo 'i' abarca 'np' cores simultáneamente
			Fut = Array{Future}( undef, RangosProc[i][end] ) ; # Arreglo de 'futuros'
			for q in RangosProc[i] # Recorre los procesadores del rango. Cada core hace una lonja de la matriz de 'Rangos[q+np(i-1)]' columnas
				Fut[ q ] = @spawnat( cores[q], Fill_Matriz_Shell_diagonal( eval(:K0), eval(:K1), eval(:rho0), eval(:rho1),
					   Rangos[ q + np*(i-1) ], eval(:selv1), eval(:vertex1), eval(:normales1), TypeNumber ) ) ;
			end
			for q in RangosProc[i]
				RangoCol = Rangos[ q + np * ( i - 1 ) ] ;
				WKSPC[ 1 : 2*N, RangoCol ] = fetch( Fut[ q ] )[ :, 1:size( RangoCol )[1] ] ;
				WKSPC[ 1 : 2*N, N .+ RangoCol ] = fetch( Fut[ q ] )[ :, size( RangoCol )[1]+1:end] ;
			end
			println("Done the chunk[A] : ", i ," / ", size( RangosProc )[1] ) ;
		end
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Bloque "B" de 2*N x 2*L 
		RangosB = BuildRangos( L, ColPerProc ) ; # 
		RangosProcB = BuildRangosProc( RangosB, np ) ; # Rangos de procesadores [1:np] salvo quizás el último
		@time for i = 1 : size( RangosProcB )[1] # Un ciclo 'i' abarca 'np' cores simultáneamente
			Fut = Array{Future}( undef, RangosProcB[i][end] ) ; # Arreglo de 'futuros'
			for q in RangosProcB[i] # Recorre los procesadores del rango. Cada core hace una lonja de la matriz de 'Rangos[q+np(i-1)]' columnas
				Fut[ q ] = @spawnat( cores[q], Fill_Matriz_Shell_offdiagonal( eval(:K1), eval(:rho1), RangosB[ q + np*(i-1) ], 
				eval(:selv1), eval(:vertex1), eval(:normales1), eval(:selv2), eval(:vertex2), eval(:normales2), TypeNumber ) ) ;
			end
			for q in RangosProcB[i]
				RangoColB = RangosB[ q + np * ( i - 1 ) ] ;
				WKSPC[ 1 : 2*N, 2*N .+ RangoColB ] = fetch( Fut[ q ] )[ :, 1:size( RangoColB )[1] ] ;
				WKSPC[ 1 : 2*N, 2*N + L .+ RangoColB ] = fetch( Fut[ q ] )[ :, size( RangoColB )[1]+1:end] ;
			end
			println( "Done the chunk[B] : ", i ," / ", size( RangosProcB )[1] ) ;
		end
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Bloque "C" de 2*L x 2*N 
		RangosC = BuildRangos( N, ColPerProc ) ; # 
		RangosProcC = BuildRangosProc( RangosC, np ) ; # Rangos de procesadores [1:np] salvo quizás el último
		@time for i = 1 : size( RangosProcC )[1] # Un ciclo 'i' abarca 'np' cores simultáneamente
			Fut = Array{Future}( undef, RangosProcC[i][end] ) ; # Arreglo de 'futuros'
			for q in RangosProcC[i] # Recorre los procesadores del rango. Cada core hace una lonja de la matriz de 'Rangos[q+np(i-1)]' columnas
				Fut[ q ] = @spawnat( cores[q], Fill_Matriz_Shell_offdiagonal( eval(:K1), eval(:rho1), RangosC[ q + np*(i-1) ], 
				eval(:selv2), eval(:vertex2), eval(:normales2), eval(:selv1), eval(:vertex1), eval(:normales1), TypeNumber ) ) ;
			end
			for q in RangosProcC[i]
				RangoColC = RangosC[ q + np * ( i - 1 ) ] ;
				WKSPC[ 2*N+1 : 2*(N+L), RangoColC ] = -fetch( Fut[ q ] )[ :, 1:size( RangoColC )[1] ] ;
				WKSPC[ 2*N+1 : 2*(N+L), N .+ RangoColC ] = -fetch( Fut[ q ] )[ :, size( RangoColC )[1]+1:end] ;
			end
			println( "Done the chunk[C] : ", i ," / ", size( RangosProcC )[1] ) ;
		end
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Bloque "D" de 2*L x 2*L 
		RangosD = BuildRangos( L, ColPerProc ) ; # 
		RangosProcD = BuildRangosProc( RangosD, np ) ; # Rangos de procesadores [1:np] salvo quizás el último
		@time for i = 1 : size( RangosProcD )[1] # Un ciclo 'i' abarca 'np' cores simultáneamente
			Fut = Array{Future}( undef, RangosProcD[i][end] ) ; # Arreglo de 'futuros'
			for q in RangosProcD[i] # Recorre los procesadores del rango. Cada core hace una lonja de la matriz de 'Rangos[q+np(i-1)]' columnas
				Fut[ q ] = @spawnat( cores[q], Fill_Matriz_Shell_diagonal( eval(:K1), eval(:K2), eval(:rho1), eval(:rho2),
				RangosD[ q + np*(i-1) ], eval(:selv2), eval(:vertex2), eval(:normales2), TypeNumber ) ) ;
			end
			for q in RangosProcD[i]
				RangoColD = RangosD[ q + np * ( i - 1 ) ] ;
				WKSPC[ 2*N+1 : 2*(N+L), 2*N .+ RangoColD ] = fetch( Fut[ q ] )[ :, 1:size( RangoColD )[1] ] ;
				WKSPC[ 2*N+1 : 2*(N+L), 2*N + L .+ RangoColD ] = fetch( Fut[ q ] )[ :, size( RangoColD )[1]+1:end] ;
			end
			println("Done the chunk[D] : ", i ," / ", size( RangosProcD )[1] ) ;
		end
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Vector con condiciones de borde
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		println( "Calculating the boundary values in parallel ... " ) ;
		BndValues = Array{ TypeNumber }( undef, 2 * ( N + L ), NPE ) ; # Array local de valores de borde
		# Se envía el cálculo a todos los cores de modo asíncrono y se espera a la finalización
		@time @sync @async for j = 1 : NPE 
			BndValues[ :, j ] = fetch( @spawn( FillBndValues_Shells( -pext[j,:], eval(:K0), eval(:rho0), 
						eval(:selv1), eval(:vertex1), eval(:normales1), L ) ) ) ;
		end
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Resolviendo el sistema 
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		println("Solving the system ... ") ;
		# Resolución del sistema para todos las incidencias
		# Se reescriben WKSPC (con autovectores) y BndValues (con la solución)
		@time LinearAlgebra.LAPACK.gesv!( WKSPC, BndValues ) ;
		# Liberación de recursos
		println("Freeing memory ... ") ;
		WKSPC = TypeNumber( 0 ) ;
		@everywhere GC.gc() ;
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Cálculo del farfield
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		println("Computing far field ... ") ;
		FarField_bs = Array{ TypeNumber }( undef, NPE ) ; # Array local de farfields
		# No vale la pena paralelizar este cálculo para valores 'NPE' no monstruosos por el overhead.
		# BndValues es ahora la "densidad" solución
		@time for j = 1 : NPE 
			FarField_bs[ j ] = FarField_Calculation( K0, pext[j,:], rho0, 
				selv1, vertex1, normales1, view( BndValues, :, j ) ) ;
		end
		return FarField_bs
	end
	

	function bem_shell_farfield_forward( pext::Array, K0::Real, K1::Real, K2::Real, rho0::Real, 
		rho1::Real, rho2::Real, selv1::Array, vertex1::Array, normales1::Array, selv2::Array, 
		vertex2::Array, normales2::Array, cores::Array )
		# Parámetros
		NPE = size( pext )[1] ; # 
		N = size( selv1 )[1] ; # Semitamaño de la mesh exterior
		L = size( selv2 )[1] ; # Semitamaño de la mesh interior
		# Configuración de la ejecución en paralelo en el cluster
		ColPerProc = 24 ; # Cada procesador calculará 'ColPerProc' en cada tarea
		TypeNumber = ComplexF64 ; # Tipo de los números usados en los entes
		np = size( cores )[1] ; # Número de procesadores (sin el proceso local "1")
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Llenado de la matriz WKSPC con los operadores
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		println( "Filling WKSPC matrix with distributed process ... " ) ;
		# Creación de la matriz total WKSPC 
		WKSPC = Array{ TypeNumber }( undef, 2*N + 2*L, 2*N + 2*L ) ;
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Bloque "A" de 2*N x 2*N 
		Rangos = BuildRangos( N, ColPerProc ) ; # 
		RangosProc = BuildRangosProc( Rangos, np ) ; # Rangos de procesadores [1:np] salvo quizás el último
		@time for i = 1 : size( RangosProc )[1] # Un ciclo 'i' abarca 'np' cores simultáneamente
			Fut = Array{Future}( undef, RangosProc[i][end] ) ; # Arreglo de 'futuros'
			for q in RangosProc[i] # Recorre los procesadores del rango. Cada core hace una lonja de la matriz de 'Rangos[q+np(i-1)]' columnas
				Fut[ q ] = @spawnat( cores[q], Fill_Matriz_Shell_diagonal( eval(:K0), 
								eval(:K1), eval(:rho0), eval(:rho1),
					   Rangos[ q + np*(i-1) ], eval(:selv1), eval(:vertex1), 
								eval(:normales1), TypeNumber ) ) ;
			end
			for q in RangosProc[i]
				RangoCol = Rangos[ q + np * ( i - 1 ) ] ;
				WKSPC[ 1 : 2*N, RangoCol ] = fetch( Fut[ q ] )[ :, 1:size( RangoCol )[1] ] ;
				WKSPC[ 1 : 2*N, N .+ RangoCol ] = fetch( Fut[ q ] )[ :, size( RangoCol )[1]+1:end] ;
			end
			println("Done the chunk[A] : ", i ," / ", size( RangosProc )[1] ) ;
		end
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Bloque "B" de 2*N x 2*L 
		RangosB = BuildRangos( L, ColPerProc ) ; # 
		RangosProcB = BuildRangosProc( RangosB, np ) ; # Rangos de procesadores [1:np] salvo quizás el último
		@time for i = 1 : size( RangosProcB )[1] # Un ciclo 'i' abarca 'np' cores simultáneamente
			Fut = Array{Future}( undef, RangosProcB[i][end] ) ; # Arreglo de 'futuros'
			for q in RangosProcB[i] # Recorre los procesadores del rango. Cada core hace una lonja de la matriz de 'Rangos[q+np(i-1)]' columnas
				Fut[ q ] = @spawnat( cores[q], Fill_Matriz_Shell_offdiagonal( eval(:K1), eval(:rho1), RangosB[ q + np*(i-1) ], 
				eval(:selv1), eval(:vertex1), eval(:normales1), eval(:selv2), eval(:vertex2), eval(:normales2), TypeNumber ) ) ;
			end
			for q in RangosProcB[i]
				RangoColB = RangosB[ q + np * ( i - 1 ) ] ;
				WKSPC[ 1 : 2*N, 2*N .+ RangoColB ] = fetch( Fut[ q ] )[ :, 1:size( RangoColB )[1] ] ;
				WKSPC[ 1 : 2*N, 2*N + L .+ RangoColB ] = fetch( Fut[ q ] )[ :, size( RangoColB )[1]+1:end] ;
			end
			println( "Done the chunk[B] : ", i ," / ", size( RangosProcB )[1] ) ;
		end
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Bloque "C" de 2*L x 2*N 
		RangosC = BuildRangos( N, ColPerProc ) ; # 
		RangosProcC = BuildRangosProc( RangosC, np ) ; # Rangos de procesadores [1:np] salvo quizás el último
		@time for i = 1 : size( RangosProcC )[1] # Un ciclo 'i' abarca 'np' cores simultáneamente
			Fut = Array{Future}( undef, RangosProcC[i][end] ) ; # Arreglo de 'futuros'
			for q in RangosProcC[i] # Recorre los procesadores del rango. Cada core hace una lonja de la matriz de 'Rangos[q+np(i-1)]' columnas
				Fut[ q ] = @spawnat( cores[q], Fill_Matriz_Shell_offdiagonal( eval(:K1), eval(:rho1), RangosC[ q + np*(i-1) ], 
				eval(:selv2), eval(:vertex2), eval(:normales2), eval(:selv1), eval(:vertex1), eval(:normales1), TypeNumber ) ) ;
			end
			for q in RangosProcC[i]
				RangoColC = RangosC[ q + np * ( i - 1 ) ] ;
				WKSPC[ 2*N+1 : 2*(N+L), RangoColC ] = -fetch( Fut[ q ] )[ :, 1:size( RangoColC )[1] ] ;
				WKSPC[ 2*N+1 : 2*(N+L), N .+ RangoColC ] = -fetch( Fut[ q ] )[ :, size( RangoColC )[1]+1:end] ;
			end
			println( "Done the chunk[C] : ", i ," / ", size( RangosProcC )[1] ) ;
		end
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Bloque "D" de 2*L x 2*L 
		RangosD = BuildRangos( L, ColPerProc ) ; # 
		RangosProcD = BuildRangosProc( RangosD, np ) ; # Rangos de procesadores [1:np] salvo quizás el último
		@time for i = 1 : size( RangosProcD )[1] # Un ciclo 'i' abarca 'np' cores simultáneamente
			Fut = Array{Future}( undef, RangosProcD[i][end] ) ; # Arreglo de 'futuros'
			for q in RangosProcD[i] # Recorre los procesadores del rango. Cada core hace una lonja de la matriz de 'Rangos[q+np(i-1)]' columnas
				Fut[ q ] = @spawnat( cores[q], Fill_Matriz_Shell_diagonal( eval(:K1), eval(:K2), eval(:rho1), eval(:rho2),
				RangosD[ q + np*(i-1) ], eval(:selv2), eval(:vertex2), eval(:normales2), TypeNumber ) ) ;
			end
			for q in RangosProcD[i]
				RangoColD = RangosD[ q + np * ( i - 1 ) ] ;
				WKSPC[ 2*N+1 : 2*(N+L), 2*N .+ RangoColD ] = fetch( Fut[ q ] )[ :, 1:size( RangoColD )[1] ] ;
				WKSPC[ 2*N+1 : 2*(N+L), 2*N + L .+ RangoColD ] = fetch( Fut[ q ] )[ :, size( RangoColD )[1]+1:end] ;
			end
			println("Done the chunk[D] : ", i ," / ", size( RangosProcD )[1] ) ;
		end
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Vector con condiciones de borde
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		println( "Calculating the boundary values in parallel ... " ) ;
		BndValues = Array{ TypeNumber }( undef, 2 * ( N + L ), NPE ) ; # Array local de valores de borde
		# Se envía el cálculo a todos los cores de modo asíncrono y se espera a la finalización
		@time @sync @async for j = 1 : NPE 
			# La única diferencia está en que ahora la dirección no es la 'kinc == pext'
			BndValues[ :, j ] = fetch( @spawn( FillBndValues_Shells( pext[j,:], eval(:K0), eval(:rho0), 
						eval(:selv1), eval(:vertex1), eval(:normales1), L ) ) ) ;
		end
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Resolviendo el sistema 
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		println("Solving the system ... ") ;
		# Resolución del sistema para todos las incidencias
		# Se reescriben WKSPC (con autovectores) y BndValues (con la solución)
		@time LinearAlgebra.LAPACK.gesv!( WKSPC, BndValues ) ;
		# Liberación de recursos
		println("Freeing memory ... ") ;
		WKSPC = TypeNumber( 0 ) ;
		@everywhere GC.gc() ;
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Cálculo del farfield
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		println("Computing far field ... ") ;
		FarField_bs = Array{ TypeNumber }( undef, NPE ) ; # Array local de farfields
		# No vale la pena paralelizar este cálculo para valores 'NPE' no monstruosos por el overhead.
		# BndValues es ahora la "densidad" solución
		@time for j = 1 : NPE 
			FarField_bs[ j ] = FarField_Calculation( K0, pext[j,:], rho0, selv1, 
					vertex1, normales1, view( BndValues, :, j ) ) ;
		end
		return FarField_bs
	end
	
	# Función que calcula el f_\infty para un shell descripto por una mesh exterior 
	# 'selv1', 'vertex1', 'normales1' y una mesh interior 'selv2', 'vertex2', 'normales2'. Se calcula
	# para direcciones dadas por 'pext' para una dirección de incidencia 'kinc'.
	# Los cálculos se hacen en paralelo en el cluster (o en un único servidor) con los procesadores
	# dados por el arreglo 'cores' (en cores no está el proceso local "1"). 

	function bem_shell_farfield_angular( pext::Array, kinc::Array, K0::Real, K1::Real, K2::Real, 
		rho0::Real, rho1::Real, rho2::Real, selv1::Array, vertex1::Array, normales1::Array, 
		selv2::Array, vertex2::Array, normales2::Array, cores::Array )
		# Parámetros
		NPE = size( pext )[1] ; #
		N = size( selv1 )[1] ; # Semitamaño de la mesh exterior
		L = size( selv2 )[1] ; # Semitamaño de la mesh interior
		# Configuración de la ejecución en paralelo en el cluster
		ColPerProc = 24 ; # Cada procesador calculará 'ColPerProc' en cada tarea
		TypeNumber = ComplexF64 ; # Tipo de los números usados en los entes
		np = size( cores )[1] ; # Número de procesadores (sin el proceso local "1")
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Llenado de la matriz WKSPC con los operadores
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		println( "Filling WKSPC matrix with distributed process ... " ) ;
		# Creación de la matriz total WKSPC 
		WKSPC = Array{ TypeNumber }( undef, 2*N + 2*L, 2*N + 2*L ) ;
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Bloque "A" de 2*N x 2*N 
		Rangos = BuildRangos( N, ColPerProc ) ; # 
		RangosProc = BuildRangosProc( Rangos, np ) ; # Rangos de procesadores [1:np] salvo quizás el último
		@time for i = 1 : size( RangosProc )[1] # Un ciclo 'i' abarca 'np' cores simultáneamente
			Fut = Array{Future}( undef, RangosProc[i][end] ) ; # Arreglo de 'futuros'
			for q in RangosProc[i] # Recorre los procesadores del rango. Cada core hace una lonja de la matriz de 'Rangos[q+np(i-1)]' columnas
				Fut[ q ] = @spawnat( cores[q], Fill_Matriz_Shell_diagonal( eval(:K0), eval(:K1), eval(:rho0), eval(:rho1),
					   Rangos[ q + np*(i-1) ], eval(:selv1), eval(:vertex1), eval(:normales1), TypeNumber ) ) ;
			end
			for q in RangosProc[i]
				RangoCol = Rangos[ q + np * ( i - 1 ) ] ;
				WKSPC[ 1 : 2*N, RangoCol ] = fetch( Fut[ q ] )[ :, 1:size( RangoCol )[1] ] ;
				WKSPC[ 1 : 2*N, N .+ RangoCol ] = fetch( Fut[ q ] )[ :, size( RangoCol )[1]+1:end] ;
			end
			println("Done the chunk[A] : ", i ," / ", size( RangosProc )[1] ) ;
		end
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Bloque "B" de 2*N x 2*L 
		RangosB = BuildRangos( L, ColPerProc ) ; # 
		RangosProcB = BuildRangosProc( RangosB, np ) ; # Rangos de procesadores [1:np] salvo quizás el último
		@time for i = 1 : size( RangosProcB )[1] # Un ciclo 'i' abarca 'np' cores simultáneamente
			Fut = Array{Future}( undef, RangosProcB[i][end] ) ; # Arreglo de 'futuros'
			for q in RangosProcB[i] # Recorre los procesadores del rango. Cada core hace una lonja de la matriz de 'Rangos[q+np(i-1)]' columnas
				Fut[ q ] = @spawnat( cores[q], Fill_Matriz_Shell_offdiagonal( eval(:K1), eval(:rho1), RangosB[ q + np*(i-1) ], 
				eval(:selv1), eval(:vertex1), eval(:normales1), eval(:selv2), eval(:vertex2), eval(:normales2), TypeNumber ) ) ;
			end
			for q in RangosProcB[i]
				RangoColB = RangosB[ q + np * ( i - 1 ) ] ;
				WKSPC[ 1 : 2*N, 2*N .+ RangoColB ] = fetch( Fut[ q ] )[ :, 1:size( RangoColB )[1] ] ;
				WKSPC[ 1 : 2*N, 2*N + L .+ RangoColB ] = fetch( Fut[ q ] )[ :, size( RangoColB )[1]+1:end] ;
			end
			println("Done the chunk[B] : ", i ," / ", size( RangosProcB )[1] ) ;
		end
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Bloque "C" de 2*L x 2*N 
		RangosC = BuildRangos( N, ColPerProc ) ; # 
		RangosProcC = BuildRangosProc( RangosC, np ) ; # Rangos de procesadores [1:np] salvo quizás el último
		@time for i = 1 : size( RangosProcC )[1] # Un ciclo 'i' abarca 'np' cores simultáneamente
			Fut = Array{Future}( undef, RangosProcC[i][end] ) ; # Arreglo de 'futuros'
			for q in RangosProcC[i] # Recorre los procesadores del rango. Cada core hace una lonja de la matriz de 'Rangos[q+np(i-1)]' columnas
				Fut[ q ] = @spawnat( cores[q], Fill_Matriz_Shell_offdiagonal( eval(:K1), eval(:rho1), RangosC[ q + np*(i-1) ], 
			eval(:selv2), eval(:vertex2), eval(:normales2), eval(:selv1), eval(:vertex1), eval(:normales1), TypeNumber ) ) ;
			end
			for q in RangosProcC[i]
				RangoColC = RangosC[ q + np * ( i - 1 ) ] ;
				WKSPC[ 2*N+1 : 2*(N+L), RangoColC ] = -fetch( Fut[ q ] )[ :, 1:size( RangoColC )[1] ] ;
				WKSPC[ 2*N+1 : 2*(N+L), N .+ RangoColC ] = -fetch( Fut[ q ] )[ :, size( RangoColC )[1]+1:end] ;
			end
			println("Done the chunk[C] : ", i ," / ", size( RangosProcC )[1] ) ;
		end
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Bloque "D" de 2*L x 2*L 
		RangosD = BuildRangos( L, ColPerProc ) ; # 
		RangosProcD = BuildRangosProc( RangosD, np ) ; # Rangos de procesadores [1:np] salvo quizás el último
		@time for i = 1 : size( RangosProcD )[1] # Un ciclo 'i' abarca 'np' cores simultáneamente
			Fut = Array{Future}( undef, RangosProcD[i][end] ) ; # Arreglo de 'futuros'
			for q in RangosProcD[i] # Recorre los procesadores del rango. Cada core hace una lonja de la matriz de 'Rangos[q+np(i-1)]' columnas
				Fut[ q ] = @spawnat( cores[q], Fill_Matriz_Shell_diagonal( eval(:K1), eval(:K2), eval(:rho1), eval(:rho2),
					   RangosD[ q + np*(i-1) ], eval(:selv2), eval(:vertex2), eval(:normales2), TypeNumber ) ) ;
			end
			for q in RangosProcD[i]
				RangoColD = RangosD[ q + np * ( i - 1 ) ] ;
				WKSPC[ 2*N+1 : 2*(N+L), 2*N .+ RangoColD ] = fetch( Fut[ q ] )[ :, 1:size( RangoColD )[1] ] ;
				WKSPC[ 2*N+1 : 2*(N+L), 2*N + L .+ RangoColD ] = fetch( Fut[ q ] )[ :, size( RangoColD )[1]+1:end] ;
			end
			println("Done the chunk[D] : ", i ," / ", size( RangosProcD )[1] ) ;
		end
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Vector con condiciones de borde
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		println( "Calculating the boundary values ... " ) ;
		# No se justifica hacerlo en paralelo
		BndValues = FillBndValues_Shells( kinc, K0, rho0, selv1, vertex1, normales1, L ) ;
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Resolviendo el sistema 
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		println("Solving the system ... ") ;
		# Resolución del sistema para todos las incidencias
		# Se reescriben WKSPC (con autovectores) y BndValues (con la solución)
		@time LinearAlgebra.LAPACK.gesv!( WKSPC, BndValues ) ;
		# Liberación de recursos
		println("Freeing memory ... ") ;
		WKSPC = TypeNumber( 0 ) ;
		@everywhere GC.gc() ; 
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Cálculo del farfield
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		println("Computing far field ... ") ;
		# Esta paralelización trivial parece ser conveniente si 'BndValues' no es monstruosa.
		# BndValues es ahora la "densidad" solución
		for i in cores # Se copia a los cores BndValues
			SendToProc( i, BndValues = BndValues ) ;
		end
		pextV = [ pext[i,:] for i=1:NPE ] ; # Se convierte a un vector de incidencias
		# Calculo en paralelo entre los cores
		@time FarField_angular = pmap( x -> FarField_Calculation( eval(:K0), x, eval(:rho0), 
				eval(:selv1), eval(:vertex1), eval(:normales1), eval(:BndValues) ), pextV ) ;
		return FarField_angular
	end


	function bem_shell_nearfield( K0::Real, K1::Real, K2::Real, rho0::Real, rho1::Real, rho2::Real,
		pext::Array, pint1::Array, pint2::Array, kinc::Array, selv1::Array, vertex1::Array, 
		normales1::Array, selv2::Array, vertex2::Array, normales2::Array, cores::Array, TypeNumber )
		# Parámetros
		NPE = size( pext )[1] ; #
		NPI1 = size( pint1 )[1] ; # Número de puntos internos (medio 1) de campo 
		NPI2 = size( pint2 )[1] ; # Número de puntos internos (medio 2) de campo 
		A0 = 1.0 ; # Amplitud de la onda plana
		N = size( selv1 )[1] ; # Semitamaño de la mesh exterior
		L = size( selv2 )[1] ; # Semitamaño de la mesh interior
		# Configuración de la ejecución en paralelo en el cluster
		ColPerProc = 24 ; # Cada procesador calculará 'ColPerProc' en cada tarea
		TypeNumber = ComplexF64 ; # Tipo de los números usados en los entes
		np = size( cores )[1] ; # Número de procesadores (sin el proceso local "1")
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Llenado de la matriz WKSPC con los operadores
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		println( "Filling WKSPC matrix with distributed process ... " ) ;
		# Creación de la matriz total WKSPC 
		WKSPC = Array{ TypeNumber }( undef, 2*N + 2*L, 2*N + 2*L ) ;
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Bloque "A" de 2*N x 2*N 
		Rangos = BuildRangos( N, ColPerProc ) ; # 
		RangosProc = BuildRangosProc( Rangos, np ) ; # Rangos de procesadores [1:np] salvo quizás el último
		@time for i = 1 : size( RangosProc )[1] # Un ciclo 'i' abarca 'np' cores simultáneamente
			Fut = Array{Future}( undef, RangosProc[i][end] ) ; # Arreglo de 'futuros'
			for q in RangosProc[i] # Recorre los procesadores del rango. Cada core hace una lonja de la matriz de 'Rangos[q+np(i-1)]' columnas
				Fut[ q ] = @spawnat( cores[q], Fill_Matriz_Shell_diagonal( eval(:K0), eval(:K1), eval(:rho0), eval(:rho1),
					Rangos[ q + np*(i-1) ], eval(:selv1), eval(:vertex1), eval(:normales1), TypeNumber ) ) ;
			end
			for q in RangosProc[i]
				RangoCol = Rangos[ q + np * ( i - 1 ) ] ;
				WKSPC[ 1 : 2*N, RangoCol ] = fetch( Fut[ q ] )[ :, 1:size( RangoCol )[1] ] ;
				WKSPC[ 1 : 2*N, N .+ RangoCol ] = fetch( Fut[ q ] )[ :, size( RangoCol )[1]+1:end] ;
			end
			println("Done the chunk[A] : ", i ," / ", size( RangosProc )[1] ) ;
		end
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Bloque "B" de 2*N x 2*L 
		RangosB = BuildRangos( L, ColPerProc ) ; # 
		RangosProcB = BuildRangosProc( RangosB, np ) ; # Rangos de procesadores [1:np] salvo quizás el último
		@time for i = 1 : size( RangosProcB )[1] # Un ciclo 'i' abarca 'np' cores simultáneamente
			Fut = Array{Future}( undef, RangosProcB[i][end] ) ; # Arreglo de 'futuros'
			for q in RangosProcB[i] # Recorre los procesadores del rango. Cada core hace una lonja de la matriz de 'Rangos[q+np(i-1)]' columnas
				Fut[ q ] = @spawnat( cores[q], Fill_Matriz_Shell_offdiagonal( eval(:K1), eval(:rho1), RangosB[ q + np*(i-1) ], 
				eval(:selv1), eval(:vertex1), eval(:normales1), eval(:selv2), eval(:vertex2), eval(:normales2), TypeNumber ) ) ;
			end
			for q in RangosProcB[i]
				RangoColB = RangosB[ q + np * ( i - 1 ) ] ;
				WKSPC[ 1 : 2*N, 2*N .+ RangoColB ] = fetch( Fut[ q ] )[ :, 1:size( RangoColB )[1] ] ;
				WKSPC[ 1 : 2*N, 2*N + L .+ RangoColB ] = fetch( Fut[ q ] )[ :, size( RangoColB )[1]+1:end] ;
			end
			println("Done the chunk[B] : ", i ," / ", size( RangosProcB )[1] ) ;
		end
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Bloque "C" de 2*L x 2*N 
		RangosC = BuildRangos( N, ColPerProc ) ; # 
		RangosProcC = BuildRangosProc( RangosC, np ) ; # Rangos de procesadores [1:np] salvo quizás el último
		@time for i = 1 : size( RangosProcC )[1] # Un ciclo 'i' abarca 'np' cores simultáneamente
			Fut = Array{Future}( undef, RangosProcC[i][end] ) ; # Arreglo de 'futuros'
			for q in RangosProcC[i] # Recorre los procesadores del rango. Cada core hace una lonja de la matriz de 'Rangos[q+np(i-1)]' columnas
				Fut[ q ] = @spawnat( cores[q], Fill_Matriz_Shell_offdiagonal( eval(:K1), eval(:rho1), RangosC[ q + np*(i-1) ], 
				eval(:selv2), eval(:vertex2), eval(:normales2), eval(:selv1), eval(:vertex1), eval(:normales1), TypeNumber ) ) ;
			end
			for q in RangosProcC[i]
				RangoColC = RangosC[ q + np * ( i - 1 ) ] ;
				WKSPC[ 2*N+1 : 2*(N+L), RangoColC ] = -fetch( Fut[ q ] )[ :, 1:size( RangoColC )[1] ] ;
				WKSPC[ 2*N+1 : 2*(N+L), N .+ RangoColC ] = -fetch( Fut[ q ] )[ :, size( RangoColC )[1]+1:end] ;
			end
			println("Done the chunk[C] : ", i ," / ", size( RangosProcC )[1] ) ;
		end
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Bloque "D" de 2*L x 2*L 
		RangosD = BuildRangos( L, ColPerProc ) ; # 
		RangosProcD = BuildRangosProc( RangosD, np ) ; # Rangos de procesadores [1:np] salvo quizás el último
		@time for i = 1 : size( RangosProcD )[1] # Un ciclo 'i' abarca 'np' cores simultáneamente
			Fut = Array{Future}( undef, RangosProcD[i][end] ) ; # Arreglo de 'futuros'
			for q in RangosProcD[i] # Recorre los procesadores del rango. Cada core hace una lonja de la matriz de 'Rangos[q+np(i-1)]' columnas
				Fut[ q ] = @spawnat( cores[q], Fill_Matriz_Shell_diagonal( eval(:K1), eval(:K2), eval(:rho1), eval(:rho2),
					   RangosD[ q + np*(i-1) ], eval(:selv2), eval(:vertex2), eval(:normales2), TypeNumber ) ) ;
			end
			for q in RangosProcD[i]
				RangoColD = RangosD[ q + np * ( i - 1 ) ] ;
				WKSPC[ 2*N+1 : 2*(N+L), 2*N .+ RangoColD ] = fetch( Fut[ q ] )[ :, 1:size( RangoColD )[1] ] ;
				WKSPC[ 2*N+1 : 2*(N+L), 2*N + L .+ RangoColD ] = fetch( Fut[ q ] )[ :, size( RangoColD )[1]+1:end] ;
			end
			println("Done the chunk[D] : ", i ," / ", size( RangosProcD )[1] ) ;
		end
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Vector con condiciones de borde
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		println( "Calculating the boundary values ... " ) ;
		# No se justifica hacerlo en paralelo
		BndValues = FillBndValues_Shells( kinc, K0, rho0, selv1, vertex1, normales1, L ) ;
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Resolviendo el sistema 
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		println("Solving the system ... ") ;
		# Resolución del sistema para todos las incidencias
		# Se reescriben WKSPC (con autovectores) y BndValues (con la solución)
		@time LinearAlgebra.LAPACK.gesv!( WKSPC, BndValues ) ;
		# Liberación de recursos
		println("Freeing memory ... ") ;
		WKSPC = TypeNumber( 0 ) ;
		@everywhere GC.gc() ; 
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Cálculo del nearfield
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		println("Computing field in the external / internal points ... ") ;
		# Esta paralelización trivial parece ser conveniente si 'BndValues' no es monstruosa.
		# BndValues es ahora la "densidad" solución
		Density1 = BndValues[1:2*N] ;
		Density2 = BndValues[2*N+1:end] ;
		for i in cores # Se copia a los cores BndValues
			SendToProc( i, BndValues = BndValues, Density1 = Density1, Density2 = Density2 ) ;
		end
		pextV = [ pext[i,:] for i=1:NPE ] ; # Se convierte a vector
		# Calculo en paralelo entre los cores.
		# Se selecciona la parte de la densidad en la mesh '1'
		@time ScaField = pmap( x -> NearField_shell( eval(:K0), x, eval(:rho0), eval(:selv1), 
						eval(:vertex1), eval(:normales1), eval(:Density1) ), pextV ) ;
		@time IncField = pmap( x -> IncidentField_Calculation( A0, eval(:K0), eval(:kinc), x ),
				 		pextV ) ;
		pint1V = [ pint1[i,:] for i=1:NPI1 ] ; # Se convierte a vector
		pint2V = [ pint2[i,:] for i=1:NPI2 ] ; # Se convierte a vector
		# Calculo en paralelo entre los cores
		@time Int1Field = pmap( x -> NearField_shell_medio1( eval(:K1), x, eval(:rho1), 
				eval(:selv1), eval(:vertex1), eval(:normales1),
				eval(:selv2), eval(:vertex2), eval(:normales2), eval(:BndValues) ), pint1V ) ;
		@time Int2Field = pmap( x -> NearField_shell( eval(:K2), x, eval(:rho2), eval(:selv2), 
						eval(:vertex2), eval(:normales2), eval(:Density2) ), pint2V ) ;
		return IncField, ScaField, Int1Field, Int2Field ;
	end


	# Función que calcula el f_\infty de backscattering para un shell descripto por una mesh exterior 
	# 'selv1', 'vertex1', 'normales1' y una mesh interior 'selv2', 'vertex2', 'normales2'. Se calcula
	# para una dirección de incidencia fija 'kinc' y para un array de frecuencias 'FrecArray' que determina
	# tres arrays de números de onda.
	# Los cálculos se hacen en paralelo en el cluster (o en un único servidor) con los procesadores
	# dados por el arreglo 'cores' (en cores no está el proceso local "1"). 

	function bem_shell_farfield_bs_K( FrecArray::Array, kinc::Array, c0::Real, c1::Real, c2::Real, 
		rho0::Real, rho1::Real, rho2::Real, selv1::Array, vertex1::Array, normales1::Array, 
		selv2::Array, vertex2::Array, normales2::Array, cores::Array )
		# Parámetros
		NFREC = size( FrecArray )[1] ; #
		N = size( selv1 )[1] ; # Semitamaño de la mesh exterior
		L = size( selv2 )[1] ; # Semitamaño de la mesh interior
		# Configuración de la ejecución en paralelo en el cluster
		ColPerProc = 24 ; # Cada procesador calculará 'ColPerProc' en cada tarea
		TypeNumber = ComplexF64 ; # Tipo de los números usados en los entes
		FarField_bs_K = Array{ TypeNumber }( undef, NFREC ) ;
		np = size( cores )[1] ; # Número de procesadores (sin el proceso local "1")
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Ciclo en frecuencias 
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		for w = 1 : NFREC
			K0 = 2 * pi / c0 * FrecArray[ w] ; # Números de onda actuales
			K1 = 2 * pi / c1 * FrecArray[ w ] ; # Números de onda actuales
			K2 = 2 * pi / c2 * FrecArray[ w ] ; # Números de onda actuales
			println("Calculating the frequency : ", FrecArray[w], " ::: ", w," of ", NFREC ) ;
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# Llenado de la matriz WKSPC con los operadores
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# Creación de la matriz total WKSPC 
			WKSPC = Array{ TypeNumber }( undef, 2*N + 2*L, 2*N + 2*L ) ;
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# Bloque "A" de 2*N x 2*N 
			Rangos = BuildRangos( N, ColPerProc ) ; # 
			RangosProc = BuildRangosProc( Rangos, np ) ; # Rangos de procesadores [1:np] salvo quizás el último
			@time for i = 1 : size( RangosProc )[1] # Un ciclo 'i' abarca 'np' cores simultáneamente
				Fut = Array{Future}( undef, RangosProc[i][end] ) ; # Arreglo de 'futuros'
				for q in RangosProc[i] # Recorre los procesadores del rango. Cada core hace una lonja de la matriz de 'Rangos[q+np(i-1)]' columnas
					Fut[ q ] = @spawnat( cores[q], Fill_Matriz_Shell_diagonal( K0, K1, eval(:rho0), eval(:rho1),
						Rangos[ q + np*(i-1) ], eval(:selv1), eval(:vertex1), eval(:normales1), TypeNumber ) ) ;
				end
				for q in RangosProc[i]
					RangoCol = Rangos[ q + np * ( i - 1 ) ] ;
					WKSPC[ 1 : 2*N, RangoCol ] = fetch( Fut[ q ] )[ :, 1:size( RangoCol )[1] ] ;
					WKSPC[ 1 : 2*N, N .+ RangoCol ] = fetch( Fut[ q ] )[ :, size( RangoCol )[1]+1:end] ;
				end
				println("Done with the chunk [A] : ", i ," / ", size( RangosProc )[1] ) ;
			end
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# Bloque "B" de 2*N x 2*L 
			RangosB = BuildRangos( L, ColPerProc ) ; # 
			RangosProcB = BuildRangosProc( RangosB, np ) ; # Rangos de procesadores [1:np] salvo quizás el último
			@time for i = 1 : size( RangosProcB )[1] # Un ciclo 'i' abarca 'np' cores simultáneamente
				Fut = Array{Future}( undef, RangosProcB[i][end] ) ; # Arreglo de 'futuros'
				for q in RangosProcB[i] # Recorre los procesadores del rango. Cada core hace una lonja de la matriz de 'Rangos[q+np(i-1)]' columnas
					Fut[ q ] = @spawnat( cores[q], Fill_Matriz_Shell_offdiagonal( K1, eval(:rho1), RangosB[ q + np*(i-1) ], 
					eval(:selv1), eval(:vertex1), eval(:normales1), eval(:selv2), eval(:vertex2), eval(:normales2), TypeNumber ) ) ;
				end
				for q in RangosProcB[i]
					RangoColB = RangosB[ q + np * ( i - 1 ) ] ;
					WKSPC[ 1 : 2*N, 2*N .+ RangoColB ] = fetch( Fut[ q ] )[ :, 1:size( RangoColB )[1] ] ;
					WKSPC[ 1 : 2*N, 2*N + L .+ RangoColB ] = fetch( Fut[ q ] )[ :, size( RangoColB )[1]+1:end] ;
				end
				println("Done with the chunk [B] : ", i ," / ", size( RangosProcB )[1] ) ;
			end
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# Bloque "C" de 2*L x 2*N 
			RangosC = BuildRangos( N, ColPerProc ) ; # 
			RangosProcC = BuildRangosProc( RangosC, np ) ; # Rangos de procesadores [1:np] salvo quizás el último
			@time for i = 1 : size( RangosProcC )[1] # Un ciclo 'i' abarca 'np' cores simultáneamente
				Fut = Array{Future}( undef, RangosProcC[i][end] ) ; # Arreglo de 'futuros'
				for q in RangosProcC[i] # Recorre los procesadores del rango. Cada core hace una lonja de la matriz de 'Rangos[q+np(i-1)]' columnas
					Fut[ q ] = @spawnat( cores[q], Fill_Matriz_Shell_offdiagonal( K1, 
									eval(:rho1), RangosC[ q + np*(i-1) ], 
					eval(:selv2), eval(:vertex2), eval(:normales2), eval(:selv1), 
								eval(:vertex1), eval(:normales1), TypeNumber ) ) ;
				end
				for q in RangosProcC[i]
					RangoColC = RangosC[ q + np * ( i - 1 ) ] ;
					WKSPC[ 2*N+1 : 2*(N+L), RangoColC ] = -fetch( Fut[ q ] )[ :, 1:size( RangoColC )[1] ] ;
					WKSPC[ 2*N+1 : 2*(N+L), N .+ RangoColC ] = -fetch( Fut[ q ] )[ :, size( RangoColC )[1]+1:end] ;
				end
				println("Done with the chunk [C] : ", i ," / ", size( RangosProcC )[1] ) ;
			end
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# Bloque "D" de 2*L x 2*L 
			RangosD = BuildRangos( L, ColPerProc ) ; # 
			RangosProcD = BuildRangosProc( RangosD, np ) ; # Rangos de procesadores [1:np] salvo quizás el último
			@time for i = 1 : size( RangosProcD )[1] # Un ciclo 'i' abarca 'np' cores simultáneamente
				Fut = Array{Future}( undef, RangosProcD[i][end] ) ; # Arreglo de 'futuros'
				for q in RangosProcD[i] # Recorre los procesadores del rango. Cada core hace una lonja de la matriz de 'Rangos[q+np(i-1)]' columnas
					Fut[ q ] = @spawnat( cores[q], Fill_Matriz_Shell_diagonal( K1, 
								K2, eval(:rho1), eval(:rho2),
					RangosD[ q + np*(i-1) ], eval(:selv2), eval(:vertex2), 
								eval(:normales2), TypeNumber ) ) ;
				end
				for q in RangosProcD[i]
					RangoColD = RangosD[ q + np * ( i - 1 ) ] ;
					WKSPC[ 2*N+1 : 2*(N+L), 2*N .+ RangoColD ] = 
						fetch( Fut[ q ] )[ :, 1:size( RangoColD )[1] ] ;
					WKSPC[ 2*N+1 : 2*(N+L), 2*N + L .+ RangoColD ] = 
						fetch( Fut[ q ] )[ :, size( RangoColD )[1]+1:end] ;
				end
				println("Done with the chunk [D] : ", i ," / ", size( RangosProcD )[1] ) ;
			end
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# Vector con condiciones de borde
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			BndValues = Array{ TypeNumber }( undef, 2 * ( N + L ) ) ; # Array local de valores de borde
			BndValues[ : ] = FillBndValues_Shells( kinc, K0, rho0, selv1, vertex1, normales1, L ) ;
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# Resolviendo el sistema 
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# Resolución del sistema para todos las incidencias
			# Se reescriben WKSPC (con autovectores) y BndValues (con la solución)
			LinearAlgebra.LAPACK.gesv!( WKSPC, BndValues ) ;
			# Liberación de recursos
			WKSPC = TypeNumber( 0 ) ;
			@everywhere GC.gc() ;
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# Cálculo del farfield
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# BndValues es ahora la "densidad" solución
			FarField_bs_K[ w ] = FarField_Calculation( K0, -kinc, rho0, selv1, vertex1, 
						normales1, BndValues ) ;
		end
		return FarField_bs_K
	end
	
	# Función que llena un rango de columnas 'rango' = (start:end) (que va de 1 a N) en cada una de las dos submatrices
	# (left,right) de la matriz diagonal de 2*N x 2*N en el caso del sistema shell.
	# Se llenan las columnas dadas por 'rango' en la submatriz left y en la right simultáneamente. Por ello es posible
	# calcular todos los operadores utilizando un único llamado a la función.
	# 
	function Fill_Matriz_Shell_diagonal( K0::Real, K1::Real, rho0::Real, rho1::Real, rango, 
			    selv::Array, vertex::Array, normales::Array, TypeNumber )
		# Constantes y parámetros
		N = size( selv )[1] ; # El bloque tiene size 2N x 2N
		Cols = size( rango )[1] ; # Número de columnas
		Js = collect( rango ) ; # Array con los índices de las columnas entre el 1:N
		SubMat = Array{ TypeNumber }( undef, 2 * N, 2 * Cols ) ; # Son dos submatrices en realidad
		for g = 1 : Cols # Recorro cada una de las columnas con un índice local 'g'
			j = Js[ g ] ; # 'j' indice global (va de 1:N)
			QA, QB, QC = VerticesTriangle( j, selv, vertex ) ;
			Nj = normales[ j, : ] ;
			for i = 1 : N     
 				P = CentroideTriangle( i, selv, vertex ) ;
				Ni = normales[i,:] ; 
				if i == j # Lponel true
					DISKMuller, DISLK0, DISLK1, DISMK0, DISMK1, DISMK0T, DISMK1T = 
						Operadores_Fluid( K0, K1, P, Ni, QA, QB, QC, Nj, true ) ;
				else  # Lponel false
					DISKMuller, DISLK0, DISLK1, DISMK0, DISMK1, DISMK0T, DISMK1T = 
						Operadores_Fluid( K0, K1, P, Ni, QA, QB, QC, Nj, false ) ;
				end               
				SubMat[ i, g ] = 1 / 2 * ( rho0 + rho1 ) * kronecker( i, j ) + 
							( rho0 * DISMK0 - rho1 * DISMK1 ) ;
				SubMat[ i + N, g ] = - DISKMuller ;
				SubMat[ i, Cols + g ] = rho0^2 * DISLK0 - rho1^2 * DISLK1 ;  
				SubMat[ i + N, Cols + g ] = 1 / 2 * ( rho0 + rho1 ) * kronecker( i, j ) -
							( rho0 * DISMK0T - rho1 * DISMK1T ) ;
			end
		end
		return SubMat # Matriz de 2*N x size(rango)
	end

	# Función que llena submatrices off diagonal en el caso del sistema shell. Son submatrices de 'M x ' siendo
	#       M : size de filas (colocación)
	#       R : size de columnas (integración)

	# todo es off diagonal (no hay lponel true)
	function Fill_Matriz_Shell_offdiagonal( K::Real, rho::Real, rango, 
		selvColoc::Array, vertexColoc::Array, normalesColoc::Array,
		selvInteg::Array, vertexInteg::Array, normalesInteg::Array, TypeNumber )
		M = size( selvColoc )[1] ; # Size de la grilla de colocación                    
		R = size( selvInteg )[1] ; # Size de la grilla de integración
		Cols = size( rango )[1] ; # Número de columnas
		Js = collect( rango ) ; # Array con los índices de las columnas entre el 1:R
		SubMat = Array{ TypeNumber }( undef, 2 * M, 2 * Cols ) ; # Son dos submatrices en realidad
 		for g = 1 : Cols # Recorro cada una de las columnas con un índice local 'g'
			j = Js[ g ] ; # 'j' indice global (va de 1:R)
			QA, QB, QC = VerticesTriangle( j, selvInteg, vertexInteg ) ;
			Nj = normalesInteg[ j, : ] ;
			for i = 1 : M     
				P = CentroideTriangle( i, selvColoc, vertexColoc ) ;
				Ni = normalesColoc[i,:] ; 
				DISNK, DISLK, DISMK, DISMKT = Operadores_Shell_offdiagonal( K, P, Ni, QA, QB, QC, Nj ) ;
				SubMat[ i, g ] = -rho * DISMK ;    
				SubMat[ M + i , g ] = DISNK ;   
				SubMat[ i, Cols + g ] = -rho^2 * DISLK ;  
				SubMat[ M + i, Cols + g ] = rho * DISMKT ;    
			end
		end
	return SubMat # Matriz de 2*N x size(rango)   
	end
