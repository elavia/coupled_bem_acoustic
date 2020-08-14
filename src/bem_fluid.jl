# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	
#	Julia fluid scattering high-level routines
#	
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# v. 13/8/20

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Index:

# 	bem_fluid_farfield_bs( K0::Real, K1::Real, rho0::Real, rho1::Real, pext::Array,
# 		selv::Array, vertex::Array, normales::Array, cores::Array, TypeNumber )
# 	bem_fluid_farfield_forward( K0::Real, K1::Real, rho0::Real, rho1::Real, pext::Array,
# 		selv::Array, vertex::Array, normales::Array, cores::Array, TypeNumber )
# 	bem_fluid_farfield_angular( K0::Real, K1::Real, rho0::Real, rho1::Real, pext::Array, 
# 		kinc::Array, selv::Array, vertex::Array, normales::Array, cores::Array, TypeNumber )
# 	bem_fluid_nearfield( K0::Real, K1::Real, rho0::Real, rho1::Real, pext::Array, pint::Array,
# 		kinc::Array, selv::Array, vertex::Array, normales::Array, cores::Array, TypeNumber )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


	# Función que calcula f_\infty de BS para una mesh 'selv, vertex, normales' y direcciones dadas por 'pext' 
	# ( las direcciones de incidencia son -pext ) utilizando proceso distribuido (entre el array de procesadores
	# definido por 'cores'). 
	# No utiliza SharedArray para la matriz WKSPC por sus consabidos problemas 
	#
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Versión para cuerpos fluidos (penetrables)
	function bem_fluid_farfield_bs( K0::Real, K1::Real, rho0::Real, rho1::Real, pext::Array,
				selv::Array, vertex::Array, normales::Array, cores::Array, TypeNumber )
		# Parámetros
		NPE = size( pext )[1] ; # Número de direcciones de observación 
		NSE = size( selv )[1] ; # Número de elementos de superficie
		# Configuración de la ejecución en paralelo en el cluster
		np = size( cores )[1] ; # Número de procesadores (sin el proceso local "1")
		ColPerProc = 64 ; # Número de columnas que llena cada procesador cada vez
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Se copian a todos los workers las estructuras necesarias
		# Este paso hay que hacerlo aquí dentro puesto que cuando evalúo eval(:VAR) dentro
		# de la función en realidad VAR se busca primero en el host remoto y no se toma 
		# del argumento de entrada VAR.
		for i in procesadores # Se copian al procesador 'i' las estructuras necesarias
			SendToProc( i, selv = selv, vertex = vertex, normales = normales, pext = pext ) ;
			SendToProc( i, K0 = K0, K1 = K1, rho0 = rho0, rho1 = rho1 ) ;
		end
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Cálculo de operadores Lk, Mk, Mkt y Müller y llenado de Matriz WKSPC
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Creación de la estructura WKSPC 
		WKSPC = Array{ TypeNumber }( undef, 2*NSE, 2*NSE ) ;
		println( "Filling WKSPC matrix with distributed process ... " ) ;
		# Se divide la mitad de columnas a llenar (NSE) entre el nro de columnas para cada procesador y se
		# construye un vector que tiene rangos de columnas a calcular ([ nro de ColPerProc] salvo quizás
		# el último).
		# Se dividen esos rangos en ciclos de 'np' procesadores. Tendremos 'size( RangosProc )[1]' ciclos de 
		# 'np' procesadores. 'RangosProc=[1:np]' salvo quizás el último ciclo. 'Rangos[q + np * ( i - 1 )]' 
		# es el rango de columnas [j_start:j_end] ( 1 <= j <= NSE ) que hace el procesador 'q' en el ciclo 'i'. 
		# Cada procesador hace dos submatrices 
		Rangos = BuildRangos( NSE, ColPerProc ) ; 
		RangosProc = BuildRangosProc( Rangos, np ) ; # Rangos de procesadores [1:np] salvo quizás el último
		@time for i = 1 : size( RangosProc )[1] # Un ciclo 'i' abarca 'np' cores simultáneamente
			Fut = Array{Future}( undef, RangosProc[i][end] ) ; # Arreglo de 'futuros'
			for q in RangosProc[i] # Recorre los procesadores del rango. 
				# Cada core hace una lonja de la matriz de 'Rangos[q+np(i-1)]' columnas
				Fut[ q ] = @spawnat( cores[q], Fill_Matriz_Fluid( eval(:K0), eval(:K1), 
						eval(:rho0), eval(:rho1), Rangos[ q + np*(i-1) ], eval(:selv), 
						eval(:vertex), eval(:normales), TypeNumber ) ) ;
			end
			for q in RangosProc[i]
				RangoCol = Rangos[ q + np * ( i - 1 ) ] ;
				WKSPC[1:2*NSE,RangoCol] = fetch( Fut[q] )[:,1:size( RangoCol )[1]] ;
				WKSPC[1:2*NSE,NSE .+ RangoCol] = fetch( Fut[q] )[:,size( RangoCol )[1]+1:end] ;
			end
			println("Chunk done : ", i ," / ", size( RangosProc )[1] ) ;
		end
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Vector con condiciones de borde
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		println( "Calculatinb boundary values in parallel ... " ) ;
		BndValues = Array{ TypeNumber }( undef, 2 * NSE, NPE ) ; # Array local de valores de borde
		# Se envía el cálculo a todos los cores de modo asíncrono y se espera a la finalización
		@time @sync @async for j = 1 : NPE 
			BndValues[ :, j ] = fetch( @spawn( FillBndValues_Fluid( -pext[j,:], eval(:K0),
				 	eval(:rho0), eval(:selv), eval(:vertex), eval(:normales) ) ) ) ;
		end
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Resolviendo el sistema 
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		println("Solving the system ... ") ;
		# Resolución del sistema para todos las incidencias
		# Se reescriben WKSPC (con autovectores) y BndValues (con la solución)
		@time LinearAlgebra.LAPACK.gesv!( WKSPC, BndValues ) ;
		# Liberación de recursos
		println("Freeing Memory ... ") ;
		WKSPC = TypeNumber( 0 ) ;
		@everywhere GC.gc() ;
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Cálculo del farfield
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		println("Computing far field ... ") ;
		FarField_bs = Array{ TypeNumber }( undef, NPE ) ; # Array local de farfields
		# No vale la pena paralelizar este cálculo para valores 'NPE' no monstruosos por el overhead.
		# BndValues es ahora la "densidad" solución
		@time for i = 1 : NPE 
			FarField_bs[ i ] = FarField_Calculation( K0, pext[i,:], rho0, 
						selv, vertex, normales, view( BndValues, :, i ) ) ;
		end
		return FarField_bs
	end


	function bem_fluid_farfield_forward( K0::Real, K1::Real, rho0::Real, rho1::Real, pext::Array,
				selv::Array, vertex::Array, normales::Array, cores::Array, TypeNumber )
		# Parámetros
		NPE = size( pext )[1] ; # Número de direcciones de observación 
		NSE = size( selv )[1] ; # Número de elementos de superficie
		# Configuración de la ejecución en paralelo en el cluster
		np = size( cores )[1] ; # Número de procesadores (sin el proceso local "1")
		ColPerProc = 64 ; # Número de columnas que llena cada procesador cada vez
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Se copian a todos los workers las estructuras necesarias
		# Este paso hay que hacerlo aquí dentro puesto que cuando evalúo eval(:VAR) dentro
		# de la función en realidad VAR se busca primero en el host remoto y no se toma 
		# del argumento de entrada VAR.
		for i in procesadores # Se copian al procesador 'i' las estructuras necesarias
			SendToProc( i, selv = selv, vertex = vertex, normales = normales, pext = pext ) ;
			SendToProc( i, K0 = K0, K1 = K1, rho0 = rho0, rho1 = rho1 ) ;
		end
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Cálculo de operadores Lk, Mk, Mkt y Müller y llenado de Matriz WKSPC
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Creación de la estructura WKSPC 
		WKSPC = Array{ TypeNumber }( undef, 2*NSE, 2*NSE ) ;
		println( "Filling WKSPC matrix with distributed process ... " ) ;
		# Se divide la mitad de columnas a llenar (NSE) entre el nro de columnas para cada procesador y se
		# construye un vector que tiene rangos de columnas a calcular ([ nro de ColPerProc] salvo quizás
		# el último).
		# Se dividen esos rangos en ciclos de 'np' procesadores. Tendremos 'size( RangosProc )[1]' ciclos de 
		# 'np' procesadores. 'RangosProc=[1:np]' salvo quizás el último ciclo. 'Rangos[q + np * ( i - 1 )]' 
		# es el rango de columnas [j_start:j_end] ( 1 <= j <= NSE ) que hace el procesador 'q' en el ciclo 'i'. 
		# Cada procesador hace dos submatrices 
		Rangos = BuildRangos( NSE, ColPerProc ) ; 
		RangosProc = BuildRangosProc( Rangos, np ) ; # Rangos de procesadores [1:np] salvo quizás el último
		@time for i = 1 : size( RangosProc )[1] # Un ciclo 'i' abarca 'np' cores simultáneamente
			Fut = Array{Future}( undef, RangosProc[i][end] ) ; # Arreglo de 'futuros'
			for q in RangosProc[i] # Recorre los procesadores del rango. 
				# Cada core hace una lonja de la matriz de 'Rangos[q+np(i-1)]' columnas
				Fut[ q ] = @spawnat( cores[q], Fill_Matriz_Fluid( eval(:K0), eval(:K1), 
						eval(:rho0), eval(:rho1), Rangos[ q + np*(i-1) ], eval(:selv), 
						eval(:vertex), eval(:normales), TypeNumber ) ) ;
			end
			for q in RangosProc[i]
				RangoCol = Rangos[ q + np * ( i - 1 ) ] ;
				WKSPC[1:2*NSE,RangoCol] = fetch( Fut[q] )[:,1:size( RangoCol )[1]] ;
				WKSPC[1:2*NSE,NSE .+ RangoCol] = fetch( Fut[q] )[:,size( RangoCol )[1]+1:end] ;
			end
			println("Chunk done : ", i ," / ", size( RangosProc )[1] ) ;
		end
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Vector con condiciones de borde
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		println( "Calculating the boundary values in parallel ... " ) ;
		BndValues = Array{ TypeNumber }( undef, 2 * NSE, NPE ) ; # Array local de valores de borde
		# Se envía el cálculo a todos los cores de modo asíncrono y se espera a la finalización
		@time @sync @async for j = 1 : NPE 
			# El único cambio respecto de bs es que 'kinc = pext' aquí
			BndValues[ :, j ] = fetch( @spawn( FillBndValues_Fluid( pext[j,:], eval(:K0),
				 	eval(:rho0), eval(:selv), eval(:vertex), eval(:normales) ) ) ) ;
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
		@time for i = 1 : NPE 
			FarField_bs[ i ] = FarField_Calculation( K0, pext[i,:], rho0, 
						selv, vertex, normales, view( BndValues, :, i ) ) ;
		end
		return FarField_bs
	end

	# Función que calcula f_\infty bajo dirección de incidencia 'kinc' para una mesh 'selv, vertex, normales'
	# según las direcciones dadas por 'pext' utilizando proceso distribuido (entre el array de procesadores
	# definido por 'cores'). 
	# No utiliza SharedArray para la matriz WKSPC por sus consabidos problemas 
	#
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Versión para cuerpos fluidos (penetrables)
	function bem_fluid_farfield_angular( K0::Real, K1::Real, rho0::Real, rho1::Real, pext::Array, 
			kinc::Array, selv::Array, vertex::Array, normales::Array, cores::Array, TypeNumber )
		# Parámetros
		NPE = size( pext )[1] ; # Número de direcciones de observación 
		NSE = size( selv )[1] ; # Número de elementos de superficie
		# Configuración de la ejecución en paralelo en el cluster
		np = size( cores )[1] ; # Número de procesadores (sin el proceso local "1")
		ColPerProc = 64 ; # Número de columnas que llena cada procesador cada vez
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Se copian a todos los workers las estructuras necesarias
		# Este paso hay que hacerlo aquí dentro puesto que cuando evalúo eval(:VAR) dentro
		# de la función en realidad VAR se busca primero en el host remoto y no se toma 
		# del argumento de entrada VAR.
		for i in procesadores # Se copian al procesador 'i' las estructuras necesarias
			SendToProc( i, selv = selv, vertex = vertex, normales = normales, pext = pext ) ;
			SendToProc( i, K0 = K0, K1 = K1, rho0 = rho0, rho1 = rho1 ) ;
		end
		# Creación de la estructura WKSPC (Array local -- se llena por proceso distribuido--)
		WKSPC = Array{ TypeNumber }( undef, 2*NSE, 2*NSE ) ;
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Cálculo de operadores Lk, Mk, Mkt y Nk y Llenado de Matriz WKSPC
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		Rangos = BuildRangos( NSE, ColPerProc ) ; 
		RangosProc = BuildRangosProc( Rangos, np ) ; # Rangos de procesadores [1:np] salvo quizás el último
		@time for i = 1 : size( RangosProc )[1] # Un ciclo 'i' abarca 'np' cores simultáneamente
			Fut = Array{Future}( undef, RangosProc[i][end] ) ; # Arreglo de 'futuros'
			for q in RangosProc[i] # Recorre los procesadores del rango. 
				# Cada core hace una lonja de la matriz de 'Rangos[q+np(i-1)]' columnas
				Fut[ q ] = @spawnat( cores[q], Fill_Matriz_Fluid( eval(:K0), eval(:K1), 
						eval(:rho0), eval(:rho1), Rangos[ q + np*(i-1) ], eval(:selv), 
						eval(:vertex), eval(:normales), TypeNumber ) ) ;
			end
			for q in RangosProc[i]
				RangoCol = Rangos[ q + np * ( i - 1 ) ] ;
				WKSPC[1:2*NSE,RangoCol] = fetch( Fut[q] )[:,1:size( RangoCol )[1]] ;
				WKSPC[1:2*NSE,NSE .+ RangoCol] = fetch( Fut[q] )[:,size( RangoCol )[1]+1:end] ;
			end
			println("Chunk done : ", i ," / ", size( RangosProc )[1] ) ;
		end
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Vector con condiciones de borde
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		println( "Calculating the boundary values ... " ) ;
		BndValues = Array{ TypeNumber }( undef, 2 * NSE ) ; # Array local de valores de borde
		for q = 1 : NSE 
			X = CentroideTriangle( q, selv, vertex ) ;
			BndValues[ q ] = -exp( im * K0 * dot( kinc, X ) ) ;
			BndValues[ NSE + q ] = im / rho0 * K0 * dot( kinc, normales[q, :] ) * 
					exp( im * K0 * dot( kinc, X ) ) ;
		end
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Resolviendo el sistema 
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		println("Solving the system ... ") ;
		# Resolución del sistema
		# Se reescriben WKSPC (con autovectores) y BndValues (con la solución)
		@time LinearAlgebra.LAPACK.gesv!( WKSPC, BndValues ) ;
		# Liberación de recursos
		println("Freeing memory ... ") ;
		WKSPC = TypeNumber( 0 ) ;
		#@everywhere GC.gc() ;
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
				eval(:selv), eval(:vertex), eval(:normales), eval(:BndValues) ), pextV ) ;
		return FarField_angular
	end


	# Calcula el nearfield para puntos correspondientes a un mesh de nearfield
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Versión para cuerpos fluidos (penetrables)
	function bem_fluid_nearfield( K0::Real, K1::Real, rho0::Real, rho1::Real, pext::Array, pint::Array,
			kinc::Array, selv::Array, vertex::Array, normales::Array, cores::Array, TypeNumber )
		# Parámetros
		NPE = size( pext )[1] ; # Número de puntos externos de campo
		NPI = size( pint )[1] ; # Número de puntos internos de campo 
		NSE = size( selv )[1] ; # Número de elementos de superficie
		A0 = 1.0 ; # Amplitud de la onda incidente
#		# Configuración de la ejecución en paralelo en el cluster
		np = size( cores )[1] ; # Número de procesadores (sin el proceso local "1")
		ColPerProc = 64 ; # Número de columnas que llena cada procesador cada vez
		# Creación de la estructura WKSPC (Array local -- se llena por proceso distribuido--)
		WKSPC = Array{ TypeNumber }( undef, 2*NSE, 2*NSE ) ;
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Cálculo de operadores Lk, Mk, Mkt y Nk y Llenado de Matriz WKSPC
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		Rangos = BuildRangos( NSE, ColPerProc ) ; 
		RangosProc = BuildRangosProc( Rangos, np ) ; # Rangos de procesadores [1:np] salvo quizás el último
		@time for i = 1 : size( RangosProc )[1] # Un ciclo 'i' abarca 'np' cores simultáneamente
			Fut = Array{Future}( undef, RangosProc[i][end] ) ; # Arreglo de 'futuros'
			for q in RangosProc[i] # Recorre los procesadores del rango. 
				# Cada core hace una lonja de la matriz de 'Rangos[q+np(i-1)]' columnas
				Fut[ q ] = @spawnat( cores[q], Fill_Matriz_Fluid( eval(:K0), eval(:K1), 
						eval(:rho0), eval(:rho1), Rangos[ q + np*(i-1) ], eval(:selv), 
						eval(:vertex), eval(:normales), TypeNumber ) ) ;
			end
			for q in RangosProc[i]
				RangoCol = Rangos[ q + np * ( i - 1 ) ] ;
				WKSPC[1:2*NSE,RangoCol] = fetch( Fut[q] )[:,1:size( RangoCol )[1]] ;
				WKSPC[1:2*NSE,NSE .+ RangoCol] = fetch( Fut[q] )[:,size( RangoCol )[1]+1:end] ;
			end
			println("Chunk done : ", i ," / ", size( RangosProc )[1] ) ;
		end
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Vector con condiciones de borde
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		println( "Calculating the boundary values ... " ) ;
		BndValues = Array{ TypeNumber }( undef, 2 * NSE ) ; # Array local de valores de borde
		for q = 1 : NSE 
			X = CentroideTriangle( q, selv, vertex ) ;
			BndValues[ q ] = -exp( im * K0 * dot( kinc, X ) ) ;
			BndValues[ NSE + q ] = im / rho0 * K0 * dot( kinc, normales[q, :] ) * 
					exp( im * K0 * dot( kinc, X ) ) ;
		end
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Resolviendo el sistema 
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		println("Solving the system ... ") ;
		# Resolución del sistema 
		# Se reescriben WKSPC (con autovectores) y BndValues (con la solución)
		@time LinearAlgebra.LAPACK.gesv!( WKSPC, BndValues ) ;
		# Liberación de recursos
		println("Freeing memory ... ") ;
		WKSPC = TypeNumber( 0 ) ;
		#@everywhere GC.gc() ;
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Cálculo del nearfield
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		println("Computing field in the external / internal points ... ") ;
		# Esta paralelización trivial parece ser conveniente si 'BndValues' no es monstruosa.
		# BndValues es ahora la "densidad" solución
		for i in cores # Se copia a los cores BndValues
			SendToProc( i, BndValues = BndValues ) ;
		end
		pextV = [ pext[i,:] for i=1:NPE ] ; # Se convierte a un vector de incidencias
		# Calculo en paralelo entre los cores
		@time NearFieldExt = pmap( x -> NearField_Calculation( eval(:K0), x, eval(:rho0), 
				eval(:selv), eval(:vertex), eval(:normales), eval(:BndValues) ), pextV ) ;
				
		@time IncField = pmap( x -> IncidentField_Calculation( A0, eval(:K0), eval(:kinc), x ),
				 pextV )	;
				
		pintV = [ pint[i,:] for i=1:NPI ] ; # Se convierte a un vector de incidencias
		# Calculo en paralelo entre los cores
		@time NearFieldInt = pmap( x -> NearField_Calculation( eval(:K1), x, eval(:rho1), 
				eval(:selv), eval(:vertex), eval(:normales), eval(:BndValues) ), pintV ) ;
		return IncField, NearFieldExt, NearFieldInt ;
	end
