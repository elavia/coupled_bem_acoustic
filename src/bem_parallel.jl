# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	
#	BEM scattering low-level routines
#	
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# v. 13/8/20

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Index:

# 	Fill_Matriz_Fluid( K0::Real, K1::Real, rho0::Real, rho1::Real, rango, 
# 		selv::Array, vertex::Array, normales::Array, TypeNumber )
# 
# 	SharedArrays routines:
# 
# 	BuildRangos( Size, ColPerProc )
# 	BuildRangosProc( Rangos::Array, np::Int64 )
# 	FillBndValues_Shells( kinc::Array, K::Real, rho::Real, selv::Array, 
# 		vertex::Array, norms::Array, L::Int64 ; TypeNumber = ComplexF64 )
# 	FillBndValues_Fluid( kinc::Array, K::Real, rho::Real, selv::Array,
# 		vertex::Array, norms::Array )
# 	FillBndValues( BV::AbstractArray, idx::Int, kinc::Array, K::Real, MU::Number, 
# 		selv::AbstractArray, vertex::AbstractArray, normales::AbstractArray )
# 	FillBndValues_Impenetrable( kinc::Array, K::Real, MU::Number, 
# 		selv::Array, vertex::Array, norms::Array ; TypeNumber = ComplexF64 )
# 	IncidentField_Calculation( A0::Real, K::Real, kinc::Array, X::Array )
# 	FarField_Calculation( K::Real, Pext::Array, rho::Real,
# 		selv, vertex, norms, Density ) 
# 	NearField_Calculation( K::Real, Pext::Array, rho::Real,
# 		selv, vertex, norms, Density ) 
# 	FarField_Calculation( K::Real, Pext, BC::String, selv, vertex, normales, Density )
# 
# 	DistributedArrays routines:
# 
# 	Fill_NeuOper_SubMatrix( Indices, K::Real, MU, selv, vertex, normales ;
# 		TypeNumber = ComplexF64, R_changeRule = 1e6 )
# 	Fill_DirOper_SubMatrix( Indices, K::Real, MU, selv, vertex, normales;
# 		TypeNumber = ComplexF64, R_changeRule = 1e6 )
# 	Fill_BndValues_SubMatrix( Indices, K::Real, MU, selv, vertex, normales, 
# 		pext ; TypeNumber = ComplexF64 )
# 	Fill_FarField_Vector( Indices, K, pext, BC, selv, vertex, normales, 
# 		Density ; TypeNumber = ComplexF64 )
# 
# 	Shells specific routines:
# 
# 	NearField_shell_medio1( K::Real, Pext::Array, rho::Real,
# 		selv1, vertex1, norms1, selv2, vertex2, norms2, Density ) 
# 	NearField_shell( K::Real, Pext::Array, rho::Real,
# 		selv, vertex, norms, Density ) 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 	0 . Rutinas para trabajar con Arrays (cluster)

	# Función que llena un rango de columnas 'rango' = (start:end) (que va de 1 a N) en la matriz Fluida de 2*N x 2*N 
	# Se llenan las columnas dadas por 'rango' en la submatriz left y en la right simultáneamente. Por ello es posible
	# calcular todos los operadores utilizando un único llamado a la función.
# 
	function Fill_Matriz_Fluid( K0::Real, K1::Real, rho0::Real, rho1::Real, rango, 
		selv::Array, vertex::Array, normales::Array, TypeNumber )
		# Constantes y parámetros
		N = size( selv )[1] ; # Size del bloque de 2N x 2N
		Cols = size( rango )[1] ; # Número de columnas
		Js = collect( rango ) ; # Array con los índices de las columnas entre el 1:N
		SubMat = Array{ TypeNumber }( undef, 2*N, 2*Cols ) ;
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
				SubMat[ i + N, g ] = -DISKMuller ;
				SubMat[ i, Cols + g ] = rho0^2 * DISLK0 - rho1^2 * DISLK1 ;  
				SubMat[ i + N, Cols + g ] = 1 / 2 * ( rho0 + rho1 ) * kronecker( i, j ) -
							 ( rho0 * DISMK0T - rho1 * DISMK1T ) ;
			end
		end
		return SubMat # Matrices de 2*N x size(rango)
	end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   1 . Rutinas para trabajar con SharedArrays

	# Construcción del array de rangos de columnas para una matriz de 'Size' columnas.
	# Cada elemento es un range = (start:end) donde la cantidad de elementos de cada rango
	# es 'ColPerProc' salvo quizás el último, que tendrá una cantidad menor.
	function BuildRangos( Size, ColPerProc )
		chunks = div( Size, ColPerProc ) ; # Nro entero de chunks 
		resto = Size - chunks * ColPerProc ;
		Rangos = [ ( ColPerProc * ( i - 1 ) + 1 : ColPerProc * i ) for i = 1 : chunks ]
		if resto != 0
			extra = [ (ColPerProc * chunks + 1 : ColPerProc * chunks + resto) ] ;
			Rangos = vcat( Rangos, extra ) ;
		end
		return Rangos
	end

	# Construcción del array de rangos de procesadores para el vector de 'Rangos' entre 'np'
	# procesadores. Cada elementeo es un range = (1:np) salvo quizás el último, que tendrá
	# una cantidad menor.
	function BuildRangosProc( Rangos::Array, np::Int64 )
		ciclos = div( size( Rangos )[ 1 ], np ) ;
		remanente = size( Rangos )[ 1 ] - ciclos * np ;
		CiclosArray = [ (1:np) for i = 1 : ciclos ] ; # Tiene los rangos de procesadores
		if remanente > 0
			CiclosArray = vcat( CiclosArray, [ 1:remanente ] ) ;
		end
		return CiclosArray
	end

	# Se genera el vector BndValues de datos de borde (caso shells) correspondiente a la incidencia 'kinc'
	#
	#  	Bnd[1:N] = - e^( i K k * X )
	#  	Bnd[N+1:2N] = i/rho  k * N e^( i K k * X )
	#  	Bnd[2N+1:2(N+L)] = 0
	#
	# El tamaño es 2*N + 2*L (el doble de los tamaños de las meshes exterior e interior)
	function FillBndValues_Shells( kinc::Array, K::Real, rho::Real, selv::Array, 
		vertex::Array, norms::Array, L::Int64 ; TypeNumber = ComplexF64 )
		N = size( selv )[ 1 ] ;
		Bnd = fill!( Array{ TypeNumber }( undef, 2*N + 2*L ), TypeNumber( 0 ) ) ;
		for q = 1 : N 
			X = CentroideTriangle( q, selv, vertex ) ;
			Bnd[ q ] = - exp( im * K * dot( kinc, X ) ) ;
			Bnd[ q + N ] = im / rho * K * dot( kinc, norms[q,:] ) * exp( im * K * dot( kinc, X ) ) ;
		end
		return Bnd ;
	end
	
	# Se genera el vector BndValues de datos de borde (caso fluido) correspondiente a la incidencia 'kinc'
	#
	#	Bnd[1:N] = - e^( i K k * X )
	#	Bnd[N+1:2N] = i/rho  k * N e^( i K k * X )
	#
	# El tamaño es 2*N (el doble del tamaño del mesh)
	function FillBndValues_Fluid( kinc::Array, K::Real, rho::Real, selv::Array,
		vertex::Array, norms::Array )
		N = size( selv )[ 1 ] ;
		Bnd = Array{ComplexF64}( undef, 2*N ) ;
		for q = 1 : N 
			X = CentroideTriangle( q, selv, vertex ) ;
			Bnd[ q ] = - exp( im * K * dot( kinc, X ) ) ;
			Bnd[ q + N ] = im / rho * K * dot( kinc, norms[q,:] ) * exp( im * K * dot( kinc, X ) ) ;
		end
		return Bnd ;
	end
	
	# Función que genera un vector de valores de borde por cada dirección de incidencia 'kinc'
	# La entrada es la matriz 'BV' de valores de borde, el índice 'idx' correspondiente a la
	# incidencia y varios parámetros más.
	function FillBndValues( BV::AbstractArray, idx::Int, kinc::Array, K::Real, MU::Number, 
			    selv::AbstractArray, vertex::AbstractArray, normales::AbstractArray )
		NSE = size( selv )[ 1 ] ;
		for q = 1 : NSE 
			X = CentroideTriangle( q, selv, vertex ) ;
			BV[ q, idx ] = -( 1.0 + MU * im * K * dot( kinc, normales[q, :] ) ) * exp( im * K * dot( kinc, X ) ) ;
		end
	end
	
	# Se genera el vector BndValues de datos de borde (caso impenetrable) correspondiente a la incidencia 'kinc'
	#
	#  	Bnd[i] = -( 1 + mu i K (Kinc * Ni) )* e^( i K k * X )
	#  	
	function FillBndValues_Impenetrable( kinc::Array, K::Real, MU::Number, 
			selv::Array, vertex::Array, norms::Array ; TypeNumber = ComplexF64 )
		NSE = size( selv )[ 1 ] ;
		Bnd = Array{ TypeNumber }( undef, NSE ) ;
		for q = 1 : NSE 
			X = CentroideTriangle( q, selv, vertex ) ;
			Bnd[ q ] =  -( 1.0 + MU * im * K * dot( kinc, normales[q, :] ) ) * exp( im * K * dot( kinc, X ) ) ;
		end
		return Bnd ;
	end
	
	# Función que calcula el campo incidente de dirección 'kinc', número de onda 'K' y
	# amplitud 'A0' en un punto 'X'
	function IncidentField_Calculation( A0::Real, K::Real, kinc::Array, X::Array )
		return A0 * exp( im * K * dot( kinc, X ) ) ;
	end
	
	# Calcula el f_\infty de farfield de backscattering para la dirección 'Pext' ( = -kinc )
	# con una regla de cuadratura QuadRule (para puntos lponel false) utilizando los valores
	# de la densidad en el borde 'Density' ( = [ Psi Phi ]). Los otros parámetros son el
	# número de onda 'K', la densidad 'rho' del medio y la superficie exterior 'selv', 'vertex',
	# 'norms'
	# Density es un SubArray o Array de size NSE y selv, vertex y norms pueden ser SharedArrays 
	# eventualmente
	function FarField_Calculation( K::Real, Pext::Array, rho::Real,
					selv, vertex, norms, Density ) 
		FF = ComplexF64( 0 ) ;
		QuadRuleOff = QR_cqutm6( ) ; # Quadrature rule for external points
		NSE = size( selv )[1] ;
		for j = 1 : NSE
			QA, QB, QC = VerticesTriangle( j, selv, vertex ) ;
			DISLK = Operador_Lk_farfield( K, Pext, QA, QB, QC, QuadRuleOff )
			DISMK = Operador_Mk_farfield( K, Pext, QA, QB, QC, norms[j,:], QuadRuleOff )
			FF += rho * DISMK * Density[ j ] + rho^2 * DISLK * Density[ NSE + j ] ;
		end
		return FF
	end
	
	# Pext es un punto externo cercano (no una dirección)
	function NearField_Calculation( K::Real, Pext::Array, rho::Real,
					selv, vertex, norms, Density ) 
		NF = ComplexF64( 0 ) ;
		QuadRuleOff = QR_cqutm6( ) ; # Quadrature rule for external points
		NSE = size( selv )[1] ;
		for j = 1 : NSE
			QA, QB, QC = VerticesTriangle( j, selv, vertex ) ;
			DISLK = Operador_Lk( K, Pext, QA, QB, QC, false, QuadRuleOff ) ;
			DISMK = Operador_Mk( K, Pext, QA, QB, QC, norms[j,:], QuadRuleOff ) ;
			NF += rho * DISMK * Density[ j ] + rho^2 * DISLK * Density[ NSE + j ] ;
		end
		return NF
	end

	# Cálculo del f_\infty de backscattering para la dirección 'Pext' ( kinc = -Pext ) en el caso
	# de Dirichlet y Neumann conditions (impenetrables). 
	# Esta función sirve para 'bs' y para 'angular'
	# Density es el vector de BndValues correspondiente a la incidencia Pext
	function FarField_Calculation( K::Real, Pext, BC::String, selv, vertex, normales, Density ) 
		FF = ComplexF64( 0 ) ;
		QuadRuleOff = QR_cquts7() ; # Quadrature rule for external points
		NSE = size( selv )[1] ; 
		if BC == "Dir"
			for j = 1 : NSE
			QA, QB, QC = VerticesTriangle( j, selv, vertex ) ;
			DISLK = Operador_Lk_farfield( K, Pext, QA, QB, QC, QuadRuleOff ) ;
			FF -= DISLK * Density[ j ] ;
			end
			return FF ;    
		elseif BC == "Neu"
			for j = 1 : NSE
			QA, QB, QC = VerticesTriangle( j, selv, vertex ) ;
			DISMK = Operador_Mk_farfield( K, Pext, QA, QB, QC, normales[j,:], QuadRuleOff ) ;
			FF += DISMK * Density[ j ] ;
			end
			return FF ;
		end
	end
	
	
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   2 . Rutinas para trabajar con Distributed Arrays (llenado de matrices)


	# Función que rellena una submatriz (dada por 'Indices') en la matriz de Operadores
	# del sistema a resolver en BEM para el caso de Neumann BC. 
	# Esta es una matriz de NSE x NSE ( nro. de triángulos )
	# Exporta 'Arreglo' en Complex64 (Float32 en real e imaginaria)
	function Fill_NeuOper_SubMatrix( Indices, K::Real, MU, selv, vertex, normales ;
		TypeNumber = ComplexF64, R_changeRule = 1e6 )
		RangeRow = Indices[1] ;
		RangeCol = Indices[2] ;
		QuadRuleOn = QR_cqutm9() ; # diagonal
		M = ( size(RangeRow)[1], size(RangeCol)[1] ) ; # Size de la matriz local
		Arreglo = Array{ TypeNumber }( undef, M ) ;
		for Idx in RangeCol
			QA, QB, QC = VerticesTriangle( Idx, selv, vertex ) ;
			Nj = normales[ Idx, : ] ;
			for i in RangeRow
				P = CentroideTriangle( i, selv, vertex ) ;
				Ni = normales[ i, : ] ;
				QuadRuleOff = QuadRule_selector( P, QA, QB, QC, R_changeRule ) ; # Off diagonal
				if i == Idx
					DISMK, DISNK = Operadores_Neu( K, P, Ni, true, QuadRuleOn, QA, QB, QC, Nj ) ;
					Arreglo[ i - Indices[1][1] + 1, Idx - Indices[2][1] + 1 ] = DISMK + MU * DISNK - 0.5 ; 
				else 
					DISMK, DISNK = Operadores_Neu( K, P, Ni, false, QuadRuleOff, QA, QB, QC, Nj ) ;
					Arreglo[ i - Indices[1][1] + 1, Idx - Indices[2][1] + 1 ] = DISMK + MU * DISNK ;
				end
			end
		end
		return Arreglo
	end
    
	# Función que rellena una submatriz (dada por 'Indices') en la matriz de Operadores
	# del sistema a resolver en BEM para el caso de Dirichlet BC. 
	# Esta es una matriz de NSE x NSE ( nro. de triángulos )
	# Exporta 'Arreglo' en Complex64 (Float32 en real e imaginaria)
	function Fill_DirOper_SubMatrix( Indices, K::Real, MU, selv, vertex, normales;
			TypeNumber = ComplexF64, R_changeRule = 1e6 )
		RangeRow = Indices[1] ;
		RangeCol = Indices[2] ;
		QuadRuleOn = QR_cqutm9() ; # diagonal
		M = ( size(RangeRow)[1], size(RangeCol)[1] ) ; # Size de la matriz local
		Arreglo = Array{ TypeNumber }( undef, M ) ;
		for Idx in RangeCol
			QA, QB, QC = VerticesTriangle( Idx, selv, vertex ) ;
			Nj = normales[ Idx, : ] ;
			for i in RangeRow
				P = CentroideTriangle( i, selv, vertex ) ;
				Ni = normales[ i, : ] ;
				QuadRuleOff = QuadRule_selector( P, QA, QB, QC, R_changeRule ) ; # Off diagonal
				if i == Idx
					DISLK, DISMKT = Operadores_Dir( K, P, Ni, true, QuadRuleOn, QA, QB, QC, Nj ) ;
					Arreglo[ i - Indices[1][1] + 1, Idx - Indices[2][1] + 1 ] = DISLK + MU * DISMKT + MU / 2 ;
				else
					DISLK, DISMKT = Operadores_Dir( K, P, Ni, false, QuadRuleOff, QA, QB, QC, Nj ) ;
					Arreglo[ i - Indices[1][1] + 1, Idx - Indices[2][1] + 1 ] = DISLK + MU * DISMKT ;
				end
			end
		end
		return Arreglo
	end

	# Función que rellena una submatriz (dada por 'Indices') en la matriz de valores en el borde
	# del sistema a resolver en BEM. Esta es una matriz de NSE x NPE ( nro. de triángulos por nro. de direcciones )
	# Exporta 'Arreglo' en Complex64 (Float32 en real e imaginaria)
	function Fill_BndValues_SubMatrix( Indices, K::Real, MU, selv, vertex, normales, pext ; TypeNumber = ComplexF64 )
		RangeRow = Indices[1] ; # Sobre triángulos
		RangeCol = Indices[2] ; # Sobre direcciones de incidencia
		M = ( size(RangeRow)[1], size(RangeCol)[1] ) ; # Size de la matriz local
		Arreglo = Array{ TypeNumber }( undef, M ) ;
		for Idx in RangeCol # Recorro direcciones de incidencia
			kinc = -pext[ Idx, : ] ;
			for i in RangeRow
				X = CentroideTriangle( i, selv, vertex ) ;
				Arreglo[ i - Indices[1][1] + 1, Idx - Indices[2][1] + 1 ] = -( 1.0 + 
					MU * im * K * dot( kinc, normales[ i, : ] ) ) * exp( im * K * dot( kinc, X ) ) ;
			end
		end
		return Arreglo
	end
    
	# Función que rellena el vector de farfield a partir de la matriz 'Density' donde
	#	Density = BndValues[:,Indices]
	# es decir que Density es la submatriz correspondiente a las incidencias que se calcularán, cuya identificación
	# dada por el rango 'Indices' pasa a ser un rango local en cada ejecución. Es decir, por ejemplo :
	#	Density[1:100] = BndValues[:,(245:345,)]
	# La variable 'Indices' tiene la forma genérica :
	#	(245:345,)
	function Fill_FarField_Vector( Indices, K, pext, BC, selv, vertex, normales, Density ; TypeNumber = ComplexF64 )
		RangeRow = Indices[1] ;
		M = size( RangeRow )[1] ; # Tamaño vector local
		Arreglo = Array{ TypeNumber }( undef, M ) ; # Arreglo tiene indice local 1,2,... pero Indices son globales
		for i in RangeRow # i es índice global
			local_idx = i - Indices[1][1] + 1 ;
			Arreglo[ local_idx ] = FarField_Calculation( K, pext[i,:], BC, selv, vertex, normales,
							Density[:,local_idx] ) ;
		end
		return Arreglo
	end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   3 . Rutinas para el caso de shells


	# Función que calcula el campo para el medio intermedio ('rho', 'K') rodeado de una
	# frontera externa '1' y una interna '2'.
	# Pext es un punto externo cercano (no una dirección)
	function NearField_shell_medio1( K::Real, Pext::Array, rho::Real,
					selv1, vertex1, norms1, selv2, vertex2, norms2, Density ) 
		NF = ComplexF64( 0 ) ;
		QuadRuleOff = QR_cqutm6( ) ; # Quadrature rule for external points
		NSE1 = size( selv1 )[1] ;
		NSE2 = size( selv2 )[1] ;
		for j = 1 : NSE1 # Se recorre la mesh 1
			QA, QB, QC = VerticesTriangle( j, selv1, vertex1 ) ;
			DISLK = Operador_Lk( K, Pext, QA, QB, QC, false, QuadRuleOff ) ;
			DISMK = Operador_Mk( K, Pext, QA, QB, QC, norms1[j,:], QuadRuleOff ) ;
			NF += rho * DISMK * Density[ j ] + rho^2 * DISLK * Density[ NSE1 + j ] ;
		end
		for j = 1 : NSE2 # Se recorre la mesh 2
			QA, QB, QC = VerticesTriangle( j, selv2, vertex2 ) ;
			DISLK = Operador_Lk( K, Pext, QA, QB, QC, false, QuadRuleOff ) ;
			DISMK = Operador_Mk( K, Pext, QA, QB, QC, norms2[j,:], QuadRuleOff ) ;
			NF += rho * DISMK * Density[ j + 2 * NSE1 ] + rho^2 * DISLK * Density[ 2 * NSE1 + NSE2 + j ] ;
		end
		return NF
	end

	# Función para calcular el campo cercano en un medio 'rho', 'K' que está
	# acotado por una única frontera. Tener cuidado de la 'density', debería
	# tener tamaño 2*size(selv) (es la parte de la densidad correspondiente
	# a la integración sobre la malla descripta por 'selv', 'vertex')
	function NearField_shell( K::Real, Pext::Array, rho::Real,
					selv, vertex, norms, Density ) 
		NF = ComplexF64( 0 ) ;
		QuadRuleOff = QR_cqutm6( ) ; # Quadrature rule for external points
		NSE = size( selv )[1] ;
		for j = 1 : NSE
			QA, QB, QC = VerticesTriangle( j, selv, vertex ) ;
			DISLK = Operador_Lk( K, Pext, QA, QB, QC, false, QuadRuleOff ) ;
			DISMK = Operador_Mk( K, Pext, QA, QB, QC, norms[j,:], QuadRuleOff ) ;
			NF += rho * DISMK * Density[ j ] + rho^2 * DISLK * Density[ NSE + j ] ;
		end
		return NF
	end
