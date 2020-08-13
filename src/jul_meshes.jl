# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	
#	Julia mesh routines
#	
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# v. 13/8/20

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Index:

# 	GetMesh_fast( file::String )
# 	CleanRepetidos( VertList::AbstractArray )
# 	GetIdxVertex( VertList::AbstractArray, VertList1::AbstractArray, V::AbstractArray )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	# Read a STL mesh in 'file' and return the structures:
	#	normales, selv, vertex
	#
	function GetMesh_fast( file::String )
		# Lectura del archivo .stl #
		if file[end-3:end] != ".stl"
			println("Error : Not a STL file?") ;
		end
		line = 2
		valorN = zeros(Float64,1,3) ;
		# Primera apertura para contar las líneas
		f = open(  file ) ;
		NL = countlines( f ) ;
		close( f ) ;
		NE = Int64( (NL-2)/7 ); # Se descuentan las 'solid' y 'endsolid'
		indicesCaras = zeros( Int64, NE, 3 ) ;
		Normales = zeros( Float64, NE, 3 ) ;
		valorA = zeros(Float64,1,3) ;
		valorB = zeros(Float64,1,3) ;
		valorC = zeros(Float64,1,3) ;
		nro_vert_A = Int64( 0 ) ;
		nro_vert_B = Int64( 0 ) ;
		nro_vert_C = Int64( 0 ) ;
		VerticesLista = Array{Array{Float64,2}}(undef,0) ;
		# Segunda apertura del archivo para construir la lista de vértices y normales
		f = open( file ) ;
		if ( split( readline( f ) )[ 1 ] != "solid" ) # Leo 'solid'
			println( "Error: No 'solid' line found." ) ;
		end
		cara = 1 ;
		# Se llena la lista de vértices sin chequear los repetidos
		while !eof( f ) && cara <= NE # Se lee mientras no se llega al fondo
			subline = 1 ;
			while subline < 8
				linea = readline( f ) ; # Se lee una línea #line += 1 ;
				if subline == 1 # Normal del elemento
					valorN = hcat( parse.( Float64, split( linea )[ 3: end] )... ) ;
				end
				if subline == 3 # Vértice 1
					valorA = hcat( parse.( Float64, split( linea )[ 2: end] )... ) ;
					push!( VerticesLista, valorA ) ;
				end
				if subline == 4 # Vértice 2
					valorB = hcat( parse.( Float64, split( linea )[ 2: end] )... ) ;
					push!( VerticesLista, valorB ) ;
				end
				if subline == 5 # Vértice 3
					valorC = hcat( parse.( Float64, split( linea )[ 2: end] )... ) ;
					push!( VerticesLista, valorC ) ;
				end
				subline += 1 ;
			end
			Normales[ cara, : ] = valorN; # Normales de los elementos
			cara += 1 ;
		end
		close( f ) ;
		VerticesListaClean = CleanRepetidos( VerticesLista ) ; # Cleaning duplicates
		VLCcolumn1 = vcat(VerticesListaClean...)[:,1] ; # Convert to matrix
		println("Cleaning") ;
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Third opening of the file
		f = open( file ) ;
		if ( split( readline( f ) )[ 1 ] != "solid" ) # Leo 'solid'
			println( "Error: No line 'solid' found." ) ;
		end
		cara = 1 ;
		# The vertex of each element are identified
		while !eof( f ) && cara <= NE # Se lee mientras no se llega al fondo
			subline = 1 ;
			while subline < 8
				linea = readline( f ) ; # Se lee una línea #line += 1 ;
				if subline == 3 # Vértice 1
					valorA = hcat( parse.( Float64, split( linea )[ 2: end] )... ) ;
					nro_vert_A = GetIdxVertex( VerticesListaClean, VLCcolumn1, valorA ) ;
				end
				if subline == 4 # Vértice 2
					valorB = hcat( parse.( Float64, split( linea )[ 2: end] )... ) ;
					nro_vert_B = GetIdxVertex( VerticesListaClean, VLCcolumn1, valorB ) ;
				end
				if subline == 5 # Vértice 3
					valorC = hcat( parse.( Float64, split( linea )[ 2: end] )... ) ;
					nro_vert_C = GetIdxVertex( VerticesListaClean, VLCcolumn1, valorC ) ;
				end
				subline += 1 ;
			end
			indicesCaras[ cara, : ] = [ nro_vert_A, nro_vert_B, nro_vert_C ];
			cara += 1 ;
			if mod( cara, 10000 ) == 0
				println("Faces :", cara )
			end
		end
		close( f ) ;
		# Cierre del archivo
		println("Number of: Triangles / Vertex: ", size(indicesCaras)[1]," / ",size(VerticesListaClean)[1] ) ;
		# Falta la parte de chequeo de orientación y normales
		return Normales, indicesCaras, vcat( VerticesListaClean... )
	end
	
	# Take a vertex list 'VertList' (with duplicate vertex) and make a cleaning,
	# exporting a list of only the non duplicate vertex.
	function CleanRepetidos( VertList::AbstractArray )
		epsil = 1E-12 ; # Hardcoded precision (maybe 0 works anyway)
		VL = sortslices( vcat(VertList...), dims=1 ) ; # Convert to matrix and sort by x
		VLC = Array{ Array{ Float64,2 } }( undef, 0 ) ; # Vertex list 
		push!( VLC, hcat( VL[1,:] )' ) ; 
		for i = 2 : size( VL )[1] # Cycle in the list
			if norm( VL[i,:] - VL[i-1,:] ) < epsil # If it is the same that the predecessor
				; # Already in the list
			else # 
				push!( VLC, hcat( VL[i,:] )' ) ; 
			end 
		end 
		return VLC
	end
	
	# Return the numerical 'id' inside of the list 'VertList' of a given
	# vertex 'V' if the input is the vector 'VertList1' that corresponds
	# to the first column in 'VertList'.
	# 'VertList' is sorted by the first component of the vertex (the "x" component).
	function GetIdxVertex( VertList::AbstractArray, VertList1::AbstractArray,
		V::AbstractArray )
		v1 = V[1] ; # Component 'x' from V
		L = length( VertList1 ) ;
		idx = 1 ;
		while v1 > VertList1[idx] # Search for the first ocurrence
			idx += 1 ;
		end
		first = idx ;
		final = idx ; # Default
		while ( final < L && v1 == VertList1[final+1] )
			final +=1 ;
		end
		for j = first:final 
			if norm( V - VertList[j]) == 0
				return j
			end 
		end 
		return 0 # If error then zero
	end
