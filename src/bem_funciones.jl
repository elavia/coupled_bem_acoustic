# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	
#	Geometrical operations with triangles
#	
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# v. 13/8/20

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Index:

# 	normal3( VECA, VECB, VECC )
# 	area( VECA, VECB, VECC )
# 	VerticesTriangle( q, selv, vertex )
# 	CentroideTriangle( q, selv, vertex )
# 	NormalesTriangle( q::Int, MP::AbstractArray )
# 	VerticesTriangle( q::Int, MP::AbstractArray )
# 	CentroideTriangle( q::Int, MP::AbstractArray )
# 	kronecker( i::Int , j::Int )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	# Función que da la normal a tres vectores. Es la normal externa
	# porque los vectores se esperan en CCW order. Traducción rutina de S. Kirkup
	function normal3( VECA, VECB, VECC )
		VBMVA = VECB - VECA ; 
		VCMVA = VECC - VECA ; 
		CR1 = VBMVA[2]*VCMVA[3] - VBMVA[3]*VCMVA[2] ;
		CR2 = VCMVA[1]*VBMVA[3] - VBMVA[1]*VCMVA[3] ;
		CR3 = VBMVA[1]*VCMVA[2] - VBMVA[2]*VCMVA[1] ;
		ACR = sqrt( CR1*CR1 + CR2*CR2 + CR3*CR3 ) ;
		return [ CR1/ACR, CR2/ACR, CR3/ACR ] ;
	end
	
	# Area del triángulo entre tres vectores
	function area( VECA, VECB, VECC )
		VECD = VECB - VECA ; 
		VECE = VECC - VECA ; 
		RQAQB = norm( VECD )^2 ; 
		RQAQC = norm( VECE )^2 ;
		DP = dot( VECD, VECE ) ;
		return sqrt( RQAQB*RQAQC - DP*DP )/2
	end

	# Función que devuelve los tres vértices A,B,C del triángulo 'q' en la mesh 'selv', 'vertex'
	function VerticesTriangle( q, selv, vertex )
		return vertex[ selv[ q, 1 ], 1:3 ], vertex[ selv[ q, 2 ], 1:3 ], vertex[ selv[ q, 3 ], 1:3 ] ;
	end	

	# Cálculo del centroide del triángulo 'q' dado por las estructuras 'selv' 'vertex'
	# centroide = ( Pa + Pb + Pc ) /3    siendo Pa, Pb, Pc los vértices del triángulo	 
	function CentroideTriangle( q, selv, vertex )
		return ( vertex[ selv[ q, 1 ], : ] + vertex[ selv[ q, 2 ], : ] + vertex[ selv[ q, 3 ], : ] ) / 3 ; 
	end

	# Función que da la normal del triángulo 'q' a partir de una matriz 'MP'
	function NormalesTriangle( q::Int, MP::AbstractArray )
		return  MP[10:12,q];
	end
	
	# Función que da los vértices del triángulo 'q' a partir de una matriz 'MP'
	function VerticesTriangle( q::Int, MP::AbstractArray )
		return MP[1:3,q], MP[4:6,q], MP[7:9,q] ;
	end
	
	# Función que da el centroide del triángulo 'q' a partir de una matriz 'MP'
	function CentroideTriangle( q::Int, MP::AbstractArray )
		return ( MP[ 1:3, q ] .+ MP[ 4:6, q ] .+ MP[ 7:9, q ] ) ./ 3 ;
	end

	# Delta de Kronecker para dos índices 'ij' tal que kronecker(i,j) = 0.0 si i\neq j y 1.0 si i=j
	function kronecker( i::Int , j::Int )
		if i == j
			return 1.0
		else
			return 0.0
		end
	end
