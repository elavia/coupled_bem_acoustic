# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	
#	Julia auxiliary routines
#	
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# v. 13/8/20

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Index:

# 	SinC( x )
#	TS( finf::Number )
#	Linspace( Start::Number, Stop::Number, Size::Number )
# 	BiG( numero::Real )
# 	BiG( numero::Complex )
# 	BiG( numero::Irrational )
# 	BiG( numero::Complex{Bool} )
# 	BiG( numero::String )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	
	function SinC( x )
		return sinc( x/pi ) ;
	end
	
	# Función para calcular el 'TS' a partir de un patrón de campo lejano 'finf' 
	function TS( finf::Number )
		return 20 * log10( abs( finf ) ) ;
	end
	
	function Linspace( Start::Number, Stop::Number, Size::Number )
		return range( Start, stop=Stop, length=Size )
	end
	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
	# Arbitrary Precision Routines
	
	# Diferentes métodos para convertir un número a Bigfloat tomando correctamente decimales
	# Para literales conviene big"nro", que evita la conversión previa a string
	function BiG( numero::Real )
		parse( BigFloat,"$(numero)" );
	end
	function BiG( numero::Complex )
		parse(BigFloat,"$(real(numero))") + im*parse(BigFloat,"$(imag(numero))") ;
	end
	# Aparentemente (para pi, e funciona) la implementación "big" lo convierte del modo correcto
	function BiG( numero::Irrational ) # Para entradas como pi, e
		big( numero ) ; # BiG( 1.0 )*numero ;
	end
	function BiG( numero::Complex{Bool} ) # Para entradas como im
		big"1" * numero ;
	end 
	function BiG( numero::String ) # Cuando mandamos un 'string' para no perder las cifras 
	# en la conversión.
		parse( BigFloat, numero );
	end	
	# Por comprehension
	BiG( numero::AbstractArray ) = [ BiG( n ) for n in numero ] ;
	BiG( numero::Tuple ) = [ BiG( n ) for n in numero ] ;
