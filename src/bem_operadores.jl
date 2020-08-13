# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	
#	Operator routines
#	
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# v. 13/8/20

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Index:

#	Punto_Q( q, QA, QuadRule, QBMQA, QCMQA )
#	Operador_Lk( K::Real, P::Array, QA::Array, QB::Array, QC::Array, Lponel::Bool, QuadRule::Array )
#	Operador_Lk_farfield( K::Real, p::Array, QA::Array, QB::Array, QC::Array, QuadRule::Array )
#	Operador_L0( P::Array, QA::Array, QB::Array, QC::Array, QBMQA::Array, QCMQA::Array, QCMQB::Array )
#	Integrando_Lk_regular( k::Real, r::Real )
#	Integrando_Lk_robust( k::Real, r::Real )
#	Operador_Mkt( K::Real, P::Array, Np::Array, QA::Array, QB::Array, QC::Array, QuadRule::Array )
#	Operador_Mk( K::Real, P::Array, QA::Array, QB::Array, QC::Array, Nq::Array, QuadRule::Array )
#	Operador_Mk_farfield( K::Real, p::Array, QA::Array, QB::Array, QC::Array, normal_Q::Array, QuadRule::Array )
#	Integrando_Mk_robust_old( k::Real, r::Real )
#	Integrando_Mk_robust( k::Real, r::Real )
#	Operador_Nk_old( K::Real, P::Array, Np::Array, QA::Array, QB::Array, QC::Array, Nq::Array, 
#					Lponel::Bool, QuadRule::Array )
#	Operador_N0_L0( K::Real, P::Array, QA::Array, QB::Array, QC::Array, QBMQA::Array, QCMQA::Array, QCMQB::Array )
#	Operador_Nk( K::Real, P::Array, Np::Array, QA::Array, QB::Array, QC::Array,
#					 Nq::Array, Lponel::Bool, QuadRule::Array )
#	Integrando_Nk_regular( k::Real, z::Real, dotrNprNq::Real, dotNpNq::Real )
#	Integrando_Nk_robust( k::Real, z::Real, dotrNprNq::Real, dotNpNq::Real )
#	Operador_Muller( k0::Real, k1::Real, P::Array, Np::Array, QA::Array, QB::Array, QC::Array, 
#			Nq::Array, Lponel::Bool, QuadRule::Array )
#	SumaTerminosMuller1( r::Real, k0::Real, k1::Real, M::Int64 )
#	SumaTerminosMuller1( r::Real, k0::Real, k1::Real )
#	SumaTerminosMuller2( r::Real, k0::Real, k1::Real, M::Int64 )
#	SumaTerminosMuller2( r::Real, k0::Real, k1::Real )
#	SumaTerminosMuller2_bis( r::Real, k0::Real, k1::Real, z0::Real, z1::Real, M::Int64 )
#	Integrando_Muller_robust( k0::Real, k1::Real, R::Real, RNPRNQ::Real, DNPNQ::Real )
#	Integrando_Muller_robust( k0::Real, k1::Real, R::Real, z0::Real, z1::Real, RNPRNQ::Real, DNPNQ::Real )
#	Operadores_Fluid( k0::Real, k1::Real, P::Array, Np::Array, QA::Array, QB::Array, QC::Array, 
#			Nq::Array, Lponel::Bool ; TypeNumber = ComplexF64, RadioCambioQuadRule = 0.25 )
#	QuadRule_selector( P, QA, QB, QC, Radio )
#	Operadores_Shell_offdiagonal( K::Real, P::Array, Np::Array, QA::Array, QB::Array, QC::Array, 
#			Nq::Array ; TypeNumber = ComplexF64, RadioCambioQuadRule = 0.25 )
#	Operadores_Dir( K::Real, P::Array, Np::Array, Lponel::Bool, QuadRule::Array, 
#			QA::Array, QB::Array, QC::Array, Nq::Array ; TypeNumber = ComplexF64 )
#	Operadores_Neu( K::Real, P::Array, Np::Array, Lponel::Bool, QuadRule::Array, 
#			QA::Array, QB::Array, QC::Array, Nq::Array ; TypeNumber = ComplexF64 )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	Rutinas accesorias comunes para todos los operadores
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	# Función que da el punto 'q'-esimo de una regla de cuadratura 'QuadRule' para un triángulo
	# con vértice principal 'QA' siendo 'QBMQA' y 'QCMQA' los vectores a los otros vértices.
	# Esta implementación es la más rápida hasta ahora
	function Punto_Q( q, QA, QuadRule, QBMQA, QCMQA )
		return [ QA[ 1 ] + QuadRule[ q, 1 ] * QBMQA[ 1 ] + QuadRule[ q, 2 ] * QCMQA[ 1 ],
			QA[ 2 ] + QuadRule[ q, 1 ] * QBMQA[ 2 ] + QuadRule[ q, 2 ] * QCMQA[ 2 ],
			QA[ 3 ] + QuadRule[ q, 1 ] * QBMQA[ 3 ] + QuadRule[ q, 2 ] * QCMQA[ 3 ] ] ;
	end	

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	OPERADOR Lk
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Función que da el operador discreto 'Lk' para el número de onda 'K' en el triángulo con
	# vértices 'QA', 'QB', 'QC' y centroide 'P'. Si 'P' pertenece al triángulo 'Q' es 'Lponel' = true,
	# de lo contrario es false. La regla de cuadratura es una matriz [ x y w ]: 
	function Operador_Lk( K::Real, P::Array, QA::Array, QB::Array, QC::Array, Lponel::Bool, QuadRule::Array )
		# WARNING : La QuadRule cambiará según Lponel. Esto debe asegurarse desde fuera	
		if Lponel # Hay que regularizar el integrando puesto que P pertenece al triángulo Q 
			QBMQA = QB - QA ; QCMQA = QC - QA ; QCMQB = QC - QB ;
			# El integrando será : Lk - L0, entonces pre-sumo L0 a Lk
			Qarea = area( QA, QB, QC ) ; 
			Lk = Operador_L0( P, QA, QB, QC, QBMQA, QCMQA, QCMQB ) / Qarea ;
			for q = 1 : size( QuadRule )[1]
				R = norm( P - Punto_Q( q, QA, QuadRule, QBMQA, QCMQA ) ) ;
				@inbounds Lk += QuadRule[ q, 3 ] * Integrando_Lk_regular( K, R ) ; # ( FPG - FPG0 ) ;
			end
			return Qarea / ( 4*pi ) * Lk ;
		else # "No Lponel" : No hay que regularizar 
			Lk = zero( ComplexF64 ) ;
			QBMQA = QB - QA ; QCMQA = QC - QA ;
			for q = 1 : size( QuadRule )[1]
				R = norm( P - Punto_Q( q, QA, QuadRule, QBMQA, QCMQA ) ) ;
				@inbounds Lk += QuadRule[ q, 3 ] * Integrando_Lk_robust( K, R ) ; # exp( im * K * R ) / ( K * R ) ; # FPG
			end
			return  area( QA, QB, QC ) / ( 4*pi ) * Lk ;
		end
	end
    
	# Función que da el operador discreto 'Lk' para el número de onda 'K' en el triángulo con vértices 
	# 'QA', 'QB', 'QC' y en la dirección dada por 'p' (farfield).
	# La regla de cuadratura es una matriz [ x y w ]: 
	function Operador_Lk_farfield( K::Real, p::Array, QA::Array, QB::Array, QC::Array, QuadRule::Array )
		Lk_ff = ComplexF64( 0 ) ;
		QBMQA = QB - QA ; # 
		QCMQA = QC - QA ; # 
		for q = 1 : size( QuadRule )[1]
			Q = Punto_Q( q, QA, QuadRule, QBMQA, QCMQA ) ;
			Lk_ff += QuadRule[ q, 3 ] * exp( - im * K * dot( p, Q ) ) ;
		end
		return  area( QA, QB, QC ) / ( 4*pi ) * Lk_ff ;
	end

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Rutinas accesorias para el operador 'Lk'

	# Le saqué la división sobre la QAREA
	function Operador_L0( P::Array, QA::Array, QB::Array, QC::Array, QBMQA::Array, QCMQA::Array, QCMQB::Array )
		L0 = ComplexF64( 0 ) ;
		RQAP = norm( P - QA ) ; 
		RQBP = norm( P - QB ) ;
		RQCP = norm( P - QC ) ;
		AR0 = [ RQAP, RQBP , RQCP ] ;
		ARA = [ RQBP, RQCP, RQAP ] ; # podría ser así : ARA = [ AR0[2], AR0[3], AR0[1] ] ;
		AOPP = [ norm( QBMQA ), norm( QCMQB ), norm( QCMQA ) ] ; # AOPP = [ RQAQB, RQBQC, RQCQA ] ;
		for i = 1 : 3 
			R0 = AR0[ i ] ;
			RA = ARA[ i ] ;
			OPP = AOPP[ i ] ;
			if R0 < RA 
				TEMP = RA ;
				RA = R0 ;
				R0 = TEMP ;
			end # END IF
			A = acos( ( RA * RA + R0 * R0 - OPP * OPP ) / 2 / RA / R0 ) ;
			B = atan( RA * sin( A ) / ( R0 - RA * cos( A ) ) ) ; 
			L0 += ( R0 * sin( B ) * ( log( tan( ( B + A ) / 2 ) ) - log( tan( B / 2 ) ) ) ) ; 
		end 
		return L0
	end

	# Función que da un integrando mejor comportado para '(exp( im*K*R ) - 1 )/R' cuando R es pequeño
	# Es para elementos en la diagonal. Integrando regularizado.
	function Integrando_Lk_regular( k::Real, r::Real )
		z = im * k * r ;
		if abs( z ) > 1E-5
			return im * k * ( exp( z ) - 1 ) / z ; #  ( exp( z ) - 1 ) / r
		else
			return im * k * ( 1 + z / 2 + z^2 / 6 + z^3 / 24 + z^4 / 120 ) ;
		end
	end	

	# Integrando del operador 'Lk' fuera de la diagonal
	function Integrando_Lk_robust( k::Real, r::Real )
		z = k * r ;
		if abs( z ) > 1E-4
			return k * exp( im * z ) /  z ;
		else
			K = BiG( k ) ;
			Z = K * BiG( r ) ;
			return ComplexF64( K * ( cos( Z ) / Z + im * SinC( Z ) ) ) ;
		end
	end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 	OPERADORES Mkt y Mk
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Función que da el operador discreto 'Mkt' para el número de onda 'K' integrando en el triángulo con
	# vértices 'QA', 'QB', 'QC' con respecto al triángulo de centroide 'P'. Si 'P' pertenece al triángulo 'Q' el operador da 
	# exactamente 0 allí (la normal del triángulo es ortogonal al vector R = P - Q ).
	# La regla de cuadratura es una matriz [ x y w ]:
	function Operador_Mkt( K::Real, P::Array, Np::Array, QA::Array, QB::Array, QC::Array, QuadRule::Array )
		Mkt = ComplexF64( 0 ) ;
		QBMQA = QB - QA ; # 
		QCMQA = QC - QA ; # QCMQB = QC - QB ; # 
		for q = 1 : size( QuadRule )[1]
			Q = Punto_Q( q, QA, QuadRule, QBMQA, QCMQA ) ;
			R = norm( P - Q ) ;
			Mkt += QuadRule[ q, 3 ] * Integrando_Mk_robust( K, R ) * dot( P - Q, Np )  ;
		end
		return  area( QA, QB, QC ) / ( 4*pi ) * Mkt ;
	end

	# Función que da el operador discreto 'Mk' para el número de onda 'K' integrando en el triángulo con
	# vértices 'QA', 'QB', 'QC' y normal 'Nq' con respecto al triángulo de centroide 'P'. 
	# Si 'P' pertenece al triángulo 'Q' el operador da exactamente 0 allí (la normal del triángulo es 
	# ortogonal al vector R = P - Q ); entonces NO TIENE TÉRMINO EN LA DIAGONAL para triángulos planos.
	# La regla de cuadratura es una matriz [ x y w ]: 
	function Operador_Mk( K::Real, P::Array, QA::Array, QB::Array, QC::Array, Nq::Array, QuadRule::Array )
		Mk = ComplexF64( 0 ) ;
		QBMQA = QB - QA ; # 
		QCMQA = QC - QA ; #  QCMQB = QC - QB ; # 
		for q = 1 : size( QuadRule )[1]
			Q = Punto_Q( q, QA, QuadRule, QBMQA, QCMQA ) ;
			R = norm( P - Q ) ;
			Mk += QuadRule[ q, 3 ] * Integrando_Mk_robust( K, R ) * -dot( P - Q, Nq )  ; # WFPGR*RNQ ; # FPGR
		end
		return  area( QA, QB, QC ) / ( 4*pi ) * Mk ;
	end

	# Función que da el operador discreto 'Mk' para el número de onda 'K' en el triángulo con vértices 
	# 'QA', 'QB', 'QC' y en la dirección dada por 'p' (farfield).
	# La regla de cuadratura es una matriz [ x y w ]: 
	function Operador_Mk_farfield( K::Real, p::Array, QA::Array, QB::Array, QC::Array, normal_Q::Array, QuadRule::Array )
		Mk_ff = ComplexF64( 0 ) ;
		QBMQA = QB - QA ;  
		QCMQA = QC - QA ; 
		pdotN =  dot( p, normal_Q ) ;
		for q = 1 : size( QuadRule )[1]
			Q = Punto_Q( q, QA, QuadRule, QBMQA, QCMQA ) ;
			Mk_ff += QuadRule[ q, 3 ] * -im * K * exp( - im * K * dot( p, Q ) ) * pdotN ; 
		end
		return  area( QA, QB, QC ) / ( 4*pi ) * Mk_ff ;
	end

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# 	Rutinas accesorias para los operadores 'Mk' y 'Mkt'

	# Función que da un integrando mejor comportado para 'exp( im*K*R )/R^2*( im*K*R - 1 )' cuando R es pequeño
	# Sirve para regularizar los operadores M_k y M_kt (fuera de la diagonal)
	function Integrando_Mk_robust_old( k::Real, r::Real )
		z = k * r ;
		if abs( z ) > 1E-5
			return k^2 * exp( im * z ) / z^2 * ( im * z - 1 ) ;
		else
			return k^2 * ( - 1 / z^2 - 1 / 2 - im * z / 3 + z^2 / 8 + im * z^3 / 30 - z^4 / 144 - im * z^5 / 840 ) ;
		end
	end

	# Función que da un integrando mejor comportado para 'exp( im*K*R )/R^3*( im*K*R - 1 )' cuando R es pequeño
	# Sirve para regularizar los operadores M_k y M_kt (fuera de la diagonal)
	# A diferencia del anterior, hemos multiplicado por un factor 1/R
	function Integrando_Mk_robust( k::Real, r::Real )
		z = k * r ;
		if abs( z ) > 1E-3
			return k^3 * exp( im * z ) / z^3 * ( im * z - 1 ) ;
		else
			return k^3 * ( - 1 / z^3 - 1/( 2*z ) - im / 3 + z / 8 + im * z^2 / 30 - z^3 / 144 - im * z^4 / 840 + z^5 / 5760 ) ;
		end
	end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 	OPERADOR Nk
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Función que da el operador discreto 'Nk' para el número de onda 'K' integrando en el triángulo de
	# vértices 'QA', 'QB', 'QC' y normal 'Nq' con respecto al triángulo de centroide 'P' y normal 'Np'.
	# Si 'P' pertenece al triángulo 'Q' es 'Lponel' = true, de lo contrario es false. 
	# La regla de cuadratura es una matriz [ x y w ]: 
	function Operador_Nk_old( K::Real, P::Array, Np::Array, QA::Array, QB::Array, QC::Array, Nq::Array, 
					Lponel::Bool, QuadRule::Array )
		# WARNING : La QuadRule cambiará según Lponel. Esto debe asegurarse desde fuera	
		Nk = ComplexF64( 0 ) ;
		QBMQA = QB - QA ;  
		QCMQA = QC - QA ;  
		QCMQB = QC - QB ; 
		Qarea = area( QA, QB, QC ) ; 
		DNPNQ = dot( Np, Nq ) ;
		if Lponel # Hay que regularizar el integrando puesto que P pertenece al triángulo Q 
			# El integrando será : Nk - N0 - k^2/2*L0, entonces pre-sumo N0 + k^2/2*L0 a Nk
			Nk += Operador_N0_L0( K, P, QA, QB, QC, QBMQA, QCMQA, QCMQB ) / Qarea ;
			for q = 1 : size( QuadRule )[1]
				Q = Punto_Q( q, QA, QuadRule, QBMQA, QCMQA ) ;
				R = norm( P - Q ) ;
				E = exp( im * K * R ) ;
				FPGR = E / R^2 * ( im * K * R - 1 ) ; 
				FPGRR = E / R^3 * ( 2 - 2 * im * K * R - K^2 * R^2 ) ;
				RNP = dot( P - Q, Np ) / R ;
				RNQ = -dot( P - Q, Nq ) / R ;
				RNPRNQ = RNP * RNQ ;
				RNPNQ = -( DNPNQ + RNPRNQ ) / R ;
				Nk += QuadRule[ q, 3 ] * ( (FPGR - (- 1 / R^2 ) ) * RNPNQ + ( FPGRR - 2 / R^3 ) * RNPRNQ - K^2/2 * 1 / R ) ;
			end
		else # "No Lponel" : No hay que regularizar 
			for q = 1 : size( QuadRule )[1]
				Q = Punto_Q( q, QA, QuadRule, QBMQA, QCMQA ) ;
				R = norm( P - Q ) ;
				E = exp( im * K * R ) ;
				FPGR = E / R^2 * ( im * K * R - 1 ) ; 
				FPGRR = E / R^3 * ( 2 - 2 * im * K * R - K^2 * R^2 ) ;
				RNP = dot( P - Q, Np ) / R ;
				RNQ = -dot( P - Q, Nq ) / R ;
				RNPRNQ = RNP * RNQ ;
				RNPNQ = -( DNPNQ + RNPRNQ ) / R ;
				Nk += QuadRule[ q, 3 ] * ( FPGR * RNPNQ + FPGRR * RNPRNQ ) ;
			end
		end
		return  Qarea / ( 4*pi ) * Nk ;
	end

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# 	Rutinas accesorias para el operador 'Nk'

	# Función que da los valores de los operadores N0 + k^2/2 * L0 para la regularización
	function Operador_N0_L0( K::Real, P::Array, QA::Array, QB::Array, QC::Array, 
		QBMQA::Array, QCMQA::Array, QCMQB::Array )
		L0 = Float64( 0 ) ;
		N0 = Float64( 0 ) ;
		RQAP = norm( P - QA ) ; 
		RQBP = norm( P - QB ) ;
		RQCP = norm( P - QC ) ;
		AR0 = [ RQAP, RQBP , RQCP ] ;
		ARA = [ RQBP, RQCP, RQAP ] ; # podría ser así : ARA = [ AR0[2], AR0[3], AR0[1] ] ;
		AOPP = [ norm( QBMQA ), norm( QCMQB ), norm( QCMQA ) ] ; # AOPP = [ RQAQB, RQBQC, RQCQA ] ;
		for i = 1 : 3 
			R0 = AR0[ i ] ;
			RA = ARA[ i ] ;
			OPP = AOPP[ i ] ;
			if R0 < RA 
				TEMP = RA ;
				RA = R0 ;
				R0 = TEMP ;
			end # END IF
			A = acos( ( RA * RA + R0 * R0 - OPP * OPP ) / 2 / RA / R0 ) ;
			B = atan( RA * sin( A ) / ( R0 - RA * cos( A ) ) ) ; 
			L0 += ( R0 * sin( B ) * ( log( tan( ( B + A ) / 2 ) ) - log( tan( B / 2 ) ) ) ) ; 
			N0 += ( ( cos(B+A) - cos(B) ) / R0 / sin( B ) ) ;
		end 
		return N0 + K^2/2 * L0
	end
 
	function Operador_Nk( K::Real, P::Array, Np::Array, QA::Array, QB::Array, 
			QC::Array, Nq::Array, Lponel::Bool, QuadRule::Array )
		# WARNING : La QuadRule cambiará según Lponel. Esto debe asegurarse desde fuera	
		Nk = ComplexF64( 0 ) ;
		QBMQA = QB - QA ;  
		QCMQA = QC - QA ;  
		QCMQB = QC - QB ; 
		Qarea = area( QA, QB, QC ) ; 
		DNPNQ = dot( Np, Nq ) ;
		if Lponel # Hay que regularizar el integrando puesto que P pertenece al triángulo Q 
			# El integrando será : Nk - N0 - k^2/2*L0, entonces pre-sumo N0 + k^2/2*L0 a Nk
			Nk += Operador_N0_L0( K, P, QA, QB, QC, QBMQA, QCMQA, QCMQB ) / Qarea ;
			for q = 1 : size( QuadRule )[1]
				Q = Punto_Q( q, QA, QuadRule, QBMQA, QCMQA ) ;
				R = norm( P - Q ) ;
		r = ( P - Q ) / R ;
		z = K * R ;
		dotrNprNq = dot( r, Np ) * dot( r, Nq ) ;
		Nk += QuadRule[ q, 3 ] * Integrando_Nk_regular(K, z, dotrNprNq, DNPNQ ) ;
			end
		else # "No Lponel" : No hay que regularizar 
			for q = 1 : size( QuadRule )[1]
				Q = Punto_Q( q, QA, QuadRule, QBMQA, QCMQA ) ;
				R = norm( P - Q ) ;
		r = ( P - Q ) / R ;
		z = K * R ;
		dotrNprNq = dot( r, Np ) * dot( r, Nq ) ;
		Nk += QuadRule[ q, 3 ] * Integrando_Nk_robust( K, z, dotrNprNq, DNPNQ ) ;
			end
		end
		return  Qarea / ( 4*pi ) * Nk ;
	end
    
    
	# Función que da un integrando mejor comportado para Nk cuando 'z' = 'K'*'R' es pequeño
	# Es para elementos en la diagonal. Integrando regularizado. No se regulariza teniendo en cuenta
	# el producto escalar dot(r,Nx) el cual será en esos casos tendiendo a nulo.
	# Cuando la norma de z es pequeña se utiliza una aproximación hasta orden 5
	function Integrando_Nk_regular( k::Real, z::Real, dotrNprNq::Real, dotNpNq::Real )
		if abs( z ) > 1E-2
			return  ( k^3 / z^3 ) * ( ( exp( im * z ) * ( 1 - im * z ) - 1 ) * dotNpNq + 
				( exp( im * z ) * ( 3 * im * z - 3 + z^2 ) + 3 ) * dotrNprNq - 1 / 2 * z^2 )
		else
			z2 = z*z ;
			z3 = z*z2 ;
			return k^3 * ( 
				( ( ( - z2/5760 + 1/144 ) * z3 - z/8 +  1/(2*z) ) + 
				im * ( z2 * ( z2/840 - 1/30 ) + 1/3 ) ) * dotNpNq +
				( ( ( - z2/1152 + 1/48 ) * z3 - z/8 - 1/(2*z) ) + 
				im * z2 * ( z2/210 - 1/15 ) ) * dotrNprNq - 1 / 2 * z2 ) ;
		end
	end

	# Integrando del operador 'Nk' fuera de la diagonal
	# Cuando la norma de z es pequeña se utiliza una aproximación hasta orden 5  
	function Integrando_Nk_robust( k::Real, z::Real, dotrNprNq::Real, dotNpNq::Real )
		if abs( z ) > 1E-2
			return k^3 * exp( im * z ) / z^3 * ( ( 1 - im * z ) * dotNpNq + 
					( 3 * ( im * z - 1 ) + z^2 ) * dotrNprNq ) ;
		else
			z2 = z*z ;
			z3 = z*z2 ;
			return k^3 * ( 
				( ( ( - z2/5760 + 1/144 ) * z3 - z/8 +  1/(2*z) + 1/z3 ) + 
					im * ( z2 * ( z2/840 - 1/30 ) + 1/3 ) ) * dotNpNq +
				( ( ( - z2/1152 + 1/48 ) * z3 - z/8 - 1/(2*z) - 3/z3 ) + 
					im * z2 * ( z2/210 - 1/15 ) ) * dotrNprNq ) ;
		end
	end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 	OPERADOR de Müller
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Función que da el operador discreto 'Muller' para los números de onda 'K0' y 'K1' integrando en el triángulo de
	# vértices 'QA', 'QB', 'QC' y normal 'Nq' con respecto al triángulo de centroide 'P' y normal 'Np'.
	# Si 'P' pertenece al triángulo 'Q' es 'Lponel' = true, de lo contrario es false.
	# La regla de cuadratura es una matriz [ x y w ]: 
	function Operador_Muller( k0::Real, k1::Real, P::Array, Np::Array, QA::Array, QB::Array, 
			QC::Array, Nq::Array, Lponel::Bool, QuadRule::Array )
		# WARNING : La QuadRule cambiará según Lponel. Esto debe asegurarse desde fuera	
		Mk0k1 = ComplexF64( 0 ) ;
		QBMQA = QB - QA ; # 
		QCMQA = QC - QA ; # 
		QCMQB = QC - QB ; # 
		DNPNQ = dot( Np, Nq ) ;
		M = 25 ;
		if Lponel # Hay que regularizar el integrando puesto que P pertenece al triángulo Q 
			# El integrando será : Nk0 - Nk1 - (k0^2 - k1^2) / 2 * L0 , entonces presumo +(k0^2-k1^2)/2*L0 a Mk0k1
			Qarea = area( QA, QB, QC ) ;
			Mk0k1 += ( k0^2 - k1^2 ) / 2 * Operador_L0( P, QA, QB, QC, QBMQA, QCMQA, QCMQB ) / Qarea ;
			for q = 1 : size( QuadRule )[1]
				R = norm( P - Punto_Q( q, QA, QuadRule, QBMQA, QCMQA ) ) ;
				Mk0k1 += QuadRule[ q, 3 ] * ( im / 2 * ( k0^3 - k1^3 ) +
						 SumaTerminosMuller2( R, k0, k1, M ) ) * DNPNQ ;
			end
			return Qarea / ( 4 * pi ) * Mk0k1 ;
		else # "No Lponel" : No hay que regularizar 
			for q = 1 : size( QuadRule )[1]
				Q = Punto_Q( q, QA, QuadRule, QBMQA, QCMQA )
				R = norm( P - Q ) ;
				RNP = dot( P - Q, Np ) / R ;
				RNQ = dot( P - Q, Nq ) / R ; # Aquí es sin el menos
				Mk0k1 += QuadRule[ q, 3 ] * Integrando_Muller_robust( k0, k1, R, RNP * RNQ, DNPNQ ) ;
			end
			return area( QA, QB, QC ) / ( 4 * pi ) * Mk0k1 ;
		end
	end

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Rutinas accesorias para el operador de Müller

	function SumaTerminosMuller1( r::Real, k0::Real, k1::Real, M::Int64 )
		suma = ComplexF64( 0 ) ;
		for n = 3 : M
			suma += - 3 * im^n * ( k1^n - k0^n ) * r^( n - 3 ) / SpecialFunctions.factorial( Float64( n ) ) +
				3 * im^( n + 1 ) * ( k1^( n + 1 ) - k0^( n + 1 ) ) * r^( n - 2 ) / SpecialFunctions.factorial( Float64( n ) ) +
				im^n * ( k1^( n + 2 ) - k0^( n + 2 ) ) * r^( n - 1 ) / SpecialFunctions.factorial( Float64( n ) ) ;
		end
		return suma
	end

	# El término máximo de suma es 27
	function SumaTerminosMuller1( r::Real, k0::Real, k1::Real )
		suma = ComplexF64( 0 ) ;
		Idx = [ 3(i-1) for i=2:15 ] ; # Array de índices para la suma (45 términos máximo)
		for n in Idx
		suma += - 3 * im^n * ( k1^n - k0^n ) * r^( n - 3 ) / SpecialFunctions.factorial( Float64( n ) ) +
			3 * im^( n + 1 ) * ( k1^( n + 1 ) - k0^( n + 1 ) ) * r^( n - 2 ) / SpecialFunctions.factorial( Float64( n ) ) +
			im^n * ( k1^( n + 2 ) - k0^( n + 2 ) ) * r^( n - 1 ) / SpecialFunctions.factorial( Float64( n ) ) ;
			Sum1 = suma ;
		suma += - 3 * im^(n+1) * ( k1^(n+1) - k0^(n+1) ) * r^( n - 2 ) / SpecialFunctions.factorial( Float64( n+1 ) ) +
			3 * im^( n + 2 ) * ( k1^( n + 2 ) - k0^( n + 2 ) ) * r^( n - 1 ) / SpecialFunctions.factorial( Float64( n+1 ) ) +
			im^(n+1) * ( k1^( n + 3 ) - k0^( n + 3 ) ) * r^( n ) / SpecialFunctions.factorial( Float64( n+1 ) ) ;
			Sum2 = suma ;
		suma += - 3 * im^(n+2) * ( k1^(n+2) - k0^(n+2) ) * r^( n - 1 ) / SpecialFunctions.factorial( Float64( n+2 ) ) +
			3 * im^( n + 3 ) * ( k1^( n + 3 ) - k0^( n + 3 ) ) * r^( n ) / SpecialFunctions.factorial( Float64( n+2 ) ) +
			im^(n+2) * ( k1^( n + 4 ) - k0^( n + 4 ) ) * r^( n + 1 ) / SpecialFunctions.factorial( Float64( n+2 ) ) ;
			Sum3 = suma ;
			if Sum3 == Sum2
				if Sum2 == Sum1
					return Sum1
				end
			else
				;
			end
		end
		return suma
	end

	function SumaTerminosMuller2( r::Real, k0::Real, k1::Real, M::Int64 )
		suma = ComplexF64( 0 ) ;
		for n = M : -1 : 3
			suma +=  ( im^(n+1) * ( k1^(n+1) - k0^(n+1) ) * r^(n-2) - im^n * ( k1^n - k0^n ) * r^(n-3) ) / SpecialFunctions.factorial( Float64( n ) ) ;
		end
		return suma
	end

	# El término máximo de suma es 27
	function SumaTerminosMuller2( r::Real, k0::Real, k1::Real )
		suma = ComplexF64( 0 ) ;
		Idx = [ 3(i-1) for i=2:15 ] ; # Array de índices para la suma
		for n in Idx
			suma +=  ( im^(n+1) * ( k1^(n+1) - k0^(n+1) ) * r^(n-2) - im^n * ( k1^n - k0^n ) * r^(n-3) ) / SpecialFunctions.factorial( Float64( n ) ) ;
			Sum1 = suma ;
			suma +=  ( im^(n+2) * ( k1^(n+2) - k0^(n+2) ) * r^(n-1) - im^(n+1) * ( k1^(n+1) - k0^(n+1) ) * r^(n-2) ) / SpecialFunctions.factorial( Float64( n+1 ) ) ;
			Sum2 = suma ;
			suma +=  ( im^(n+3) * ( k1^(n+3) - k0^(n+3) ) * r^(n) - im^(n+2) * ( k1^(n+2) - k0^(n+2) ) * r^(n-1) ) / SpecialFunctions.factorial( Float64( n+2 ) ) ;
			Sum3 = suma ;
			if Sum3 == Sum2
				if Sum2 == Sum1
					return Sum1
				end
			else
				;
			end
		end
		return suma
	end

	function SumaTerminosMuller2_bis( r::Real, k0::Real, k1::Real, z0::Real, z1::Real, M::Int64 )
		suma = ComplexF64( 0 ) ;
		for n = 3 : M
			suma +=  ( im^(n+1) * ( k1^3 * z1^(n-2) - k0^3 * z0^(n-2) ) - im^n * ( k1^3 * z1^(n-3) - k0^3 * z0^(n-3) ) ) / 
					SpecialFunctions.factorial( Float64( n ) ) ;
		end
		return suma
	end

	# Integrando de Müller robustizado para elementos off-diagonal (P,Q pertenecen a triángulos diferentes)
	# No obstante R = norm(P-Q) puede ser pequeño. Es la resta entre los operadores : N_k0 - N_k1 (sin regularizaciones,
	# puesto que nos hallamos fuera de la diagonal)
	function Integrando_Muller_robust( k0::Real, k1::Real, R::Real, RNPRNQ::Real, DNPNQ::Real )
		if min( abs( k0*R ), abs( k1*R ) ) > 1E-2 # 1e-2 EASY
			return ( exp( im*(k0*R) )*( -3 + 3*im*( k0*R ) + ( k0*R )^2 ) - exp( im*(k1*R) )*( -3 + 3*im*(k1*R) + ( k1*R )^2 ) ) / R^3 * RNPRNQ +
				( exp( im*(k1*R) )*( im*(k1*R) - 1 ) - exp( im*(k0*R) )*( im*(k0*R) - 1 ) ) / R^3 * DNPNQ ;
		else # DIFFICULT (versión regularizada)
			return ( 1/(2*R)*( k1^2 - k0^2 ) + im/2*( k1^3 - k0^3 ) + R/2*( k1^4 - k0^4 ) - SumaTerminosMuller1( R, k0, k1 ) )* RNPRNQ +
			( 1/(2*R)*( k0^2 - k1^2 ) + im/2*( k0^3 - k1^3 ) + SumaTerminosMuller2( R, k0, k1 ) ) * DNPNQ ;
		end
	end

	function Integrando_Muller_robust( k0::Real, k1::Real, R::Real, z0::Real, z1::Real, RNPRNQ::Real, DNPNQ::Real )
		if min( abs( z0 ), abs( z1 ) ) > 1E-2 # 1e-2 EASY
			return ( k0^3 * exp( im * z0 ) / z0^3 *( -3 + 3 * im * z0 + z0 * z0 ) - 
				k1^3 * exp( im * z1 ) / z1^3 * ( -3 + 3 * im * z1 + z1 * z1 ) ) * RNPRNQ +
				( k1^3 * exp( im * z1 ) / z1^3 *( im * z1 - 1 ) - k0^3 * exp( im * z0 ) / z0^3 *( im * z0 - 1 ) ) * DNPNQ ;
		else # DIFFICULT (versión regularizada)
			return ( 1/(2*R)*( k1^2 - k0^2 ) + im/2*( k1^3 - k0^3 ) + R/2*( k1^4 - k0^4 ) - SumaTerminosMuller1( R, k0, k1 ) )* RNPRNQ +
			( 1/(2*R)*( k0^2 - k1^2 ) + im/2*( k0^3 - k1^3 ) + SumaTerminosMuller2( R, k0, k1 ) ) * DNPNQ ;
		end
	end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 	CÁLCULO DE TODOS LOS OPERADORES
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	Desde el punto de vista de la conveniencia en tiempos de cómputo conviene realizar la evaluación de los
#	operadores toda de una vez puesto que se ahorran operaciones que de otro modo se hacen duplicadas.
#	Eso no quita la utilidad de tener cada uno de los operadores por separado puesto que en casos particulares
#	(cuando no se utilizan todos a la vez o para ver cosas como Maue) resultan necesarios.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 	Funciones para el caso fluido

	# Función que devuelve todos los operadores evaluados conjuntamente. El punto externo es el centroide 'P'
	# asociado a la normal 'Np' y la integración se lleva a cabo en el triángulo de vértices 'QA, QB, QC' para
	# los números de onda 'k0' (externo) y 'k1' (interno).
	function Operadores_Fluid( k0::Real, k1::Real, P::Array, Np::Array, QA::Array, QB::Array, QC::Array, 
				Nq::Array, Lponel::Bool ; TypeNumber = ComplexF64, RadioCambioQuadRule = 0.25 )
		# WARNING : La QuadRule cambiará según Lponel. Esto debe asegurarse desde fuera	
		Lk0 = zero( TypeNumber ) ;
		Lk1 = zero( TypeNumber ) ;
		Mk0k1 = zero( TypeNumber ) ;
		Mk0 = zero( TypeNumber ) ;
		Mk1 = zero( TypeNumber ) ;
		Mk0t = zero( TypeNumber ) ;
		Mk1t = zero( TypeNumber ) ;        
		QBMQA = QB - QA ; QCMQA = QC - QA ; QCMQB = QC - QB ;
		DNPNQ = dot( Np, Nq ) ;
		if Lponel # Hay que regularizar el integrando puesto que P pertenece al triángulo Q 
			# Regularización de los que necesitan regularizarse
			QuadRule = QR_cqutm11( ) ; # diagonal
			Qarea = area( QA, QB, QC ) ;
			L0 = Operador_L0( P, QA, QB, QC, QBMQA, QCMQA, QCMQB ) / Qarea ;
			Mk0k1 += ( k0^2 - k1^2 ) / 2 * L0 ;
			Lk0 += L0 ;
			Lk1 += L0 ;
			for q = 1 : size( QuadRule )[ 1 ]
				R = norm( P - Punto_Q( q, QA, QuadRule, QBMQA, QCMQA ) ) ;
				Mk0k1 += QuadRule[ q, 3 ] * ( im / 2 * ( k0^3 - k1^3 ) + 
						SumaTerminosMuller2( R, k0, k1 ) ) * DNPNQ ;
				Lk0 += QuadRule[ q, 3 ] * Integrando_Lk_regular( k0, R ) ;
				Lk1 += QuadRule[ q, 3 ] * Integrando_Lk_regular( k1, R ) ;
			end
			return Qarea / ( 4 * pi ) * [ Mk0k1, Lk0, Lk1, Mk0, Mk1, Mk0t, Mk1t ] ;
		else # "No Lponel" : No hay que regularizar 
			QuadRule = QuadRule_selector( P, QA, QB, QC, RadioCambioQuadRule ) ; # Off diagonal
			for q = 1 : size( QuadRule )[ 1 ]
				Q = Punto_Q( q, QA, QuadRule, QBMQA, QCMQA ) ;
				PmQ = P - Q ;
				R = norm( PmQ ) ;
				RdotNq = dot( PmQ, Nq ) ;
				RdotNp = dot( PmQ, Np ) ;
				Mk0k1 += QuadRule[ q, 3 ] * Integrando_Muller_robust( k0, k1, R, (RdotNp / R) * 
							(RdotNq / R), DNPNQ ) ;
				Lk0 += QuadRule[ q, 3 ] * Integrando_Lk_robust( k0, R ) ;
				Lk1 += QuadRule[ q, 3 ] * Integrando_Lk_robust( k1, R ) ;
				Mk0 += QuadRule[ q, 3 ] * Integrando_Mk_robust( k0, R ) * -RdotNq ;
				Mk1 += QuadRule[ q, 3 ] * Integrando_Mk_robust( k1, R ) * -RdotNq ;
				Mk0t += QuadRule[ q, 3 ] * Integrando_Mk_robust( k0, R ) * RdotNp ;
				Mk1t += QuadRule[ q, 3 ] * Integrando_Mk_robust( k1, R ) * RdotNp ;
			end
			return area( QA, QB, QC ) / ( 4 * pi ) * [ Mk0k1, Lk0, Lk1, Mk0, Mk1, Mk0t, Mk1t ] ;
		end
	end
    
	# Función que devuelve una regla de cuadratura para integrar un triángulo de centroide 'C' con
	# respecto a un centroide 'P' dependiendo si la distancia entre centroides es mayor o menor que
	# un cierto valor umbral. Esto posibilita utilizar reglas sencillas cuando los elementos están
	# muy alejados entre sí y todos las distancias |P - Q| entre el centroide los puntos Q de la
	# cuadratura son parecidas.
	function QuadRule_selector( P, QA, QB, QC, Radio )
		# Centroide del triángulo 'Q' es C = (QA+QB+QC)/3
		if norm( P - ( QA + QB + QC ) / 3 ) > Radio
			return QR_cqutm6() ; # Regla de cuadratura de menor precisión
		else
			return QR_cqutm9() ; # Regla de cuadratura de mayor precisión
		end 
	end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 	Funciones para el caso de shells

	# Función que llena los bloques off diagonal en el sistema resultante para un shell
	function Operadores_Shell_offdiagonal( K::Real, P::Array, Np::Array, QA::Array, QB::Array, QC::Array, 
				Nq::Array ; TypeNumber = ComplexF64, RadioCambioQuadRule = 0.25 )
		# WARNING : La QuadRule cambiará según Lponel. Esto debe asegurarse desde fuera	
		Nk = zero( TypeNumber ) ;
		Lk = zero( TypeNumber ) ;
		Mk = zero( TypeNumber ) ;
		Mkt = zero( TypeNumber ) ;
		QBMQA = QB - QA ; QCMQA = QC - QA ; QCMQB = QC - QB ;
		DNPNQ = dot( Np, Nq ) ;
		QuadRule = QuadRule_selector( P, QA, QB, QC, RadioCambioQuadRule ) ; # Off diagonal
		for q = 1 : size( QuadRule )[ 1 ]
				Q = Punto_Q( q, QA, QuadRule, QBMQA, QCMQA ) ;
				PmQ = P - Q ;
				R = norm( PmQ ) ;
		RdotNq = dot( PmQ, Nq ) ;
		RdotNp = dot( PmQ, Np ) ;
				FPGR = exp( im * K * R ) / R^2 * ( im * K * R - 1 ) ; 
				FPGRR = exp( im * K * R ) / R^3 * ( 2 - 2 * im * K * R - K^2 * R^2 ) ;
				RNPNQ = -( DNPNQ - RdotNp/R * RdotNq/R ) / R ;
				Nk += QuadRule[ q, 3 ] * ( FPGR * RNPNQ + FPGRR * RdotNp/R * -RdotNq/R ) ;
		Lk += QuadRule[ q, 3 ] * Integrando_Lk_robust( K, R ) ;
		Mk += QuadRule[ q, 3 ] * Integrando_Mk_robust( K, R ) * -RdotNq ;
		Mkt += QuadRule[ q, 3 ] * Integrando_Mk_robust( K, R ) * RdotNp ;
			end
		return area( QA, QB, QC ) / ( 4 * pi ) * [ Nk, Lk, Mk, Mkt ] ;
	end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 	Funciones para el caso impenetrable

	# Función que calcula los operadores necesarios para el caso de Dirichlet BC
	# Reemplaza a h3lc
	function Operadores_Dir( K::Real, P::Array, Np::Array, Lponel::Bool, QuadRule::Array, 
			QA::Array, QB::Array, QC::Array, Nq::Array ; TypeNumber = ComplexF64 )
		# WARNING : La QuadRule cambiará según Lponel. Esto debe asegurarse desde fuera	
		Lk = zero( TypeNumber ) ;
		Mkt = zero( TypeNumber ) ;   
		QBMQA = QB - QA ; QCMQA = QC - QA ; QCMQB = QC - QB ;
		# DNPNQ = dot( Np, Nq ) ;
		if Lponel # Hay que regularizar el integrando (la regla debería ser una "QuadRuleOn")
			Qarea = area( QA, QB, QC ) ;
			Lk += Operador_L0( P, QA, QB, QC, QBMQA, QCMQA, QCMQB ) / Qarea ;
			for q = 1 : size( QuadRule )[ 1 ]
				R = norm( P - Punto_Q( q, QA, QuadRule, QBMQA, QCMQA ) ) ;
				Lk += QuadRule[ q, 3 ] * Integrando_Lk_regular( K, R ) ;
			end
			return Qarea / ( 4 * pi ) * [ Lk, Mkt ] ;
		else # "No Lponel" : No hay que regularizar 
			for q = 1 : size( QuadRule )[ 1 ]
				PmQ = P - Punto_Q( q, QA, QuadRule, QBMQA, QCMQA ) ;
				R = norm( PmQ ) ;
				RdotNp = dot( PmQ, Np ) ;
				Lk += QuadRule[ q, 3 ] * Integrando_Lk_robust( K, R ) ;
				Mkt += QuadRule[ q, 3 ] * Integrando_Mk_robust( K, R ) * RdotNp ;
			end
			return area( QA, QB, QC ) / ( 4 * pi ) * [ Lk, Mkt ] ;
		end
	end

	# Función que calcula los operadores necesarios para el caso de Neumann BC
	# Reemplaza a h3lc
	function Operadores_Neu( K::Real, P::Array, Np::Array, Lponel::Bool, QuadRule::Array, 
		    QA::Array, QB::Array, QC::Array, Nq::Array ; TypeNumber = ComplexF64 )
		# WARNING : La QuadRule cambiará según Lponel. Esto debe asegurarse desde fuera	
		Mk = zero( TypeNumber ) ;
		Nk = zero( TypeNumber ) ;
		DNPNQ = dot( Np, Nq ) ;
		if Lponel # Hay que regularizar el integrando (la regla debería ser una "QuadRuleOn")
			Qarea = area( QA, QB, QC ) ;
			QBMQA = QB - QA ; QCMQA = QC - QA ; QCMQB = QC - QB ;
			Nk += Operador_N0_L0( K, P, QA, QB, QC, QBMQA, QCMQA, QCMQB ) / Qarea ;
			for q = 1 : size( QuadRule )[ 1 ]
				PmQ = P - Punto_Q( q, QA, QuadRule, QBMQA, QCMQA ) ;
				R = norm( PmQ ) ;
				r = PmQ / R ;
				z = K * R ;
				dotrNprNq = dot(r,Np) * dot(r,Nq) ;
				Nk += QuadRule[ q, 3 ] * Integrando_Nk_regular( K, z, dotrNprNq, DNPNQ ) ;
			end
			return Qarea / ( 4 * pi ) * [ Mk, Nk ] ;
		else # "No Lponel" : No hay que regularizar 
			QBMQA = QB - QA ; QCMQA = QC - QA ; 
			for q = 1 : size( QuadRule )[1]
				PmQ = P - Punto_Q( q, QA, QuadRule, QBMQA, QCMQA ) ;
				R = norm( PmQ ) ;
				r = PmQ / R ;
				z = K * R ;
				dotrNprNq = dot( r, Np ) * dot( r, Nq ) ;
				Mk += QuadRule[ q, 3 ] * Integrando_Mk_robust( K, R ) * -dot( PmQ, Nq )  ; 
				Nk += QuadRule[ q, 3 ] * Integrando_Nk_robust( K, z, dotrNprNq, DNPNQ ) ;
			end
			return area( QA, QB, QC ) / ( 4 * pi ) * [ Mk, Nk ] ;
		end
	end

