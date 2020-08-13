# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	
#	Julia parallel routines
#	
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# v. 13/8/20

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Index:

# 	SendToProc( p::Int ; args... )
# 	GetFromProc( p::Int, nm::Symbol ; mod=Main )
	
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	# Defines locally a relation 'nm = val' given in 'args' and
	# send to process 'p' the variable 'nm'
	function SendToProc( p::Int ; args... )
	    for ( nm, val ) in args
		@spawnat( p, Core.eval( Main, Expr( :(=), nm, val) ) ) ;
	    end
	end

	# Brings back from process 'p' a variable 'nm' (it is entered as ':nm')
	function GetFromProc( p::Int, nm::Symbol ; mod=Main )
		fetch( @spawnat( p, getfield( mod, nm ) ) ) ;
	end
