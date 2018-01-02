#ifndef NODE_H__
#define NODE_H__

using namespace std;


class NODE{
	int		m;
	double	*ub;
	double	*lb;
	double	*warm;
	double	*sumfixed;
	double	relax_objval;
	double	*relax_solval;
	int		dpt;
	bool		zero;
	int		index;
	bool		solved;
	
	public:

		NODE();											// default constructor	
		NODE( const NODE &source );			// copy constructor
		NODE& operator=( const NODE& );		// assignment operator
		~NODE();											// destructor
		bool operator<(const NODE &rhs) const
		{ return relax_objval < rhs.relax_objval; }

		void	set_vals( int s_m, double *s_ub, double *s_lb,
								double *s_warm, double s_relax_objval,
								int s_dpt, bool s_zero, int s_index);

		double	get_lowerbound(){ return relax_objval; }
		int		get_index(){ return index; }
		int		get_dpt(){ return dpt; }
		bool		get_zero(){ return zero; }
		double*	get_ub(){ return ub; }
		double*	get_lb(){ return lb; }
		double*	get_warm(){ return warm; }
		double*	get_relaxsolval(){ return relax_solval; }
		bool		get_solved(){ return solved; }

		void	set_lowerbound( double s_ro ){
			if( relax_objval < s_ro ) relax_objval = s_ro;
		}
		void	set_relaxsolval( double *solval );
		void	set_solved( bool s_solved ){ solved = s_solved; }
		bool	alloc_sumfixed();
		void	set_sumfixed( double c, double *s_sumfixed );
		void	add_sumfixed( double c, double *s_sumfixed );
		double*	get_sumfixed(){ return sumfixed; }
};

#endif
