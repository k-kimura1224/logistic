#ifndef CUT_PLANE_H__
#define CUT_PLANE_H__

using namespace std;


class CUT_PLANE{
	int		m;
	double	*coef;	// [m]
	double	*lb;		// [1]
	double	*ub;		// [1]
	
	public:

		CUT_PLANE();											// default constructor	
		CUT_PLANE( const CUT_PLANE &source );			// copy constructor
		CUT_PLANE& operator=( const CUT_PLANE& );		// assignment operator
		~CUT_PLANE();											// destructor

		void set_cut(int s_m, double *s_coef, double s_lb, double s_ub);

		double*	get_coef(){ return coef; }
		double*	get_lb(){ return lb; }
		double*	get_ub(){ return ub; }
};

#endif
