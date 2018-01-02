#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <time.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <iomanip>

#include <vector>
#include <algorithm>

#define ep 1.0e-10
using namespace std;

// return true if given values are equal
template <typename T>
bool Equal(
	   T	a,
		T  b
	)
{
	bool rslt = false;
	auto max = fabs( a );

	if( max < fabs( b ) )	max = fabs( b );
	if( max < 1.0 ) 	max = 1.0;

	if( fabs( a - b ) < max * ep ) rslt = true;

	return rslt;
}

class POINT{
   double point;
   double a; // y = a x + b
   double b;

   public:
      POINT();
		POINT( const POINT &source );			// copy constructor
		POINT& operator=( const POINT& );		// assignment operator
		~POINT();											// destructor
		bool operator<(const POINT &rhs) const
		{ return point < rhs.point; }

      void set_pab( double s_p, double s_a, double s_b)
      {
         point = s_p;
         a = s_a;
         b = s_b;
      }
      auto get_p() { return point; }
      auto get_a() { return a; }
      auto get_b() { return b; }

};

POINT::POINT()
{ // default constructor
   point = 0.0;
   a = 0.0;
   b = 0.0;
}

POINT::POINT( const POINT &source )
{	// copy constructor
   point = source.point;
   a = source.a;
   b = source.b;
}

POINT& POINT::operator=( const POINT& source )
{ //assignment operator
   if( this != &source )
   {
      point = source.point;
      a = source.a;
      b = source.b;
   }

	return *this;
}

POINT::~POINT()
{ // destructor

}

double FindIntersectionPoint(
      double   a1,
      double   b1,
      double   a2,
      double   b2
      )
{
   assert( !Equal( a1 - a2, 0.0 ) );
   return ( b2 - b1 ) / ( a1 - a2 );
}

POINT GeneratePoint(
      POINT    input1,
      POINT    input2
      )
{
   double  x1, y1; // an intersection point

   auto  p1 = input1.get_p();
   auto  a1 = input1.get_a();
   auto  b1 = input1.get_b();

   auto  p2 = input2.get_p();
   auto  a2 = input2.get_a();
   auto  b2 = input2.get_b();

   // find an intersection point
   x1 = FindIntersectionPoint( a1, b1, a2, b2 );
   y1 = a1 * x1 + b1;

   //cout << "( x1, y1 ) = ( " << x1 << ", " << y1 << " )" << endl;

   auto  lb = p1;
   auto  ub = p2;

   assert( lb <= ub );
   assert( ub - lb > 1.0e-5 );

   double   gap = ub - lb;
   auto     k = 3;
   double   t;
   double   testpoint;
   double   grad;
   double   func;
   double   buf;
   double   a, b;
   double   x2, y2;
   double   x3, y3;
   double   S;
   double   S_max;
   double   p_max;

   do
   {
      t = gap / (double) (k + 1);

      testpoint = lb;
      S_max = 0.0;
      p_max = lb;
      for( auto i=0; i<k; i++ )
      {
         testpoint += t;
         buf = 1 + exp( -testpoint );
         func = log( buf );
         grad = 1.0/buf - 1.0;
         a = grad;
         b = - grad * testpoint + func;

         //cout << "testpoint = " << testpoint << endl;
         //cout << "buf = " << buf << endl;
         //cout << "func = " << func << endl;
         //cout << "grad = " << grad << endl;
         //cout << "a = " << a << endl;
         //cout << "b = " << b << endl;

         x2 = FindIntersectionPoint( a1, b1, a, b );
         y2 = a1 * x2 + b1;

         x3 = FindIntersectionPoint( a2, b2, a, b );
         y3 = a2 * x3 + b2;

         S = 0.5 * fabs( ( x1 - x3 ) * ( y2 - y3 ) - ( x2 - x3 ) * ( y1 - y3 ) );

         //cout << "S = " << S << endl;

         if( S > S_max )
         {
            S_max = S;
            p_max = testpoint;
         }
         else
         {
            break;
         }
      }

      ub = p_max + t;
      lb = p_max - t;
      gap = ub - lb;

      //cout << "(ub, lb) = (" << ub << ", " << lb << ")" << endl;
      //cout << "gap = " << gap << endl;
      //cout << "----------------------------" << endl;

   } while( gap > 1.0e-6 );

   POINT output;

   buf = 1 + exp( - p_max );
   func = log( buf );
   grad = 1.0/buf - 1.0;
   a = grad;
   b = - grad * testpoint + func;
   output.set_pab( p_max, a, b );

   assert( p_max > p1 );
   assert( p_max < p2 );

   return output;
}


int main( int argc, char** argv )
{
   auto infinity = 1.0e+1;
   vector<POINT> points_5;
   vector<POINT> points_9;
   vector<POINT> points_17;

   // initialize points_5
   POINT inf;
   POINT minf;

   minf.set_pab( -infinity, -1.0, 0.0 );
   inf.set_pab( infinity, 0.0, 0.0 );

   POINT ps3[3];
   ps3[0] = GeneratePoint( minf, inf );
   ps3[1] = GeneratePoint( minf, ps3[0] );
   ps3[2] = GeneratePoint( ps3[0], inf );

   points_5.push_back( minf );
   points_5.push_back( ps3[1] );
   points_5.push_back( ps3[0] );
   points_5.push_back( ps3[2] );
   points_5.push_back( inf );

   // initialize points_9
   points_9 = points_5;

   POINT ps4[4];
   for( auto i = 0; i < 4; i++ )
      ps4[i] = GeneratePoint( points_5[i], points_5[i+1] );

   for( auto i = 0; i < 4; i++ )
      points_9.push_back( ps4[i] );

   sort( points_9.begin(), points_9.end() );

   // initialize points_17
   points_17 = points_9;

   POINT ps8[8];
   for( auto i = 0; i < 8; i++ )
      ps8[i] = GeneratePoint( points_9[i], points_9[i+1] );

   for( auto i = 0; i < 8; i++ )
      points_17.push_back( ps8[i] );

   sort( points_17.begin(), points_17.end() );


   cout << "points_5 = { ";
   for( auto p : points_5 )
      cout << p.get_p() << ", ";
   cout << "}" << " --> " << points_5.size() << endl;

   cout << "points_9 = { ";
   for( auto p : points_9 )
      cout << p.get_p() << ", ";
   cout << "}" << " --> " << points_9.size() << endl;

   cout << "points_17 = { ";
   for( auto p : points_17 )
      cout << p.get_p() << ", ";
   cout << "}" << " --> " << points_17.size() << endl;

   return 0;
}
