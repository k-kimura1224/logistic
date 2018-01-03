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
   b = - grad * p_max + func;
   output.set_pab( p_max, a, b );

   assert( p_max > p1 );
   assert( p_max < p2 );

   return output;
}

/** exit a program if a input file has error */
static
bool readFail(
   int                   status
   )
{
   if( !status )
   {
      printf("Reading failed.\n");
      return false;
   }

   return true;
}


/** read dimension from the input file */
static
bool readDataDim(
   const char            *filename,          /**< name of the input file */
   int*                  n,                  /**< pointer to store the number of data points*/
   int*                  p,                  /**< pointer to store the  number of explanatory variables */
   int*                  i_ex                /**< pointer to store the index of explained varaible */
   )
{
   FILE *file;
   int status;
   char s[1000];

   assert(filename != NULL);
   assert(n != NULL);
   assert(p != NULL);
   assert(i_ex != NULL);

   /* open file */
   file = fopen(filename, "r");
   if( file == NULL )
   {
      printf("Could not open file <%s>.\n", filename);
      return false;
   }

   /* skip one line */
   if( fgets(s, 1000, file) == NULL )
   {
      printf("Error reading file <%s>.\n", filename);
      return false;
   }

   /* read n */
   status = fscanf(file, "%d", n);
   readFail(status);

   /* read p */
   status = fscanf(file, "%d", p);
   readFail(status);

   /* read i_ex */
   status = fscanf(file, "%d", i_ex);
   readFail(status);

   /* close file */
   fclose(file);

   return true;
}


/** read data points */
static
bool readData(
   const char*           filename,           /**< name of the input file */
   int                   n,                  /**< the number of data points */
   int                   p,                  /**< the number of explanatory varaibles */
   double*               data                /**< array to store data points */
   )
{
   int i;
   FILE *file;
   int buf;
   int status;
   char s[1000];

   /* open file */
   file = fopen(filename, "r");
   if( file==NULL )
   {
      printf("Could not open file <%s>.\n", filename);
      return false;
   }

   /* skip one line */
   if( fgets(s, 1000, file) == NULL )
   {
      printf("Error reading file <%s>.\n", filename);
      return false;
   }

   /* skip n */
   status = fscanf(file, "%d", &buf);
   readFail( status );

   /* skip p */
   status = fscanf(file, "%d", &buf);
   readFail( status );

   /* skip i_ex */
   status = fscanf(file, "%d", &buf);
   readFail( status );

   /* read data points */
   for( i = 0; i < n * ( p + 1 ); i++ )
   {
      status = fscanf(file, "%lf", (data + i));
      readFail( status );
   }

   /* close file */
   fclose(file);

   return true;
}


/** calculate the mean values for normalization */
static
bool calcMean(
   int                   n,                  /**< the number of data points */
   int                   p,                  /**< the number of explanatory varaibles */
   double*               explanatory,               /**< data points */
   double*               mean                /**< array to store the mean values */
   )
{
   int i;
   int j;

   assert(n > 0);
   assert(p > 0);
   assert(explanatory != NULL);
   assert(mean != NULL);

   for( i = 0; i < p; i++ )
   {
      mean[i] = 0.0;

      for( j = 0; j < n; j++ )
         mean[i] += explanatory[(j*p)+i];

      mean[i] /= (double) n;
   }

   return true;
}


/** calculate the variance values for normalization */
static
bool calcVariance(
   int                   n,                  /**< the number of data points */
   int                   p,                  /**< the number of explanatory varaibles */
   double*            explanatory,               /**< data points */
   double*            mean,               /**< the mean values */
   double*            variance            /**< array to store the variance values */
   )
{
   int i;
   int j;

   assert(n > 1);
   assert(p > 0);
   assert(explanatory != NULL);
   assert(mean != NULL);
   assert(variance != NULL);

   for( i = 0; i < p; i++ )
   {
      variance[i] = 0.0;

      for( j = 0; j < n; j++ )
         variance[i] += pow( explanatory[j*p + i] - mean[i], 2.0);

      variance[i] *= 1.0 / ((double) n - 1.0);
   }

   return true;
}



/** normalize data */
static
bool normalization(
   int                   n,                  /**< the number of data points */
   int                   p,                  /**< the number of explanatory variables */
   double*               explanatory         /**< data points */
   )
{
   int i;
   int j;
   double* mean;
   double* variance;

   assert(n > 0);
   assert(p > 0);
   assert(explanatory != NULL);

   /* allocate memory for mean and variance */
   mean = new double[p];
   variance = new double[p];

   /* calculate the mean values and the variance values */
   calcMean(n, p, explanatory, mean );
   calcVariance(n, p, explanatory, mean, variance);

   for( i = 0; i < n; i++ )
   {
      for( j = 0; j < p; j++ )
      {
         assert( !Equal( variance[j], 0.0 ) );
         explanatory[i*p + j] = (explanatory[i*p + j] - mean[j] ) / sqrt(variance[j]);
      }
   }

   /* free */
   delete[] mean;
   delete[] variance;

   return true;
}


/** divide data into explained variable and explanatory variables */
static
bool divideData(
   int                   n,                  /**< the number of data points */
   int                   p,                  /**< the number of explanatory variables */
   int                   i_ex,               /**< index of the explained variable */
   double*               data,               /**< data points */
   bool*                 explained,          /**< array to store the explained variable */
   double*               explanatory         /**< array to store the explanatory variables */
   )
{
   int i;
   int j;
   int ct = 0;

   assert(n > 0);
   assert(p > 0);
   assert(i_ex > 0 && i_ex <= p+1);
   assert(data != NULL);
   assert(explained != NULL);
   assert(explanatory != NULL);

   for( i = 0; i < n; i++ )
   {
      *(explained + i) = (bool) *(data + (i * (p + 1)) + (i_ex - 1));
   }

   for( i = 0; i < n; i++ )
   {
      for( j = 0; j < p + 1; j++ )
      {
         if( j != i_ex - 1 )
         {
            *(explanatory + ct) = *(data + (i * (p + 1)) + j);
            ct++;
         }
      }
   }

   return true;
}

static string sign( double x )
{
   if ( x >= 0.0 )
      return " + ";
   else
      return " - ";
}

static
void objectfunction(
   int      n,
   int      p,
   bool*    explained,
   double*  explanatory
   )
{
   cout << "minimize" << endl;

   cout << " 2 t_1";
   for ( auto i = 2; i <= n; i++ )
      cout << " + 2 t_" << i;

   for ( auto i = 0; i < p+1; i++ )
      cout << " + 2 z_" << i;

   cout << endl;
}

static
void ct_const(
   int*   ct
   )
{
   printf("c%d:  ", *ct);
   *ct += 1;
}

static
void constraints(
   int      n,
   int      p,
   bool*    y,
   double*  X,
   POINT*   v,
   int      k
   )
{
   auto ct = 1;

   printf("subject to\n");

   POINT point;
   double a;
   double b;
   double c;
   double d;

   for ( auto s = 0; s < k; s++ )
   {
      point = v[s];
      a = point.get_a();
      b = point.get_b();

      for ( auto i = 0; i < n; i++ )
      {
         ct_const( &ct );
         cout << "t_" << i;

         if( y[i] == true )
            d = - a;
         else
            d = a;

         cout << sign( d ) << fabs( d ) <<" b_0";
         for( auto j = 0; j < p; j++ )
         {
            c = d * X[i*p + j];
            cout << sign( c ) << fabs( c ) << " b_" << j+1;
         }
        cout << " >= " << b << endl;
      }
   }

   // indicator constraints
   for ( auto j = 0; j < p+1; j++ )
   {
      ct_const( &ct );
      cout << "z_" << j << " = 0 -> " << "b_" << j << " = 0" << endl;
   }
}

static
void boundary(
      int   n,
      int   p
      )
{
   cout << "bounds" << endl;

   // free
   for ( auto j = 0; j < p+1; j++ )
      cout << "b_" << j << " free" << endl;

   for ( auto i = 0; i < n; i++ )
      cout << "t_" << i << " free" << endl;

   // binary
   cout << "binary" << endl;
   for ( auto j = 0; j < p+1; j++ )
      cout << "z_" << j << endl;
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


   //cout << "points_5 = { ";
   //for( auto p : points_5 )
   //   cout << p.get_p() << ", ";
   //cout << "}" << " --> " << points_5.size() << endl;

   //cout << "points_9 = { ";
   //for( auto p : points_9 )
   //   cout << p.get_p() << ", ";
   //cout << "}" << " --> " << points_9.size() << endl;

   //cout << "points_17 = { ";
   //for( auto p : points_17 )
   //   cout << p.get_p() << ", ";
   //cout << "}" << " --> " << points_17.size() << endl;



   int   n, p, i_ex;
   double*  data;
   bool*  explained;
   double*  explanatory;

   if( argc != 3 )
   {
      cout << "Error : commandline arguments" << endl;
      return -1;
   }

   int V = atoi( argv[2] );
   assert( V == 5 || V == 9 || V == 17 );

   auto result = readDataDim( argv[1], &n, &p, &i_ex );

   if ( !result )
   {
      cout << "Error: read dim of data" << endl;
      return -1;
   }

   //cout << endl;
   //cout << "n = " << n << endl;
   //cout << "p = " << p << endl;
   //cout << "i_ex = " << i_ex << endl;
   //cout << endl;

   data = new double[ n * ( p + 1 ) ];
   result = readData( argv[1], n, p, data );

   if ( !result )
   {
      cout << "Error: read data" << endl;
      return -1;
   }

   //cout << "reading .." << endl;
   //cout << data[0] << ", " << data[1] << " .. " << data[n*p + n - 2];
   //cout << ", " << data[n*p + n - 1] << endl;

   explained = new bool[n];
   explanatory = new double[n*p];

   divideData( n, p, i_ex, data, explained, explanatory );
   normalization( n, p, explanatory );

   //for ( auto i = 0; i < p; i++ )
   //{
   //   auto sum = 0.0;
   //   auto ct = i;
   //   for ( auto j = 0; j < n; j++ )
   //   {
   //      cout << explanatory[ct] << " ";
   //      sum += explanatory[ct];
   //      ct += p;
   //   }
   //   cout << "\n " << sum << endl;
   //}
   //auto ct = 0;
   //for ( auto i = 0; i < n; i++ )
   //{
   //   if( explained[i] == true ) ct++;
   //}
   //cout << ct << endl;

   objectfunction( n, p, explained, explanatory );
   if ( V == 5 )
   {
      assert( V == (int)points_5.size() );
      constraints( n, p, explained, explanatory, &points_5[0], V );
   }
   else if ( V == 9 )
   {
      assert( V == (int)points_9.size() );
      constraints( n, p, explained, explanatory, &points_9[0], V );
   }
   else
   {
      assert( V == (int)points_17.size() );
      constraints( n, p, explained, explanatory, &points_17[0], V );
   }

   boundary( n, p );

   cout << "end" << endl;

   delete[] data;
   delete[] explained;
   delete[] explanatory;

   return 0;
}
