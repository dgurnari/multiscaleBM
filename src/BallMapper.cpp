// #define rcpp_code

#ifdef rcpp_code
  #include <Rcpp.h>
  using namespace Rcpp;
#endif
#include <vector>
#include <map>
#include <iostream>
#include <cmath>
#include <tuple>
#include <algorithm>

// using namespace std;


// For standard non symmetric data
template <typename point>
inline double compute_distance
  ( const std::vector< point >& pts , size_t f , size_t s ,  double p = 2 )
{
  double result = 0;
  bool dbg = false;

  for ( size_t i = 0 ; i != pts[f].size() ; ++i )
  {
     result += pow( ( pts[f][i]-pts[s][i] ) , p );
  }

  return pow(result,1/p);
}//compute_distance



template <typename point>
void compute_landmarks
(
 const std::vector< point >& points ,
 std::vector< std::vector<size_t> >& coverage,
 std::vector< size_t >& landmarks ,
 double epsilon , int number_of_points
)
{
  bool dbg = false;
  int first_not_covered_point = 0;
  int number_of_landmark = -1;
  while ( first_not_covered_point < number_of_points )
  {
     if (dbg)
     {
			std::cout << "first_not_covered_point : " << first_not_covered_point << std::endl;
     }

     landmarks.push_back( first_not_covered_point+1 );//!!!
     ++number_of_landmark;
     //check what is covered by points_vec[first_not_covered_point]:
     for ( int i = 0 ; i != number_of_points ; ++i )
     {
       if ( compute_distance<point>( points , first_not_covered_point , i ) <= epsilon )
       {
         coverage[i].push_back( number_of_landmark+1 );
       }
     }
     //now find the next not covered point:
     while ( (first_not_covered_point!=number_of_points) && (coverage[first_not_covered_point].size() != 0) )
     {
       ++first_not_covered_point;
     }
  }
} //compute_landmarks



//internal procedure
template < typename vect >
void internal_procedure_fill_points_covered_by_landmarks
(
vect& numer_of_covered_points ,
std::vector< std::vector< int > >& points_covered_by_landmarks ,
const std::vector< size_t >& landmarks ,
const std::vector< std::vector<size_t> >& coverage
)
{
	bool dbg = true;
	//Now we will compute points_covered_by_landmarks. Firstly let us initialize all the structures:

	for ( size_t i = 0 ; i != coverage.size() ; ++i )
	{
		for ( size_t j = 0 ; j != coverage[i].size() ; ++j )
		{
			numer_of_covered_points[ coverage[i][j]-1 ]++;
		}
	}

	points_covered_by_landmarks = std::vector< std::vector< int > > ( landmarks.size() );
	for ( size_t i = 0 ; i != landmarks.size() ; ++i )
	{
		std::vector<int> aa;
		aa.reserve( numer_of_covered_points[i] );
		points_covered_by_landmarks[i] = aa;
	}
	//now when the variables are initialized, we can fill in the numer_of_covered_points list:
	for ( size_t i = 0 ; i != coverage.size() ; ++i )
	{
		for ( size_t j = 0 ; j != coverage[i].size() ; ++j )
		{
			points_covered_by_landmarks[ coverage[i][j]-1 ].push_back( i+1 );
		}
	}

	if (dbg)
	{
		#ifdef rcpp_code
			Rstd::cerr << "points_covered_by_landmarks.size() : " << points_covered_by_landmarks.size() << std::endl;
		#endif
		#ifndef rcpp_code
			std::cerr << "points_covered_by_landmarks.size() : " << points_covered_by_landmarks.size() << std::endl;
		#endif
	}
}//internal_procedure_fill_points_covered_by_landmarks



template <typename vect>
void internal_procedure_fill_coloring
(
std::vector< double >& coloring ,
const std::vector< std::vector< int > >& points_covered_by_landmarks ,
const vect& values
)
{
	double dbg = false;
	coloring = std::vector< double >( points_covered_by_landmarks.size() , 0 );
	for ( size_t i = 0 ; i != points_covered_by_landmarks.size() ; ++i )
	{
		double av = 0;
		for ( size_t j = 0 ; j != points_covered_by_landmarks[i].size() ; ++j )
		{
			av += values[ points_covered_by_landmarks[i][j]-1 ];
		}
		av = av / (double)points_covered_by_landmarks[i].size();
		coloring[i] = av;
	}

	if (dbg)
	{
		#ifdef rcpp_code
			Rstd::cerr << "Here is the coloring : \n";
			for ( size_t i = 0 ; i != coloring.size() ; ++i )
			{
				std::cerr << coloring[i] << " , ";
			}
		#endif
		#ifndef rcpp_code
		    std::cerr << "Here is the coloring : \n";
			for ( size_t i = 0 ; i != coloring.size() ; ++i )
			{
				std::cerr << coloring[i] << " , ";
			}
		#endif
	}
}//internal_procedure_fill_coloring



template < typename vect >
void
internal_procedure_build_graph
(
std::vector< std::vector< int > >& graph_incidence_matrix,
vect& from,
vect& to,
vect& strength_of_edges,
const std::vector< size_t >& landmarks,
const std::vector< std::vector<size_t> > coverage
)
{
  bool dbg = false;
  graph_incidence_matrix = std::vector< std::vector< int > >( landmarks.size() );
  for ( size_t i = 0 ; i != landmarks.size() ; ++i )
  {
    graph_incidence_matrix[i] = std::vector<int>( i );
  }

  if (dbg)
  {
	  #ifdef rcpp_code
			Rstd::cerr << "coverage.size() : " << coverage.size() << std::endl;
	  #endif
	  #ifndef rcpp_code
			std::cerr << "coverage.size() : " << coverage.size() << std::endl;
	  #endif
   }

  for ( size_t i = 0 ; i < coverage.size() ; ++i )
  {
    for ( size_t j = 0 ; j < coverage[i].size() ; ++j )
    {
      for ( size_t k = j+1 ; k < coverage[i].size() ; ++k )
      {
        //note that the landmarks in coverage are sorted from smallest to largest
        //therefore we only need this:
        //Rstd::cerr << coverage[i][k] << " " << coverage[i][j] << std::endl;
        graph_incidence_matrix[ coverage[i][k]-1 ][ coverage[i][j]-1  ]++;
      }
    }
  }


  //first let us count the number of edges in the graph:
  int number_of_edges = 0;
  for ( size_t i = 0 ; i != graph_incidence_matrix.size() ; ++i )
  {
    for ( size_t j = 0 ; j != graph_incidence_matrix[i].size() ; ++j )
    {
       if ( graph_incidence_matrix[i][j] != 0 )++number_of_edges;
    }
  }

  if (dbg)
  {
	  #ifdef rcpp_code
		Rstd::cerr << "Number of edges in the graph : " << number_of_edges << std::endl;
	  #endif
	  #ifndef rcpp_code
		std::cerr << "Number of edges in the graph : " << number_of_edges << std::endl;
	  #endif
  }


  from = vect(number_of_edges);
  to = vect(number_of_edges);
  strength_of_edges = vect(number_of_edges);

  int edg_no = 0;

  for ( size_t i = 0 ; i != graph_incidence_matrix.size() ; ++i )
  {
    for ( size_t j = 0 ; j != graph_incidence_matrix[i].size() ; ++j )
    {
      if ( graph_incidence_matrix[i][j] != 0 )
      {
         from[edg_no] = j+1;
         to[edg_no] = i+1;
         strength_of_edges[edg_no] = graph_incidence_matrix[i][j];
         ++edg_no;
      }
    }
  }
}//internal_procedure_build_graph



std::tuple< std::vector< int > , std::vector<std::pair< int, int> > , std::vector< int > , std::vector< std::vector< int > > , std::vector< size_t > , std::vector< double > , std::vector< std::vector<size_t> > >
BallMapperCppInterfacePython( const std::vector< std::vector<double> >& points , const std::vector<double>& values , double epsilon )
{
  bool dbg = true;
  if (dbg) std::cerr << "+++++++++BallMapperCppInterfacePython+++++++++\n";
	int number_of_points = points.size();

  if (dbg)
  {
  std::cerr << "Number of points : " << number_of_points << std::endl;
  std::cerr << "values.size() : " << values.size() << std::endl;
  }

	std::vector< std::vector<size_t> > coverage( number_of_points );
	std::vector< size_t > landmarks;
	landmarks.reserve( (size_t)(0.2*number_of_points) );

	//here we outsource computations of landmark points:
  if (dbg) std::cerr << "Entering compute_landmarks.\n";
	compute_landmarks< std::vector<double> >( points , coverage, landmarks , epsilon , number_of_points );
  if (dbg)
  {
  std::cerr << "landmarks.size() : " << landmarks.size() << std::endl;
    std::cerr << "Here are the landmarks: \n";
    for ( size_t i = 0 ; i != landmarks.size() ; ++i )
    {
      std::cerr << landmarks[i] << " , ";
    }
    std::cerr << "coverage.size() : " << coverage.size() << std::endl;
  }

	std::vector< std::vector< int > > points_covered_by_landmarks;
	std::vector< int > numer_of_covered_points( landmarks.size() , 2 );
  if (dbg) std::cerr << "Entering internal_procedure_fill_points_covered_by_landmarks.\n";
	internal_procedure_fill_points_covered_by_landmarks< std::vector< int > >( numer_of_covered_points , points_covered_by_landmarks ,  landmarks , coverage );

	//now let us deal with the coloring of the vertices:
	//change2
	std::vector< double > coloring;
  if (dbg) std::cerr << "Entering internal_procedure_fill_coloring.\n";
	internal_procedure_fill_coloring< std::vector<double> >( coloring , points_covered_by_landmarks , values );

	//Now let us create the graph and record the strength of each edge. To measure this, we will create the incidence matrix of the graph.
	//change3
	std::vector< std::vector< int > > graph_incidence_matrix;
	std::vector<int> from;
	std::vector<int> to;
	std::vector<int> strength_of_edges;
  if (dbg) std::cerr << "Entering internal_procedure_build_graph.\n";
	internal_procedure_build_graph< std::vector<int> >( graph_incidence_matrix, from, to, strength_of_edges, landmarks, coverage );

  // create the vector of edges
  std::vector<std::pair<int, int> > edges;
  if (dbg) std::cerr << "creating vector of edges.\n";
  for (int i = 0; i < from.size(); i++)
  {
      edges.push_back(std::make_pair(from[i], to[i]));
  }


	if (dbg) std::cerr << "+++++++++ THE END +++++++++\n\n";
	return std::make_tuple( numer_of_covered_points , edges , strength_of_edges , points_covered_by_landmarks , landmarks , coloring , coverage );

}//BallMapperCppInterfacePython




///////////////////////////////////////////////
// For symmetric data

std::vector< std::vector< std::pair<unsigned , double > > > dense_to_sparse_vector( const std::vector< std::vector<double> >& dense_points )
{
    std::vector< std::vector< std::pair<unsigned , double > > > sparse_points( dense_points.size() );

    for ( int i = 0 ; i != dense_points.size() ; ++i )
    {
        std::vector< std::pair<unsigned , double > > pt;
        pt.reserve( dense_points[i].size() );
        for ( size_t j = 0 ; j != dense_points[i].size() ; ++j )
        {
        	if ( dense_points[j][i] != 0 )
        	{
      		   pt.push_back( std::make_pair( j , dense_points[i][j] ) );
      		}
        }
        sparse_points[i] = pt;
    }
    return sparse_points;
}// dense_to_sparse_vector



inline double compute_distance_standard_points_sparse_points
  ( const std::vector< std::pair<unsigned,double> >& pt1 , const std::vector< std::pair<unsigned,double> >& pt2 , double p = 2 )
{
  double result = 0;
  size_t pt1_it = 0;
  size_t pt2_it = 0;

  bool dbg = false;
  if ( dbg )
  {
    #ifndef rcpp_code
        std::cerr << pt1.size() << " " << pt2.size() << std::endl;
    #endif

    std::cerr << "pt1: " ;
    for (size_t i=0; i < pt1.size(); i++) {
      std::cerr << "(" << pt1[i].first << " " <<  pt1[i].second << "), " ;
    }
    std::cerr << std::endl;

    std::cerr << "pt2: " ;
    for (size_t i=0; i < pt2.size(); i++) {
      std::cerr << "(" << pt2[i].first << " " <<  pt2[i].second << "), " ;
    }
    std::cerr << std::endl;
  }

  while ( ( pt1_it < pt1.size() ) && ( pt2_it < pt2.size() ) )
  {
  	if ( pt1[pt1_it].first < pt2[pt2_it].first )
  	{
	        //we move pt1:
	        result += pow( ( pt1[pt1_it].second ) , p );
		++pt1_it;
	}
	else
	{
		if ( pt1[pt1_it].first > pt2[pt2_it].first )
		{
			//we move pt1:
	        	result += pow( ( pt2[pt2_it].second ) , p );
			++pt2_it;
		}
		else
		{
			//in this case pt1[pt1_it].first == pt2[pt2_it]
			result += pow( ( pt1[pt1_it].second - pt2[pt2_it].second ) , p );
			++pt1_it;
			++pt2_it;
		}
	}
  }

  //if there is anything left, we count it here:
  while ( pt1_it < pt1.size() )
  {
	  result += pow( ( pt1[pt1_it].second ) , p );
	  ++pt1_it;
  }
  while ( pt2_it < pt2.size() )
  {
	result += pow( ( pt2[pt2_it].second ) , p );
	pt2_it++;
  }


  if ( dbg )
  {
    std::cerr << "distance: " << pow(result,1/p) << std::endl <<  std::endl;
  }

  return pow(result,1/p);
}//compute_distance_standard_points_sparse_points



void compute_landmarks_not_transposed_pts_group_action_sparse_points
                                         ( std::vector< std::vector<size_t> >& coverage ,
                                           std::vector< size_t > & landmarks ,
                                           const std::vector< std::vector< std::pair<unsigned,double> > >& points ,
                                           const std::vector< std::vector<int> >& orbit ,
                                           double epsilon
                                         )
{
   bool dbg = false;

   //here we outsource computations of landmark points:
  size_t current_point = 0;
  size_t current_landmark = 0;

  if ( dbg )
  {
     std::cerr << "orbit.size() : " << orbit.size() << std::endl;
  }

  while ( true )
  {
      while ( (current_point != points.size()) && (coverage[current_point].size() != 0) )
      {
          ++current_point;
      }

      if ( dbg ) std::cerr << "Current point : " << current_point << std::endl;

      if ( current_point == points.size() )break;
      for ( size_t i = 0 ; i != orbit[ current_point ].size() ; ++i )
      {
           if ( dbg ) std::cerr << "orbit["<<current_point<<"][" << i << "] : " << orbit[ current_point ][i] << std::endl;

           landmarks.push_back( orbit[ current_point ][i] );
           for ( size_t j = 0 ; j != points.size() ; ++j )
           {
              if ( compute_distance_standard_points_sparse_points( points[j] , points[ orbit[ current_point ][i]-1 ] ) <= epsilon )
              {
                  coverage[j].push_back( current_landmark+1);
              }
           }
           current_landmark++;
      }
      if ( dbg ) std::cerr << "Out of the internal while loop. \n";
  }
}


std::tuple< std::vector< int > , std::vector<std::pair< int, int> > , std::vector< int > , std::vector< std::vector< int > > , std::vector< size_t > , std::vector< double > , std::vector< std::vector<size_t> > >
SimplifiedBallMapperCppInterfaceGroupActionAndSparseRepresentationPython( const std::vector< std::vector<double> >& dense_points , const std::vector<double>& values ,
                                                                          double epsilon , const std::vector< std::vector<int> > orbit )
{
  bool dbg = true;
  if (dbg) std::cerr << "+++++++++SimplifiedBallMapperCppInterfaceGroupActionAndSparseRepresentationPython+++++++++\n\n";
  if ( dense_points.size() == 0 )
  {
    std::cerr << "No points in the BallMapperCpp procedure, the program will now terminate";
    throw "No points in the BallMapperCpp procedure, the program will now terminate";
  }

  std::vector< std::vector< std::pair<unsigned , double > > > points = dense_to_sparse_vector(dense_points);

  int number_of_points = points.size();
  if (dbg)
  {
	std::cerr << "Number of points : " << number_of_points << std::endl;
	std::cerr << "orbit.size() : " << orbit.size() << std::endl;
	std::cerr << "values.size() : " << values.size() << std::endl;
  }

  std::vector< std::vector<size_t> > coverage( number_of_points );
  std::vector< size_t > landmarks;
  landmarks.reserve( (size_t)(0.2*number_of_points) );


  if (dbg) std::cerr << "Entering compute_landmarks_not_transposed_pts_group_action.\n";
  compute_landmarks_not_transposed_pts_group_action_sparse_points( coverage , landmarks ,  points ,  orbit , epsilon );

  if (dbg)
  {
	std::cerr << "landmarks.size() : " << landmarks.size() << std::endl;
    std::cerr << "Here are the landmarks: \n";
    for ( size_t i = 0 ; i != landmarks.size() ; ++i )
    {
      std::cerr << landmarks[i] << " , ";
    }
    std::cerr << "coverage.size() : " << coverage.size() << std::endl;
  }


  std::vector< std::vector< int > > points_covered_by_landmarks;
	std::vector< int > numer_of_covered_points( landmarks.size() , 2 );
  if (dbg) std::cerr << "Entering internal_procedure_fill_points_covered_by_landmarks.\n";
	internal_procedure_fill_points_covered_by_landmarks< std::vector< int > >( numer_of_covered_points , points_covered_by_landmarks ,  landmarks , coverage );

  std::vector< double > coloring;
  if (dbg) std::cerr << "Entering internal_procedure_fill_coloring.\n";
  internal_procedure_fill_coloring< std::vector< double > >( coloring , points_covered_by_landmarks , values );


  std::vector< std::vector< int > > graph_incidence_matrix;
  std::vector< int > from;
  std::vector< int > to;
  std::vector< int > strength_of_edges;
  if (dbg) std::cerr << "Entering internal_procedure_build_graph.\n";
  internal_procedure_build_graph< std::vector< int > >( graph_incidence_matrix, from, to, strength_of_edges, landmarks, coverage );

  // create the vector of edges
  std::vector<std::pair<int, int> > edges;
  if (dbg) std::cerr << "creating list of edges.\n";
  for (int i = 0; i < from.size(); i++)
  {
      edges.push_back(std::make_pair(from[i], to[i]));
  }

  if (dbg) std::cerr << "+++++++++ THE END +++++++++\n\n";
  return std::make_tuple(numer_of_covered_points , edges, strength_of_edges, points_covered_by_landmarks, landmarks, coloring, coverage);
}//SimplifiedBallMapperCppInterfaceGroupActionAndSparseRepresentationPython
