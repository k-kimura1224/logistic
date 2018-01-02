#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <assert.h>

#include "svpsolver.h"
#include "node.h"

using namespace std;

#define debug	0


int	SVPsolver::select_node(
	int index,
	int disp
	)
{
	int	rslt = -1;
	switch ( NODESELECTION ) {
		case 0:
		{
			double	min;
			list<NODE>::iterator it = NodeList.begin();
			//min = NodeList[0].get_lowerbound();
			min = it->get_lowerbound();
			rslt = 0;
			for(int i=1; i<listsize; i++){
				++it;
				if( min > it->get_lowerbound() ){
					min = it->get_lowerbound();
					rslt = i;
				}
			}
			break;
		}
		case 1:
		{
			rslt = 0;
			break;
		}
		case 2:
		{
			
			if( (index-1) % 10000 == 0 && disp < index ){
				clock_t start = clock();
				NodeList.sort();
				clock_t end = clock();
				cout << "Sorting Time: ";
				cout << (double)(end-start)/CLOCKS_PER_SEC;
				cout << "s" << endl;
			}
			rslt = 0;
			break;
		}
		default:
		{
			cout << "error: select_node" << endl;
			exit(-1);
			break;
		}
	}
	
	assert( rslt >= 0 );
	assert( rslt < listsize );

	return rslt;
}
