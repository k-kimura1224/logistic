/* read_data.h */
/*********************************************/
/*	read data ( dim of matrix )					*/
/*	read data ( data )								*/
/*********************************************/

#include <stdio.h>
#include <stdlib.h>
#include "read_data.h"

void ReadFail(	int status
					)
{
	if( !status ){
		printf("Reading failed.\n");
		exit(0);
	}
}

void ReadDim(	const char	*filename,
					int			*m			/* size of data matrix */ 
					)
{
	FILE	*file;
	int	status;
	
	/* open file */
	file = fopen(filename, "r");
	if( file==NULL ){
		printf("Could not open file <%s>.\n", filename);
		exit(0);
	}

	/* read m */ 
	status = fscanf( file, "%d", m);
	ReadFail(status);

	fclose( file);
}

void ReadData(	const char	*filename,
					int			m,
					double		*data
					)
{
	int 	i;

	FILE	*file;
	int	buf;
	int	status;
	
	/* open file */
	file = fopen(filename, "r");
	if( file==NULL ){
		printf("Could not open file <%s>.\n", filename);
		exit(0);
	}

	/* skip m */
	status = fscanf( file, "%d", &buf);
	ReadFail(status);

	/* read data */
	for(i=0; i<(m*m); ++i){
		status = fscanf( file, "%lf", (data+i));
		ReadFail(status);
	}

	fclose( file);

}
