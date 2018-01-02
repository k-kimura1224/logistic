/* read_data.h */
/*********************************************/
/*	read data ( dim of matrix )					*/
/*	read data ( data )								*/
/*********************************************/
#ifndef READ_DATA_H
#define READ_DATA_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

extern
void ReadFail(	int status
					);

extern
void ReadDim(	const char	*filename,
					int			*m			/* size of data matrix */ 
					);

extern
void ReadData(	const char	*filename,
					int			m,
					double		*data
					);

#ifdef __cplusplus
}
#endif

#endif
