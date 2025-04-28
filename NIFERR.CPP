// Filename: message.c
// Function: NIF system message
//

#include <iostream.h>
#include "nif.h"
#include "message.h"


void nif_msg(char **m, int n)
{
	if (INTERACT)
	{
		cout<<m[n]<<endl;
	}
}


int report_err(int err)
{
	char e1[5]={'P','T','V','H','S'};

	if (err<=5)
	{
		if (INTERACT) 
		{
			cout<<"error"<<err<<",  "<<e1[err-1]<<" overflow!\n";
			return 1;
		}
		else return (err);
	}
	else return 0;
}
