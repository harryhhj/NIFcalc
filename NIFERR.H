// Filename: message.h
// Function: system messages
//

static char msgNIF[6][20]={
	"Out of NIF range!\n",
	"P out of range!\n",
	"T out of range!\n",
	"V out of range!\n",
	"H out of range!\n",
	"S out of range!\n",
	};

void nif_msg(char **m, int n);
