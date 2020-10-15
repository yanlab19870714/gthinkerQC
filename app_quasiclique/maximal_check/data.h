#pragma once 

#pragma warning (disable:4996)

#include <stdio.h>
#include <stdlib.h>

#define INIT_TRANS_LEN  200

struct Transaction
{
	int length;
	int *t;
};

class Data
{
public:
	Transaction* mptransaction;
	int micapacity;

	
	Data(char *filename);
	~Data();

	int isOpen();
	Transaction *getNextTransaction();
  
private:
  
	FILE *in;
};

