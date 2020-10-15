/*----------------------------------
File     : data.cpp
Contents : data set management
----------------------------------*/
#include <string.h>

#include "data.h"

Data::Data(char *filename)
{
	in = fopen(filename,"rt");
	if(in==NULL)
	{
		printf("Error: cannot open file %s for read\n", filename);
		exit(-1);
	}

	mptransaction = new Transaction;
	mptransaction->t = new int[INIT_TRANS_LEN]; // transaction container; one transaction is the adj-list of one node
	mptransaction->length = 0;
	micapacity = INIT_TRANS_LEN;

}

Data::~Data()
{
	if(in)
		fclose(in);
	if(mptransaction!=NULL)
	{
		delete []mptransaction->t;
		delete mptransaction;
	}
}

int Data::isOpen()
{
	if(in) 
		return 1;
	else 
		return 0;
}

Transaction *Data::getNextTransaction()
{
	char c;

	mptransaction->length = 0;

	// read list of items
	do {
		int item=0, pos=0;
		c = getc(in);
		while((c >= '0') && (c <= '9')) 
		{
			item *=10;
			item += int(c)-int('0');
			c = getc(in);
			pos++;
		}
		if(pos)
		{
			if(mptransaction->length < micapacity)
			{
				mptransaction->t[mptransaction->length] = item;
				mptransaction->length++;
			}
			else 
			{
				int *ptrans;
				micapacity = 2*micapacity; // like std::vector, capacity reached so double the capacity
				ptrans = new int[micapacity];
				memcpy(ptrans, mptransaction->t, sizeof(int)*mptransaction->length);
				delete []mptransaction->t;
				mptransaction->t = ptrans;
				mptransaction->t[mptransaction->length] = item;
				mptransaction->length++;
			}
		}
	}while(c != '\n' && !feof(in)); // one line per transaction

	if(feof(in) && mptransaction->length==0)
	{
		rewind(in);
		return 0;
	}
	else 
		return mptransaction;
	
}


