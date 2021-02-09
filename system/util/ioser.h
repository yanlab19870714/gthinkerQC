//########################################################################
//## Copyright 2018 Da Yan http://www.cs.uab.edu/yanda
//##
//## Licensed under the Apache License, Version 2.0 (the "License");
//## you may not use this file except in compliance with the License.
//## You may obtain a copy of the License at
//##
//## //http://www.apache.org/licenses/LICENSE-2.0
//##
//## Unless required by applicable law or agreed to in writing, software
//## distributed under the License is distributed on an "AS IS" BASIS,
//## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//## See the License for the specific language governing permissions and
//## limitations under the License.
//########################################################################

//Acknowledgements: the operator overloading is implemented based on pregel-mpi (https://code.google.com/p/pregel-mpi/) by Chuntao Hong.

#ifndef IOSER_H
#define IOSER_H

#include <stdio.h> //for FILE pointer
#include <string.h> //for memcpy
#include <sys/stat.h> //for file size
#include "serialization.h"

#include <vector>
#include <set>
#include <string>
#include <map>

#include "global.h"

//### newly added to make the file more IO-robust
#include <cassert>
#include <errno.h>
#include <iostream>

#define MAX_FWRITE_TRIALS 8 //how many times fwrite(.) can try to write "membuf"
//######

using namespace std;

#define STREAM_MEMBUF_SIZE 65536 //64k

//-------------------------------------

class ifbinstream {///////

private:
    char* membuf;
    size_t bufpos;
    size_t totpos;
    FILE * file;

public:

    ifbinstream()//empty
	{
		file = NULL;
		bufpos = 0;
		totpos = 0;
		membuf = new char[STREAM_MEMBUF_SIZE];
	}

    ifbinstream(const char* path)
    {
    	file = fopen(path, "wb");
        //### newly added to make the file more IO-robust
        if (file == NULL) {
            cout<<"Error opening file: "<<path<<endl;
            perror("Error printed by perror");
        }
        assert(file != NULL);
        //######
    	bufpos = 0;
    	totpos = 0;
    	membuf = new char[STREAM_MEMBUF_SIZE];
    }

    /* //old non-robust implementation
    inline void bufflush()
    {
    	fwrite(membuf, 1, bufpos, file);
    }
    */
    
    //### newly added to make the file more IO-robust
    inline void bufflush()
    {
        size_t bytes_written, curpos = 0, rest = bufpos;
        bool count = 0;
        do{
            if(count > 0)
            {
                cout<<_my_rank<<": should flush "<<rest<<" bytes, but flushes only "<<bytes_written<<" bytes !!!"<<endl;
                perror("Error printed by perror");
                assert(count == MAX_FWRITE_TRIALS);
                //------
                rest -= bytes_written;
                curpos += bytes_written;
                cout<<_my_rank<<": try to write the remaining "<<rest<<" bytes ..."<<endl;
            }
            bytes_written = fwrite(membuf + curpos, 1, rest, file);
            count++;
        }while(bytes_written != rest);
    }
    //######

    ~ifbinstream()
	{
        if(file == NULL)
        {
            delete[] membuf;
            return; //already closed
        }
        if(bufpos > 0) bufflush();
        int re = fclose(file);
        //### newly added to make the file more IO-robust
        if(re == EOF)
        {
            cout<<"Error closing file."<<endl;
            perror("Error printed by perror");
        }
        assert(re == 0);
        //######
        delete[] membuf;
	}

    inline size_t size()
    {
        return totpos;
    }

    void raw_byte(char c)
    {
    	if(bufpos == STREAM_MEMBUF_SIZE)
    	{
    		bufflush();
    		bufpos = 0;
    	}
    	membuf[bufpos] = c;
    	bufpos++;
        totpos++;
    }

    void raw_bytes(const void* ptr, size_t size)
    {
    	totpos += size;
    	size_t gap = STREAM_MEMBUF_SIZE - bufpos;
    	char * cptr = (char *)ptr;
    	if(gap < size)
    	{
    		memcpy(membuf + bufpos, cptr, gap);
    		bufpos = STREAM_MEMBUF_SIZE; //useful for correct exec of bufflush()
    		bufflush();
    		size -= gap;
    		cptr += gap;
    		while(size > STREAM_MEMBUF_SIZE)
    		{
    			memcpy(membuf, cptr, STREAM_MEMBUF_SIZE);
    			bufflush();
    			size -= STREAM_MEMBUF_SIZE;
    			cptr += STREAM_MEMBUF_SIZE;
    		}
    		memcpy(membuf, cptr, size);
    		bufpos = size;
    	}
    	else
    	{
    		memcpy(membuf + bufpos, ptr, size);
    		bufpos += size;
    	}
    }

    //below is for object reuse

    void close() //also for flushing
    {
    	if(file == NULL) return; //already closed
    	if(bufpos > 0) bufflush();
        //### newly added to make the file more IO-robust
        int re = fclose(file);
        if(re == EOF)
        {
            cout<<"Error closing file."<<endl;
            perror("Error printed by perror");
        }
        assert(re == 0);
        //######
    	file = NULL; //set status to closed
    }

    void open(const char* path) //it does not check whether you closed previous file
    {
    	file = fopen(path,"wb");
        //### newly added to make the file more IO-robust
        if (file == NULL) {
            cout<<"Error opening file: "<<path<<endl;
            perror("Error printed by perror");
        }
        //######
		bufpos = 0;
		totpos = 0;
    }

    bool is_open()
    {
    	return file != NULL;
    }

};

//make sure mm only contains one object (mm should be cleared before serializing an object)
ifbinstream & operator<<(ifbinstream & m, ibinstream mm)
{
    m.raw_bytes(mm.get_buf(), mm.size());
    return m;
}

ifbinstream & operator<<(ifbinstream & m, size_t i)
{
    m.raw_bytes(&i, sizeof(size_t));
    return m;
}

ifbinstream & operator<<(ifbinstream & m, bool i)
{
    m.raw_bytes(&i, sizeof(bool));
    return m;
}

ifbinstream & operator<<(ifbinstream & m, int i)
{
    m.raw_bytes(&i, sizeof(int));
    return m;
}

ifbinstream & operator<<(ifbinstream & m, long long i)
{
    m.raw_bytes(&i, sizeof(long long));
    return m;
}

ifbinstream & operator<<(ifbinstream & m, double i)
{
    m.raw_bytes(&i, sizeof(double));
    return m;
}

ifbinstream & operator<<(ifbinstream & m, char c)
{
    m.raw_byte(c);
    return m;
}

template <class T>
ifbinstream & operator<<(ifbinstream & m, const T* p)
{
    return m << *p;
}

template <class T>
ifbinstream & operator<<(ifbinstream & m, const vector<T>& v)
{
    m << v.size();
    for (typename vector<T>::const_iterator it = v.begin(); it != v.end(); ++it) {
        m << *it;
    }
    return m;
}

template <>
ifbinstream & operator<<(ifbinstream & m, const vector<int> & v)
{
    m << v.size();
    m.raw_bytes(&v[0], v.size() * sizeof(int));
    return m;
}

template <>
ifbinstream & operator<<(ifbinstream & m, const vector<double> & v)
{
    m << v.size();
    m.raw_bytes(&v[0], v.size() * sizeof(double));
    return m;
}

template <class T>
ifbinstream & operator<<(ifbinstream & m, const set<T> & v)
{
    m << v.size();
    for(typename set<T>::const_iterator it = v.begin(); it != v.end(); ++it) {
        m << *it;
    }
    return m;
}

ifbinstream & operator<<(ifbinstream & m, const string & str)
{
    m << str.length();
    m.raw_bytes(str.c_str(), str.length());
    return m;
}

template <class KeyT, class ValT>
ifbinstream & operator<<(ifbinstream & m, const map<KeyT, ValT> & v)
{
    m << v.size();
    for (typename map<KeyT, ValT>::const_iterator it = v.begin(); it != v.end(); ++it) {
        m << it->first;
        m << it->second;
    }
    return m;
}

template <class KeyT, class ValT>
ifbinstream & operator<<(ifbinstream & m, const hash_map<KeyT, ValT> & v)
{
    m << v.size();
    for (typename hash_map<KeyT, ValT>::const_iterator it = v.begin(); it != v.end(); ++it) {
        m << it->first;
        m << it->second;
    }
    return m;
}

template <class T>
ifbinstream & operator<<(ifbinstream & m, const hash_set<T> & v)
{
    m << v.size();
    for (typename hash_set<T>::const_iterator it = v.begin(); it != v.end(); ++it) {
        m << *it;
    }
    return m;
}

//-------------------------------------

class ofbinstream {

private:
	char* membuf;
	size_t bufpos;
	size_t bufsize; //membuf may not be full (e.g. last batch)
	size_t totpos;
	size_t filesize;
	FILE * file;

public:
	inline void fill()
	{
		bufsize = fread(membuf, 1, STREAM_MEMBUF_SIZE, file);
		bufpos = 0;
	}

	ofbinstream()
	{
		membuf = new char[STREAM_MEMBUF_SIZE];
		file = NULL; //set status to closed
	}

	ofbinstream(const char* path)
	{
		membuf = new char[STREAM_MEMBUF_SIZE];
		file = fopen(path, "rb");
        //### newly added to make the file more IO-robust
        if (file == NULL) {
            cout<<"Error opening file: "<<path<<endl;
            perror("Error printed by perror");
        }
        //######
		//get file size
		filesize = -1;
		struct stat statbuff;
		if(stat(path, &statbuff) == 0) filesize = statbuff.st_size;
		//get first batch
		fill();
		totpos = 0;
	}

	bool open(const char* path) //return whether the file exists
	{
		file = fopen(path, "rb");
        //### newly added to make the file more IO-robust
        if (file == NULL) {
            cout<<"Error opening file: "<<path<<endl;
            perror("Error printed by perror");
        }
        //######
		if(file == NULL) return false;
		//get file size
		filesize = -1;
		struct stat statbuff;
		if(stat(path, &statbuff) == 0) filesize = statbuff.st_size;
		//get first batch
		fill();
		totpos = 0;
		return true;
	}

	inline size_t size()
	{
		return filesize;
	}

	inline bool eof()
	{
		return totpos >= filesize;
	}

    ~ofbinstream()
    {
    	delete[] membuf;
    	if(file == NULL) return; //already closed
		fclose(file);
    }

    char raw_byte()
    {
    	totpos++;
    	if(bufpos == bufsize) fill();
        return membuf[bufpos++];
    }

    void* raw_bytes(size_t n_bytes)
    {
    	totpos += n_bytes;
    	size_t gap = bufsize - bufpos;
    	if(gap >= n_bytes)
    	{
    		char* ret = membuf + bufpos;
    		bufpos += n_bytes;
			return ret;
    	}
    	else
    	{
    		//copy the last gap-batch to head of membuf
    		//!!! require that STREAM_MEMBUF_SIZE >= n_bytes !!!
    		memcpy(membuf, membuf + bufpos, gap);
    		//gap-shifted refill
    		bufsize = gap + fread(membuf + gap, 1, STREAM_MEMBUF_SIZE - gap, file);
    		bufpos = n_bytes;
    		return membuf;
    	}
    }

    void close()
	{
    	if(file == NULL) return; //already closed
        int re = fclose(file);
		//### newly added to make the file more IO-robust
        if(re == EOF)
        {
            cout<<"Error closing file."<<endl;
            perror("Error printed by perror");
        }
        assert(re == 0);
        //######
		file = NULL; //set status to closed
	}

    //=============== add skip function ===============
    void skip(size_t num_bytes)
    {
    	totpos += num_bytes;
    	if(totpos >= filesize) return; //eof
    	bufpos += num_bytes; //done if bufpos < bufsize
    	if(bufpos >= bufsize)
    	{
    		fseek(file, bufpos - bufsize, SEEK_CUR);
    		fill();
    	}
    }
};

ofbinstream & operator>>(ofbinstream & m, size_t & i)
{
    i = *(size_t*)m.raw_bytes(sizeof(size_t));
    return m;
}

ofbinstream & operator>>(ofbinstream & m, bool & i)
{
    i = *(bool*)m.raw_bytes(sizeof(bool));
    return m;
}

ofbinstream & operator>>(ofbinstream & m, int & i)
{
    i = *(int*)m.raw_bytes(sizeof(int));
    return m;
}

ofbinstream & operator>>(ofbinstream & m, double & i)
{
    i = *(double*)m.raw_bytes(sizeof(double));
    return m;
}

ofbinstream & operator>>(ofbinstream & m, long long & i)
{
    i = *(long long*)m.raw_bytes(sizeof(long long));
    return m;
}

ofbinstream & operator>>(ofbinstream & m, char & c)
{
    c = m.raw_byte();
    return m;
}

template <class T>
ofbinstream & operator>>(ofbinstream & m, T* & p)
{
    p = new T;
    return m >> (*p);
}

template <class T>
ofbinstream & operator>>(ofbinstream & m, vector<T> & v)
{
    size_t size;
    m >> size;
    v.resize(size);
    for (typename vector<T>::iterator it = v.begin(); it != v.end(); ++it) {
        m >> *it;
    }
    return m;
}

template <>
ofbinstream & operator>>(ofbinstream & m, vector<int> & v)
{
    size_t size;
    m >> size;
    vector<int>::iterator it = v.begin();
    size_t len = STREAM_MEMBUF_SIZE / 2 / sizeof(int);
    size_t bytes = len * sizeof(int);
    while(size > len)
    {
        int* data = (int*)m.raw_bytes(bytes);
        v.insert(it, data, data + len);
        it = v.end();
        size -= len;
    }
    int* data = (int*)m.raw_bytes(sizeof(int) * size);
    v.insert(it, data, data + size);
    return m;
}

template <>
ofbinstream & operator>>(ofbinstream & m, vector<double> & v)
{
    size_t size;
    m >> size;
    vector<double>::iterator it = v.begin();
    size_t len = STREAM_MEMBUF_SIZE / 2 / sizeof(double);
    size_t bytes = len * sizeof(double);
    while(size > len)
    {
        double* data = (double*)m.raw_bytes(bytes);
        v.insert(it, data, data + len);
        it = v.end();
        size -= len;
    }
    double* data = (double*)m.raw_bytes(sizeof(double) * size);
    v.insert(it, data, data + size);
    return m;
}

template <class T>
ofbinstream & operator>>(ofbinstream & m, set<T> & v)
{
    size_t size;
    m >> size;
    for (size_t i = 0; i < size; i++) {
        T tmp;
        m >> tmp;
        v.insert(v.end(), tmp);
    }
    return m;
}

ofbinstream & operator>>(ofbinstream & m, string & str)
{
    size_t length;
    m >> length;
    str.clear();

    size_t HALF_BUF = STREAM_MEMBUF_SIZE/2;
    while(length > HALF_BUF)
    {
        char* data = (char*)m.raw_bytes(HALF_BUF); //raw_bytes cannot accept input > STREAM_MEMBUF_SIZE
        str.append(data, HALF_BUF);
        length -= HALF_BUF;
    }
    char* data = (char*)m.raw_bytes(length);
    str.append(data, length);

    return m;
}

template <class KeyT, class ValT>
ofbinstream & operator>>(ofbinstream & m, map<KeyT, ValT> & v)
{
    size_t size;
    m >> size;
    for (size_t i = 0; i < size; i++) {
        KeyT key;
        m >> key;
        m >> v[key];
    }
    return m;
}

template <class KeyT, class ValT>
ofbinstream & operator>>(ofbinstream & m, hash_map<KeyT, ValT> & v)
{
    size_t size;
    m >> size;
    for (size_t i = 0; i < size; i++) {
        KeyT key;
        m >> key;
        m >> v[key];
    }
    return m;
}

template <class T>
ofbinstream & operator>>(ofbinstream & m, hash_set<T> & v)
{
    size_t size;
    m >> size;
    for (size_t i = 0; i < size; i++) {
        T key;
        m >> key;
        v.insert(key);
    }
    return m;
}

#endif
