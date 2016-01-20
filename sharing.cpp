//

// sharing.cpp
//
// Copyright (C) 2013 - 2015 jones@scss.tcd.ie
//
// This program is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software Foundation;
// either version 2 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software Foundation Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
//
// 19/11/12 first version
// 19/11/12 works with Win32 and x64
// 21/11/12 works with Character Set: Not Set, Unicode Character Set or Multi-Byte Character
// 21/11/12 output results so they can be easily pasted into a spreadsheet from console
// 24/12/12 increment using (0) non atomic increment (1) InterlockedIncrement64 (2) InterlockedCompareExchange
// 12/07/13 increment using (3) RTM (restricted transactional memory)
// 18/07/13 added performance counters
// 27/08/13 choice of 32 or 64 bit counters (32 bit can oveflow if run time longer than a couple of seconds)
// 28/08/13 extended struct Result
// 16/09/13 linux support (needs g++ 4.8 or later)
// 21/09/13 added getWallClockMS()
// 12/10/13 Visual Studio 2013 RC
// 12/10/13 added FALSESHARING
// 14/10/14 added USEPMS
//

//
// NB: hints for pasting from console window
// NB: Edit -> Select All followed by Edit -> Copy
// NB: paste into Excel using paste "Use Text Import Wizard" option and select "/" as the delimiter
//

#include "stdafx.h"                             // pre-compiled headers
#include <iostream>                             // cout
#include <iomanip>                              // setprecision
#include "helper.h"                             //


#include <random>	
#include <assert.h>

#include <fstream>
using namespace std;                            // cout

#define K           1024                        //
#define GB          (K*K*K)                     //
#define NOPS        10000                    //
#define NSECONDS    2                           // run each test for NSECONDS

#define COUNTER64                               // comment for 32 bit counter
												//#define FALSESHARING                          // allocate counters in same cache line
												//#define USEPMS                                // use PMS counters

#ifdef COUNTER64
#define VINT    UINT64                          //  64 bit counter
#else
#define VINT    UINT                            //  32 bit counter
#endif

#define ALIGNED_MALLOC(sz, align) _aligned_malloc(sz, align)

#ifdef FALSESHARING
#define GINDX(n)    (g+n)                       //
#else
#define GINDX(n)    (g+n*lineSz/sizeof(VINT))   //
#endif


#pragma region MCS_INIT
 /*
template <class T>
	class ALIGNEDMA
	{
	public:
		void* operator new(size_t); //override new
		void operator delete(void*); //override delete
	};

template<class T> 
void* ALIGNEDMA<T>::operator new(size_t sz) //aligned memory allocator
{
	sz = (sz + lineSz - 1) / lineSz * lineSz; //Make sz a multiple of lineSz
	return _aligned_malloc(sz, lineSz);	//Allocate on a lineSz boundary
}

template <class T>
void ALIGNEDMA<T>::operator delete(void* p)
{
	_aligned_free(p);
}

//MCS Lock
class QNode : public ALIGNEDMA<QNode>
{
public:
	volatile int waiting;
	volatile QNode *next;

	QNode()
	{
		this->waiting = 0;
		this->next = NULL;
	}
};
*/
#pragma endregion MCS_INIT


int bound = 0;

int upperBounds[] =
{
	16,	
	256,
	4096,
	65536,
	1048576
};


/**
	TREE_TYPE = 0   -> Standard BST
				1	-> BST w/ TTAS
				2	-> BST w/ HLE
				3	-> BST w/ RTM
*/
#define TREE_TYPE	0		

/**
	RTM_TYPE	= 0	-> Standard RTM (No Backoff/retries)
				= 1 -> RTM w/ backoff
*/

#define RTM_TYPE	0

UINT64 nabort;

#pragma region NUMBER_GENERATOR

/**
Random Number Generator
*/

class RandomNumberGenerator {
public:
	int generate(int upper);
	UINT64 generate64(int upper);
};

/**
	Generates random numbers using a uniform distributon model
	Upper limits for the assignment are
		16 [0..15]
		256
		4096
		65536
		1048576
*/
int RandomNumberGenerator::generate(int upper)
{
	//Set up generator
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<int> dis(0, upper);

	//Generate Random No & return
	return dis(gen);
}

/**
	Generates 64bit numbers, just testing it
*/
UINT64 RandomNumberGenerator::generate64(int upper)
{
	//Set up generator
	random_device rd;
	mt19937_64 gen(rd());
	uniform_int_distribution<UINT64> dis(0, upper);

	//I shoud probably do something about all of these instances
	//Are they automatically destroyed, or do I have to do some manual garbage collection?

	//Generate Random No & return
	return dis(gen);
}

RandomNumberGenerator rn;

#pragma endregion NUMBER_GENERATOR

#pragma region BST

//Binary Search Trees
class Node
{
public:
	UINT64 volatile key;
	Node* volatile left;
	Node* volatile right;

	Node()
	{
		key = 0;
		right = left = NULL;
	}
};

class BST
{
public:
	Node* volatile root;

	BST()
	{
		root = NULL;
		nodeCount = 0;
	}

	void initialize(int upper);
	//int contains(INT64 key);
	void add(Node *nm);
	Node* createNode(UINT64 value);
	void createAndInsert(UINT64 value);
	Node* volatile remove(UINT64 key);
	int size(Node* root);

	int getNodeCount()
	{
		return this->nodeCount;
	}

	void incNodeCount()
	{
		this->nodeCount++;
	}

	void decNodeCount()
	{
		this->nodeCount--;
	}

private:
	int nodeCount;
};

//int BST::contains(INT64 key) {}

void BST::add(Node *n)
{
	if (root == NULL)
	{
		root = n;
		return;
	}

	Node* volatile* pp = &root;
	Node* volatile p = root;
	while (p)
	{
		if (n->key < p->key)
		{
			pp = &p->left;
		}
		else if (n->key > p->key)
		{
			pp = &p->right;
		}
		else
		{
			return;
		}
		p = *pp;
	}
	*pp = n;
}

Node* BST::createNode(UINT64 value)
{
	Node* volatile n = new Node();
	n->key = value;
	return n;
}

void BST::createAndInsert(UINT64 value)
{
	//Creates a new node and inserts it into the BST
	Node* volatile n = new Node();
	n->key = value;
	this->add(n);
}

Node* volatile BST::remove(UINT64 key)
{
	Node* volatile *pp = &root;
	Node* volatile p = root;
	while (p)
	{
		if (key < p->key)
		{
			pp = &p->left;
		}
		else if (key > p->key)
		{
			pp = &p->right;
		}
		else {
			break;
		}
		p = *pp;
	}

	if (p == NULL)
		return NULL;

	if (p->left == NULL && p->right == NULL)
	{
		*pp = NULL; // NO children	
	}
	else if (p->left == NULL)
	{
		*pp = p->right; // ONE child
	}
	else if (p->right == NULL)
	{
		*pp = p->left; // ONE child
	}
	else
	{
   		Node* volatile r = p->right; // TWO children
		Node* volatile *ppr = &p->right; // find min key in right sub tree

		while (r->left) {
			ppr = &r->left;
			r = r->left;
		}
		p->key = r->key; // could move...
		p = r; // node instead
		*ppr = r->right;
	}
	return p; // return removed node

}

/**
	Returns the size of the BST
	Uses the standard recursive alogrithm for finding the size of a binary tree
*/

int BST::size(Node* root)
{
	if (root == NULL) 
	{
		return 0; 
	}
	else 
	{
		return 1 + size(root->left) + size(root->right);
	}
}

/**
	Pre-fills the tree with random 64-bit values
	Depends on the bounds
*/
void BST::initialize(int upper)
{
	int amount = upper / 2;

	for (int i = 0; i < amount; i++)
	{
		UINT64 value = rn.generate64(upper);
		createAndInsert(value);
	}
}

BST* volatile tree;

#pragma endregion BST

												//
												// OPTYP
												//
												// 0:inc
												// 1:InterlockedIncrement
												// 2:InterlockedCompareExchange
												// 3:RTM (restricted transactional memory)
												// 4:Bakery Lock
												// 5:TestAndSet Lock
												// 6:TestAndTestAndSet Lock
												// 7:MCS Lock
#define MAXTHREAD 8
int number[MAXTHREAD];
volatile int choosing[MAXTHREAD];
int pid = 0;

//DWORD tlsIndex = TlsAlloc();
//DWORD tlsIndex;
//QNode* mcsLock;	//MCS Lock
//volatile UINT64 lock = 0;

volatile long long lock = 0;		//Everything Else

#define MINBACKOFF 4				//RTM retry value
#define LOCKBACKOFF 16

#define OPTYP   0                    // set op type

#pragma region LOCK_MACROS
#if OPTYP == 0

#define OPSTR       "inc"
#define INC(g)      (*g)++;

#elif OPTYP == 1

#ifdef COUNTER64
#define OPSTR       "InterlockedIncrement64"
#define INC(g)      InterlockedIncrement64((volatile LONG64*) g)
#else
#define OPSTR       "InterlockedIncrement"
#define INC(g)      InterlockedIncrement(g)
#endif

#elif OPTYP == 2

#ifdef COUNTER64
#define OPSTR       "InterlockedCompareExchange64"
#define INC(g)      do {                                                                        \
                        x = *g;                                                                 \
                    } while (InterlockedCompareExchange64((volatile LONG64*) g, x+1, x) != x);
#else
#define OPSTR       "InterlockedCompareExchange"
#define INC(g)      do {                                                                        \
                        x = *g;                                                                 \
                    } while (InterlockedCompareExchange(g, x+1, x) != x);
#endif

#elif OPTYP == 3

#define OPSTR       "RTM (restricted transactional memory)"
#define INC(g)      {                                                                           \
                        UINT status = _xbegin();                                                \
                        if (status == _XBEGIN_STARTED) {                                        \
                            (*g)++;                                                             \
                            _xend();                                                            \
                        } else {                                                                \
                            nabort++;                                                           \
                            InterlockedIncrement64((volatile LONG64*)g);                        \
                        }                                                                       \
                    }
										
#endif
#pragma endregion LOCK_MACROS

#pragma region MCSLOCK
/*
DWORD tlsIndex = TlsAlloc();
void mcsAcquire(QNode **lock)
{
	volatile QNode *qn = (QNode*)TlsGetValue(tlsIndex);
	qn->next = NULL;
	volatile QNode *pred = (QNode*)InterlockedExchangePointer((PVOID*)lock, (PVOID)qn);
	if (pred == NULL)
	{
		return;
	}
	qn->waiting = 1;
	pred->next = qn;
	while (qn->waiting)
	{
		//MCS goes absolutely nuts without this Sleep operation
		Sleep(0);
	}
}

inline void mcsRelease(QNode **lock)
{
	volatile QNode *qn = (QNode*)TlsGetValue(tlsIndex); 
	volatile QNode *succ; 
	if (!(succ = qn->next)) 
	{ 
		if (InterlockedCompareExchangePointer((PVOID*)lock, NULL, (PVOID)qn) == qn)
			return; 
		do 
		{ 
			succ = qn->next; 
		} while (!succ); 
	}
	succ->waiting = 0;
}
*/
#pragma endregion MCSLOCK

#pragma region TTAS_LOCK

void ttas_acquire()
{
#ifdef WIN32
	do																			\
	{																			\
		while (lock == 1)														\
			_mm_pause(); 														\
	} while (InterlockedExchange64(&lock, 1));
#elif __linux__
	do																			\
	{																			\
		while (lock == 1)														\
			_mm_pause(); 														\
	} while (__sync_lock_test_and_set(&lock, 1));
#endif
}

void ttas_release()
{
	lock = 0;
}

#pragma endregion TTAS_LOCK

#pragma region HLE

void hle_acquire()
{
#ifdef WIN32
	while (_InterlockedExchange64_HLEAcquire(&lock, 1))							\
	{																			\
		do																		\
		{																		\
			_mm_pause();														\
		} while (lock == 1);													\
	}
#elif __linux__
	while (__atomic_exchange_n(&lock, 1, __ATOMIC_ACQUIRE | __ATOMIC_HLE_ACQUIRE))\
	{																			\
		do																		\
		{																		\
			_mm_pause();														\
		} while (lock == 1);													\
	}
#endif
}

void hle_release()
{
#ifdef WIN32
	_Store64_HLERelease(&lock, 0);

#elif __linux__
	__atomic_store_n(&lock, 0, __ATOMIC_RELEASE | __ATOMIC_HLE_RELEASE);
#endif
}

#pragma endregion HLE


#pragma region BST_OPERATIONS

Node* createNode(UINT64 value)
{
	Node* volatile n = new Node();
	n->key = value;
	return n;
}

void standardInsert(UINT64 value)
{			
	tree->createAndInsert(value);																						
}

void standardInsert(Node* n)
{
	tree->add(n);
}

void standardDelete(UINT64 value)
{
	tree->remove(value);
}

void ttasInsert(UINT64 value)
{
	//ttas_acquire();
	//tree->createAndInsert(value);
	//ttas_release();

	Node* n = createNode(value);

	/// NOTE TO LATER SELF - CREATE NODES IN WORKER INSTEAD
	ttas_acquire();
	tree->add(n);
	ttas_release();
}

void ttasInsert(Node* n)
{
	ttas_acquire();
	tree->add(n);
	ttas_release();
}

void ttasRemove(UINT64 value)
{
	ttas_acquire();
	tree->remove(value);
	ttas_release();
}

void hleInsert(UINT64 value)
{
	Node* n = createNode(value);

	hle_acquire();
	tree->add(n);
	hle_release();
}

void hleInsert(Node* n)
{
	hle_acquire();
	tree->add(n);
	hle_release();
}

void hleRemove(UINT64 value)
{
	hle_acquire();
	tree->remove(value);
	hle_release();
}

//Because my laptop doesn't support RTM (I wonder why HLE isn't flagged)
#if TREE_TYPE == 3
void rtmInsert(UINT64 value)
{
	Node* n = createNode(value);

	UINT status = _xbegin();
	if (status == _XBEGIN_STARTED)
	{
		tree->createAndInsert(value);
		_xend();
	}	
	else
	{
		nabort++;
		//On abort, attempt to acquire a standard lock and carry out the operation
		ttasInsert(value);
	}
}

void rtmInsert(Node* n)
{
	Node* n = createNode(value);

	UINT status = _xbegin();
	if (status == _XBEGIN_STARTED)
	{
		tree->add(n);
		_xend();
	}
	else
	{
		nabort++;
		//On abort, attempt to acquire a standard lock and carry out the operation
		ttasInsert(value);
	}
}

void rtmInsertAlt(UINT64 value)
{
	int backoff = MINBACKOFF;
	volatile UINT64 wait = 0;
	UINT status;
	while (1)
	{
		if (backoff < LOCKBACKOFF)
		{
			status = _xbegin();
		}
		else
		{
			ttas_acquire();
			status = -1;
		}

		if (status == _XBEGIN_STARTED)
		{
			if (backoff >= LOCKBACKOFF || lock == 0)
			{
				tree->createAndInsert(value);
				if (backoff < LOCKBACKOFF)
				{
					_xend();
				}
				else
				{
					lock = 0;
				}
				return;
			}
		}
		else
		{
			if (wait == 0)
			{
				//Can't get this damn intrinsic working on Linux
				//while (_rdrand64_step((UINT64*)&wait) == 0);

				//Using my generator instead
				wait = rn.generate64(MAXUINT64);
			}
			wait &= (1ULL << backoff) - 1;
			while (wait--);
			backoff++;
		}
	}
}

void rtmInsertAlt(Node* n)
{
	int backoff = MINBACKOFF;
	volatile UINT64 wait = 0;
	UINT status;
	while (1)
	{
		if (backoff < LOCKBACKOFF)
		{
			status = _xbegin();
		}
		else
		{
			ttas_acquire();
			status = -1;
		}

		if (status == _XBEGIN_STARTED)
		{
			if (backoff >= LOCKBACKOFF || lock == 0)
			{
				tree->add(n);
				if (backoff < LOCKBACKOFF)
				{
					_xend();
				}
				else
				{
					lock = 0;
				}
				return;
			}
		}
		else
		{
			if (wait == 0)
			{
				//Can't get this damn intrinsic working on Linux
				//while (_rdrand64_step((UINT64*)&wait) == 0);

				//Using my generator instead
				wait = rn.generate64(MAXUINT64);
			}
			wait &= (1ULL << backoff) - 1;
			while (wait--);
			backoff++;
		}
	}
}

void rtmRemoveAlt(UINT64 value)
{
	int backoff = MINBACKOFF;
	volatile UINT64 wait = 0;
	UINT status;
	while (1)
	{
		if (backoff < LOCKBACKOFF)
		{
			status = _xbegin();
		}
		else
		{
			ttas_acquire();
			status = -1;
		}

		if (status == _XBEGIN_STARTED)
		{
			if (backoff >= LOCKBACKOFF || lock == 0)
			{
				tree->remove(value);
				if (backoff < LOCKBACKOFF)
				{
					_xend();
				}
				else
				{
					lock = 0;
				}
				return;
			}
		}
		else
		{
			if (wait == 0)
			{
				//Can't get this damn intrinsic working on Linux
				//while (_rdrand64_step((UINT64*)&wait) == 0);

				//Using my generator instead
				wait = rn.generate64(MAXUINT64);
			}
			wait &= (1ULL << backoff) - 1;
			while (wait--);
			backoff++;
		}
	}
}

void rtmRemove(UINT64 value)
{
	UINT status = _xbegin();
	if (status == _XBEGIN_STARTED)
	{
		tree->remove(value);
		_xend();
	}
	else
	{
		nabort++;
		//On abort, attempt to acquire a standard lock and carry out the operation
		ttasRemove(value);
	}
}
#endif

#pragma endregion BST_OPERATIONS

UINT64 tstart;                                  // start of test in ms
int sharing;                                    // % sharing
int lineSz;                                     // cache line size
int maxThread;                                  // max # of threads

THREADH *threadH;                               // thread handles
UINT64 *ops;                                    // for ops per thread

#if OPTYP == 3 || TREE_TYPE == 3
UINT64 *aborts;                                 // for counting aborts
#endif

typedef struct {
	int sharing;                                // sharing
	int nt;                                     // # threads
	UINT64 rt;                                  // run time (ms)
	UINT64 ops;                                 // ops
	UINT64 incs;                                // should be equal ops
	UINT64 aborts;                              //
} Result;

Result *r;                                      // results
UINT indx;                                      // results index

volatile VINT *g;                               // NB: position of volatile

												//
												// test memory allocation [see lecture notes]
												//
ALIGN(64) UINT64 cnt0;
ALIGN(64) UINT64 cnt1;
ALIGN(64) UINT64 cnt2;
UINT64 cnt3;                                    // NB: in Debug mode allocated in cache line occupied by cnt0

												//
												// PMS
												//
#ifdef USEPMS

UINT64 *fixedCtr0;                              // fixed counter 0 counts
UINT64 *fixedCtr1;                              // fixed counter 1 counts
UINT64 *fixedCtr2;                              // fixed counter 2 counts
UINT64 *pmc0;                                   // performance counter 0 counts
UINT64 *pmc1;                                   // performance counter 1 counts
UINT64 *pmc2;                                   // performance counter 2 counts
UINT64 *pmc3;                                   // performance counter 2 counts

												//
												// zeroCounters
												//
void zeroCounters()
{
	for (UINT i = 0; i < ncpu; i++) {
		for (int j = 0; j < 4; j++) {
			if (j < 3)
				writeFIXED_CTR(i, j, 0);
			writePMC(i, j, 0);
		}
	}
}

//
// void setupCounters()
//
void setupCounters()
{
	if (!openPMS())
		quit();

	//
	// enable FIXED counters
	//
	for (UINT i = 0; i < ncpu; i++) {
		writeFIXED_CTR_CTRL(i, (FIXED_CTR_RING123 << 8) | (FIXED_CTR_RING123 << 4) | FIXED_CTR_RING123);
		writePERF_GLOBAL_CTRL(i, (0x07ULL << 32) | 0x0f);
	}

#if OPTYP == 3

	//
	// set up and enable general purpose counters
	//
	for (UINT i = 0; i < ncpu; i++) {
		writePERFEVTSEL(i, 0, PERFEVTSEL_EN | PERFEVTSEL_USR | RTM_RETIRED_START);
		writePERFEVTSEL(i, 1, PERFEVTSEL_EN | PERFEVTSEL_USR | RTM_RETIRED_COMMIT);
		writePERFEVTSEL(i, 2, PERFEVTSEL_IN_TXCP | PERFEVTSEL_IN_TX | PERFEVTSEL_EN | PERFEVTSEL_USR | CPU_CLK_UNHALTED_THREAD_P);  // NB: TXCP in PMC2 ONLY
		writePERFEVTSEL(i, 3, PERFEVTSEL_IN_TX | PERFEVTSEL_EN | PERFEVTSEL_USR | CPU_CLK_UNHALTED_THREAD_P);
	}

#endif

}

//
// void saveCounters()
//
void saveCounters()
{
	for (UINT i = 0; i < ncpu; i++) {
		fixedCtr0[indx*ncpu + i] = readFIXED_CTR(i, 0);
		fixedCtr1[indx*ncpu + i] = readFIXED_CTR(i, 1);
		fixedCtr2[indx*ncpu + i] = readFIXED_CTR(i, 2);
		pmc0[indx*ncpu + i] = readPMC(i, 0);
		pmc1[indx*ncpu + i] = readPMC(i, 1);
		pmc2[indx*ncpu + i] = readPMC(i, 2);
		pmc3[indx*ncpu + i] = readPMC(i, 3);
	}
}

#endif

//
// worker
//
WORKER worker(void *vthread)
{
	int thread = (int)((size_t)vthread);

	UINT64 n = 0;

	volatile VINT *gt = GINDX(thread);
	volatile VINT *gs = GINDX(maxThread);

#if OPTYP == 2
	VINT x;
#elif OPTYP == 3 || TREE_TYPE == 3
	UINT64 nabort = 0;
#elif OPTYP == 7
	QNode *qn = new QNode();
	TlsSetValue(tlsIndex, qn);
#endif

	pid = thread;
	runThreadOnCPU(thread % ncpu);

	while (1) {

		//
		// do some work
		//
		for (int i = 0; i < NOPS / 4; i++)
		{
			switch (TREE_TYPE)
			{
			case(0) :
			{
				//RANDOMLY Add/Delete Nodes
				int action = rn.generate(1); //Generate 0 or 1
				UINT64 value = rn.generate64(bound); //Generate value to be insterted/deleted
				Node* n = createNode(value);

				switch (action)
				{
				case(0) :
					//standardInsert(value);
					standardInsert(n);
					break;
				case(1) :
					standardDelete(value);
					break;
				}
				break;
			}


			case(1) :
			{
				//RANDOMLY Add/Delete Nodes
				int action = rn.generate(1); //Generate 0 or 1
				UINT64 value = rn.generate64(bound); //Generate value to be insterted/deleted
				Node* n = createNode(value);

				switch (action)
				{
				case(0) :
					//ttasInsert(value);
					ttasInsert(n);
					break;
				case(1) :
					ttasRemove(value);
					break;
				}
			}
			
			case(2) :
			{
				//RANDOMLY Add/Delete Nodes
				int action = rn.generate(1); //Generate 0 or 1
				UINT64 value = rn.generate64(bound); //Generate value to be insterted/deleted
				Node* n = createNode(value);

				switch (action)
				{
				case(0) :
					//hleInsert(value);
					hleInsert(n);
					break;
				case(1) :
					hleRemove(value);
					break;
				}
			}
			case(3) :
			{
			#if TREE_TYPE == 3 
				//RANDOMLY Add/Delete Nodes
				int action = rn.generate(1); //Generate 0 or 1
				UINT64 value = rn.generate64(bound); //Generate value to be insterted/deleted
				Node* n = createNode(value);

				switch (action)
				{
				case(0) :
					if (RTM_TYPE == 0)
						rtmInsert(value);
					else if (RTM_TYPE == 1)
						rtmInsertAlt(value);
					break;
				case(1) :
					if (RTM_TYPE == 0)
						rtmRemove(value);
					else if (RTM_TYPE == 1)
						rtmRemoveAlt(value);
					break;
				}
			#endif
			}
	
			}
		}
		n += NOPS;

		//
		// check if runtime exceeded
		//
		if ((getWallClockMS() - tstart) > NSECONDS * 1000)
			break;

	}

	ops[thread] = n;
#if OPTYP == 3 || TREE_TYPE == 3
	aborts[thread] = nabort;
#endif
	return 0;

}

/*
void createResultFile(Result *r)
{
	ofstream myfile;
	char* filename = "../results.csv";
	myfile.open(filename);
	//myfile << "Lock Implementation Results\n";
	//myfile << "sharing,nt,rt,ops,inc\n";
	myfile << "Sharing,No. Of Threads,Increments\n";
	for (UINT i = 0; i < indx; i++)
	{
		//myfile << r[i].sharing << "," << r[i].nt << "," << r[i].rt << "," << r[i].ops << "," << r[i].incs << endl;
		myfile << r[i].sharing << "," << r[i].nt << "," << r[i].incs << endl;
	}
	myfile.close();
}*/


//
// main
//
int main()
{
	ncpu = getNumberOfCPUs();   // number of logical CPUs
	maxThread = 2 * ncpu;       // max number of threads
	//maxThread = 2;
								//
								// get date
								//
	char dateAndTime[256];
	getDateAndTime(dateAndTime, sizeof(dateAndTime));

	//
	// console output
	//
	cout << getHostName() << " " << getOSName() << " sharing " << (is64bitExe() ? "(64" : "(32") << "bit EXE)";
#ifdef _DEBUG
	cout << " DEBUG";
#else
	cout << " RELEASE";
#endif
	cout << " [" << OPSTR << "]" << " NCPUS=" << ncpu << " RAM=" << (getPhysicalMemSz() + GB - 1) / GB << "GB " << dateAndTime << endl;
#ifdef COUNTER64
	cout << "COUNTER64";
#else
	cout << "COUNTER32";
#endif
#ifdef FALSESHARING
	cout << " FALSESHARING";
#endif
	cout << " NOPS=" << NOPS << " NSECONDS=" << NSECONDS << " OPTYP=" << OPTYP;
#ifdef USEPMS
	cout << " USEPMS";
#endif
	cout << endl;
	cout << "Intel" << (cpu64bit() ? "64" : "32") << " family " << cpuFamily() << " model " << cpuModel() << " stepping " << cpuStepping() << " " << cpuBrandString() << endl;
#ifdef USEPMS
	cout << "performance monitoring version " << pmversion() << ", " << nfixedCtr() << " x " << fixedCtrW() << "bit fixed counters, " << npmc() << " x " << pmcW() << "bit performance counters" << endl;
#endif

	//
	// get cache info
	//
	lineSz = getCacheLineSz();
	//lineSz *= 2;

	if ((&cnt3 >= &cnt0) && (&cnt3 < (&cnt0 + lineSz / sizeof(UINT64))))
		cout << "Warning: cnt3 shares cache line used by cnt0" << endl;
	if ((&cnt3 >= &cnt1) && (&cnt3 < (&cnt1 + lineSz / sizeof(UINT64))))
		cout << "Warning: cnt3 shares cache line used by cnt1" << endl;
	if ((&cnt3 >= &cnt2) && (&cnt3 < (&cnt2 + lineSz / sizeof(UINT64))))
		cout << "Warning: cnt2 shares cache line used by cnt1" << endl;

#if OPTYP == 3

	//
	// check if RTM supported
	//
	if (!rtmSupported()) {
		cout << "RTM (restricted transactional memory) NOT supported by this CPU" << endl;
		quit();
		return 1;
	}

#endif

	cout << endl;

	//
	// allocate global variable
	//
	// NB: each element in g is stored in a different cache line to stop false sharing
	//
	threadH = (THREADH*)ALIGNED_MALLOC(maxThread*sizeof(THREADH), lineSz);             // thread handles
	ops = (UINT64*)ALIGNED_MALLOC(maxThread*sizeof(UINT64), lineSz);                   // for ops per thread

#if OPTYP == 3 || TREE_TYPE == 3
	aborts = (UINT64*)ALIGNED_MALLOC(maxThread*sizeof(UINT64), lineSz);                // for counting aborts
#endif

#ifdef FALSESHARING
	g = (VINT*)ALIGNED_MALLOC((maxThread + 1)*sizeof(VINT), lineSz);                     // local and shared global variables
#else
	g = (VINT*)ALIGNED_MALLOC((maxThread + 1)*lineSz, lineSz);                         // local and shared global variables
#endif
	r = (Result*)ALIGNED_MALLOC(5 * maxThread*sizeof(Result), lineSz);                   // for results
	memset(r, 0, 5 * maxThread*sizeof(Result));                                           // zero

	indx = 0;

#ifdef USEPMS
	//
	// set up performance monitor counters
	//
	setupCounters();
#endif

	//
	// use thousands comma separator
	//
	setCommaLocale();

	//
	// header
	//
	cout << "bound";
	cout << setw(9) << "nt";
	cout << setw(10) << "rt";
	cout << setw(10) << "ops";
	cout << setw(12) << "rel";
	cout << setw(14) << "treeSize";
#if OPTYP == 3 || TREE_TYPE == 3
	cout << setw(8) << "commit";
#endif
	cout << endl;

	cout << "-----";              // sharing
	cout << setw(9) << "--";        // nt
	cout << setw(10) << "--";        // rt
	cout << setw(10) << "---";      // ops
	cout << setw(12) << "---";       // rel
	cout << setw(14) << "------";
#if OPTYP == 3 || TREE_TYPE == 3
	cout << setw(8) << "------";
#endif
	cout << endl;

	//
	// boost process priority
	// boost current thread priority to make sure all threads created before they start to run
	//
#ifdef WIN32
	//  SetPriorityClass(GetCurrentProcess(), ABOVE_NORMAL_PRIORITY_CLASS);
	//  SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_ABOVE_NORMAL);
#endif

	//
	// run tests
	//
	UINT64 ops1 = 1;

	//tree = new BST();

	for (bound = 0; bound < 5; bound++) {

		tree = new BST();
		//tree->initialize(upperBounds[bound]);

		for (int nt = 1; nt <= maxThread; nt *= 2, indx++) {

			tree->initialize(upperBounds[bound]);


			//
			//  zero shared memory
			//
			for (int thread = 0; thread < nt; thread++)
				*(GINDX(thread)) = 0;   // thread local
				*(GINDX(maxThread)) = 0;    // shared

			//
			// get start time
			//
			tstart = getWallClockMS();

			//
			// create worker threads
			//
			for (int thread = 0; thread < nt; thread++)
				createThread(&threadH[thread], worker, (void*)(size_t)thread);

			//
			// wait for ALL worker threads to finish
			//
			waitForThreadsToFinish(nt, threadH);
			UINT64 rt = getWallClockMS() - tstart;

#ifdef USEPMS
			saveCounters();             // save PMS counters
#endif

										//
										// save results and output summary to console
										//
			for (int thread = 0; thread < nt; thread++) {
				r[indx].ops += ops[thread];
				//r[indx].incs += *(GINDX(thread));
#if OPTYP == 3 || TREE_TYPE == 3
				r[indx].aborts += aborts[thread];
#endif
			}
			//r[indx].incs += *(GINDX(maxThread));
			if ((sharing == 0) && (nt == 1))
				ops1 = r[indx].ops;
			r[indx].sharing = sharing;
			r[indx].nt = nt;
			r[indx].rt = rt;

			cout << setw(8) << upperBounds[bound];
			cout << setw(8) << nt;
			cout << setw(10) << fixed << setprecision(2) << (double)rt / 1000;
			cout << setw(10) << r[indx].ops;
			cout << setw(12) << fixed << setprecision(2) << (double)r[indx].ops / ops1;
			cout << setw(14) << tree->size(tree->root);

#if OPTYP == 3 || TREE_TYPE == 3

			cout << setw(7) << fixed << setprecision(0) << 100.0 * (r[indx].ops - r[indx].aborts) / r[indx].ops << "%";

#endif

			//if (r[indx].ops != r[indx].incs)
			//	cout << " ERROR incs " << setw(3) << fixed << setprecision(0) << 100.0 * r[indx].incs / r[indx].ops << "% effective";

			cout << endl;

			//
			// delete thread handles
			//
			for (int thread = 0; thread < nt; thread++)
				closeThread(threadH[thread]);

		}

	}

	cout << endl;

	//
	// output results so they can easily be pasted into a spread sheet from console window
	//
	setLocale();
	cout << "bound/nt/rt/ops/incs";
#if OPTYP == 3 || TREE_TYPE == 3
	cout << "/aborts";
#endif
	cout << endl;
	for (UINT i = 0; i < indx; i++) {
		cout << r[i].sharing << "/" << r[i].nt << "/" << r[i].rt << "/" << r[i].ops << "/" << r[i].incs;
#if OPTYP == 3 || TREE_TYPE == 3
		cout << "/" << r[i].aborts;
#endif      
		cout << endl;
	}
	cout << endl;

	//createResultFile(r);
	quit();

	return 0;

}

// eof
