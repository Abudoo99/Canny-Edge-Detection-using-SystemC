/*
 * ECPS203 - Assignment 3
 *
 * Introduction to SystemC language and simulation
 *
 * Source code adapted from "SystemC Training: The Definitive Guide to SystemC"
 * by Doulos Ltd., 2015
 *
 * Author: Abyukth Kumar
 *
*/

#include "systemc.h"

SC_MODULE(Mult)
{
	sc_in<int> a;
	sc_in<int> b;
	sc_out<int> f;

	void multiply();

	SC_CTOR(Mult)
	{
		SC_METHOD(multiply);
			sensitive << a << b;
	}
};
