/*
 *	ECPS203 - Assignment 3
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

SC_MODULE(Mon)
{
	sc_in<int> a;
	sc_in<int> b;
	sc_in<int> f;
	sc_in<bool> clk;

	void display();

	SC_CTOR(Mon)
	{
		SC_THREAD(display)
		sensitive << clk.neg();
	}
};
