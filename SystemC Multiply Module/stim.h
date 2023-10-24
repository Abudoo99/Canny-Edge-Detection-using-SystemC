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

SC_MODULE(Stim)
{
	sc_in<bool> clk;
	sc_out<int> a;
	sc_out<int>	b;	

	void stimulus();

	SC_CTOR(Stim)
	{
		SC_THREAD(stimulus);
		sensitive << clk.pos();
	}
};