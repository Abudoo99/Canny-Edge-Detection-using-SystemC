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
#include "stim.h"
#include "mult.h"
#include "mon.h"

SC_MODULE(Top)
{
	sc_signal<int> asig, bsig, fsig;

	sc_clock testclk;

	Stim stim1;
	Mult uut;
	Mon mon1;

	SC_CTOR(Top):testclk("testclk", 10, SC_NS), stim1("stim1"), uut("uut"), mon1("mon1")
	{
		stim1.a.bind(asig);
		stim1.b.bind(bsig);
		stim1.clk.bind(testclk);	
		
		uut.a.bind(asig);
		uut.b.bind(bsig);
		uut.f.bind(fsig);

		mon1.a.bind(asig);
		mon1.b.bind(bsig);
		mon1.f.bind(fsig);
		mon1.clk.bind(testclk);
	}
};	
